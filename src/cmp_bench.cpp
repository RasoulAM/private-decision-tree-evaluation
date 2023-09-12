#include <omp.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <vector>
#include "seal/seal.h"
#include "utils.h"

using namespace seal;
using namespace std;

enum ComparisonType {
    RANGE_COVER=0,
    FOLKLORE=1,
    XCMP=2
};

inline string uint64_to_hex_string(uint64_t value) {
    return seal::util::uint_to_hex_string(&value, size_t(1));
}


 // computes range cover from [a, 2^n - 1]
void compRangeLargest(uint64_t a, uint64_t largest, uint64_t* range, int n)
{
    for (int i = 0; i < n + 1; i++) {
        range[i] = largest + 1;
    }
    // compute height of last common node in the path

    vector<int> exponents;
    uint64_t curr_largest = largest;

    // special case for 0, due to a being an unsigned int
    if (a == 0) {
        range[n] = 0;
        return;
    }

    while (curr_largest >= a) {
        int j = 0;

        while (curr_largest >= pow(2, j) && curr_largest - pow(2, j) + 1 >= a) {
            j++;
        }

        --j;

        uint64_t value = (largest + 1) / pow(2, j);

        for (auto& k : exponents) {
            value -= pow(2, k - j);
        }

        range[j] = value - 1;

        exponents.push_back(j);
        uint64_t powers_sum = 0;
        for (auto& q : exponents){
            uint64_t added = pow(2, q);
            // cout << powers_sum << " " << added << endl;
            powers_sum += added;
            // cout << powers_sum << endl;
        }
        curr_largest = largest - powers_sum;
    }
}

// Function for extended Euclidean Algorithm
int extended_euclidean(int a, int b, int *x, int *y) {
    if(a == 0) {
        *x = 0;
        *y = 1;
        return b;
    }

    int x1 = 0, y1 = 0;
    int gcd = extended_euclidean(b % a, a, &x1, &y1);

    *x = y1 - (b / a) * x1;
    *y = x1;

    return gcd;
}

int mod_inverse(int a, int m) {
    int x = 0, y = 0;
    int g = extended_euclidean(a, m, &x, &y);
    assert(g == 1);
    return (x % m + m) % m;
}

class ComparisonBenchmark{

public:
    uint64_t log_poly_mod_degree;
    uint64_t poly_modulus_degree;
    uint64_t row_count;
    int n;
    int hamming_weight;
    int code_length=0;
    uint64_t num_slots_per_element;
    uint64_t num_cmps; 

    ComparisonType comparison;

    EncryptionParameters parms;
    SEALContext* context;
    PublicKey pk;
    Encryptor *encryptor;
    Decryptor *decryptor;
    Evaluator *evaluator;
    BatchEncoder* batch_encoder;
    RelinKeys* rlk_server;

    GaloisKeys* gal_keys_server;
    uint64_t plain_modulus;

    uint64_t inv;

    uint64_t plain_op;
    uint64_t encrypted_op;

    bool correct;

    ComparisonBenchmark(){}

    void init(ComparisonType cmp, int nn, int hw){
        log_poly_mod_degree=14;
        comparison=cmp;
        n = nn;
        hamming_weight = hw;

        uint64_t prime_bitlength=20;
        switch (comparison) {
            case RANGE_COVER:
                if (2<=hamming_weight && hamming_weight<=4) this->log_poly_mod_degree=13;
                if (hamming_weight==1) this->log_poly_mod_degree=12;
                break;
            case FOLKLORE:
                log_poly_mod_degree=14;
                break;
            case XCMP:
                log_poly_mod_degree=13;
                prime_bitlength=18;
                if (n<=12){
                    log_poly_mod_degree=12;                    
                    prime_bitlength=16;
                }
        }

        poly_modulus_degree=1<<log_poly_mod_degree;
        row_count=poly_modulus_degree/2;
        num_slots_per_element=2*n+1;

        parms = EncryptionParameters(scheme_type::bfv);
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, prime_bitlength));
        parms.set_coeff_modulus(get_coeff_modulus(poly_modulus_degree));
        context = new SEALContext(parms);
        
        
        KeyGenerator keygen(*context);
        keygen.create_public_key(pk);
        SecretKey sk = keygen.secret_key();
        encryptor = new Encryptor(*context, pk);
        decryptor = new Decryptor(*context, sk);
        evaluator = new Evaluator(*context);
        batch_encoder = new BatchEncoder(*context);
        rlk_server = new RelinKeys();
        keygen.create_relin_keys(*rlk_server);

        gal_keys_server = new GaloisKeys();

        if (comparison==XCMP){
            vector<uint32_t> elts;
            for(size_t i = 2; i <= poly_modulus_degree; i *= 2) {
                elts.push_back(i + 1);
            }
            keygen.create_galois_keys(elts, *gal_keys_server);
        } else {
            keygen.create_galois_keys(*gal_keys_server);
        }

        plain_modulus = parms.plain_modulus().value();
        inv = mod_inverse(poly_modulus_degree, plain_modulus);

        srand(time(NULL));
        plain_op = rand() % ((uint64_t)1 << n);
        encrypted_op = rand() % ((uint64_t)1 << n);
    }

    void range_cover_compare(){

        uint64_t num_cmps_per_row = row_count/num_slots_per_element;
        this->num_cmps = 2*num_cmps_per_row;
        code_length = find_code_length(n, hamming_weight);
        // cout << "code_length: " << code_length << endl;

        vector<uint64_t> attribute_vector;
        for (int i = 0; i < num_cmps_per_row; i++) {
            attribute_vector.push_back(encrypted_op);
        }

        vector<uint64_t> server_values;
        for (int i = 0; i < num_cmps_per_row; i++) {
            server_values.push_back(plain_op);
        }

        // Ciphertext ct;
        // encryptor->encrypt(Plaintext("1"), ct);
        // vector<Ciphertext> client_input(code_length, ct);
        
        vector<Ciphertext> client_input(code_length);
        vector<vector<vector<uint64_t>>> codes;

        uint64_t largest = pow(2, n) - 1;
            uint64_t range[n + 1];
        for (int i = 0; i < num_cmps_per_row; i++) {
            compRangeLargest(attribute_vector[i], largest, range, n);

            std::vector<std::vector<uint64_t>> rangeCode;

            rangeCode.reserve(n + 1);

            uint64_t m = code_length;

            for (int i = n; i >= 0; --i) {
                if (range[i] == largest + 1) {
                    rangeCode.insert(rangeCode.begin(), std::vector<uint64_t>(m, 0));
                } else {
                    rangeCode.insert(rangeCode.begin(),
                        get_OP_CHaW_encoding(range[i], m, hamming_weight, true));
                }
            }
            codes.push_back(rangeCode);

        }

        /////
        // Constructing plaintexts/ciphertexts from the bits of the range codes
        for (int i = 0; i < code_length; ++i) {
            vector<uint64_t> pod_matrix(poly_modulus_degree, 0ULL);
            for (int k = 0; k < codes.size(); k++) {
                for (int j = 0; j < n + 1; ++j) {
                    if (codes[k][j].size() <= i) {
                        continue;
                    } else {
                        pod_matrix[k * num_slots_per_element + j] = codes[k][j][i];
                        pod_matrix[k * num_slots_per_element + j + row_count] = codes[k][j][i];
                    }
                }
            }
            Plaintext plain_matrix;
            Ciphertext ct;
            batch_encoder->encode(pod_matrix, plain_matrix);
            encryptor->encrypt(plain_matrix, ct);
            client_input[i] = ct;
        }

        /////

        vector<Plaintext> cts_range_bits_pathcode_vec;

        // inner is the array from previous code
        vector<vector<uint64_t>> paths(num_cmps_per_row, vector<uint64_t>());

        for (int i=0; i < num_cmps_per_row; i++) {
            uint64_t cur = server_values[i];
            for (int j = 0; j < n; j++) {
                paths[i].push_back(cur);
                cur /= 2;
            }
            paths[i].push_back(0);
        }

        uint64_t m = code_length;

        // middle is for array in old code
        std::vector<std::vector<std::vector<uint64_t>>> pathcodes;

        for (int j = 0; j < num_cmps_per_row; ++j) {
            vector<vector<uint64_t>> pathcode;
            pathcode.reserve(n + 1);
            for (int i = n; i >= 0; --i) {
                pathcode.insert(pathcode.begin(),
                    get_OP_CHaW_encoding(paths[j][i], m, hamming_weight, true));
            }
            pathcodes.push_back(pathcode);
        }

        // std::vector<Plaintext> cts_range_bits_pathcode;

        Plaintext plain_matrix;
        // Ciphertext ct;
        for (int i = 0; i < m; ++i) {
            vector<uint64_t> pod_matrix(poly_modulus_degree, 1ULL);
            for (int k = 0; k < num_cmps_per_row; k++) {
                for (int j = 0; j < n + 1; ++j) {
                    if (pathcodes[k][j].size() <= i) {
                        continue;
                    } else {
                        pod_matrix[k * num_slots_per_element + j] = pathcodes[k][j][i];
                        pod_matrix[k * num_slots_per_element + j + row_count] = pathcodes[k][j][i];
                    }
                }
            }
            batch_encoder->encode(pod_matrix, plain_matrix);
            cts_range_bits_pathcode_vec.push_back(plain_matrix);
        }

        // Constant-weight Comparison
        ////////////////

            vector<Ciphertext> operands;

            #pragma omp parallel for
            for (int i = 0; i < code_length; i++) {
                if (!cts_range_bits_pathcode_vec[i].is_zero()){
                    Ciphertext temp;
                    evaluator->multiply_plain(client_input[i], cts_range_bits_pathcode_vec[i], temp);
                    #pragma omp critical
                    {
                        operands.push_back(temp);
                    }
                }
            }

            vector<Ciphertext> cts_ops(hamming_weight);
            evaluator->add_many(operands, cts_ops[0]);

            #pragma omp parallel for
            for (int i = 1; i < hamming_weight; i++) {
                Plaintext plain_matrix;
                vector<uint64_t> pod_matrix(batch_encoder->slot_count(), i);
                batch_encoder->encode(pod_matrix, plain_matrix);
                evaluator->sub_plain(cts_ops[0], plain_matrix, cts_ops[i]);
            }

            Ciphertext result_cts;

            if (hamming_weight > 1) {
                uint64_t ceil_log_k = ceil(log2(hamming_weight));

                for (int i = 0; i < ceil_log_k; i++) {
                    #pragma omp parallel for
                    for (int j = 0; j < 1 << (ceil_log_k - 1 - i); j++) {

                        if (j + (1 << (ceil_log_k - 1 - i)) < cts_ops.size()) {
                            evaluator->multiply_inplace(cts_ops[j],
                                cts_ops[j + (1 << (ceil_log_k - 1 - i))]);
                            evaluator->relinearize_inplace(cts_ops[j],
                                *rlk_server);
                        }
                        
                    }
                }

                // find multiplicative inverse of k! mod p
                uint64_t inv = prime_mod_inverse(mod_iter_factorial(hamming_weight, plain_modulus), plain_modulus);

                Plaintext pt;
                vector<uint64_t> pod_matrix(batch_encoder->slot_count(),
                    inv);
                batch_encoder->encode(pod_matrix, pt);
                
                evaluator->multiply_plain(cts_ops[0], pt, result_cts);
            } else {
                result_cts = cts_ops[0];
            }

        ////////////////

        Ciphertext batchedComparison=result_cts;

        vector<Ciphertext> ops(n);
        
        #pragma omp parallel for
        for (int i = 1; i <= n; ++i) {
            evaluator->rotate_rows(batchedComparison, -i,
                *gal_keys_server, ops[i-1]);
        }
        ops.push_back(batchedComparison);

        Ciphertext comparisonResult;
        evaluator->add_many(ops, comparisonResult);
        
        Plaintext pt;
        decryptor->decrypt(comparisonResult, pt);
        vector<uint64_t> res;
        batch_encoder->decode(pt, res);
        this->correct = (res[n] == (encrypted_op <= plain_op));

    }

    void folklore_compare(){

        uint64_t num_cmps_per_row = row_count/num_slots_per_element;
        this->num_cmps = poly_modulus_degree/num_slots_per_element;
        Ciphertext comparisonResult;

 
        vector<uint64_t> attribute_vector;
        for (int i = 0; i < num_cmps_per_row; i++) {
            attribute_vector.push_back(encrypted_op);
        }

        Ciphertext client_input;
        Plaintext plaintext; 
        {
            vector<uint64_t> bit_reps(batch_encoder->slot_count(), plain_modulus - 1);
            for(size_t i = 0; i < attribute_vector.size(); i++) {
                int idx = i * num_slots_per_element + n;
                for(int j = 0; j < n; ++j) {
                    bit_reps[idx + j] = plain_modulus - ((attribute_vector[i] >> j) & 1);
                    bit_reps[idx + j + row_count] = bit_reps[idx + j];
                }
            }   
            batch_encoder->encode(bit_reps, plaintext);
        }
        encryptor->encrypt(plaintext, client_input);

        vector<uint64_t> server_values;
        for (int i = 0; i < num_cmps_per_row; i++) {
            server_values.push_back(plain_op);
        }

        Plaintext a; 
        {
            vector<uint64_t> bit_reps(batch_encoder->slot_count(), 1);
            for(size_t i = 0; i < server_values.size(); i++) {
                int idx = i * num_slots_per_element + n;
                for(int j = 0; j < n; ++j) {
                    bit_reps[idx + j] = ((server_values[i] + 1) >> j) & 1;
                    bit_reps[idx + j + row_count] = bit_reps[idx + j];
                }
            }   
            batch_encoder->encode(bit_reps, a);
        }

        // batch_encoder->encode(vector<uint64_t>(poly_modulus_degree, 1), a);

        const auto &neg_b = client_input;
        Plaintext one("1");
        Ciphertext eq, gt;
        evaluator->add_plain(neg_b, a, eq);
        evaluator->square_inplace(eq);
        evaluator->relinearize_inplace(eq, *rlk_server);
        evaluator->negate_inplace(eq);
        evaluator->add_plain_inplace(eq, one);
        evaluator->add_plain(neg_b, one, gt);
        evaluator->multiply_plain_inplace(gt, a);
        evaluator->relinearize_inplace(gt, *rlk_server);

        vector<Ciphertext> eqShifts;
        eqShifts.push_back(gt);
        for(int i = 1; i < n; ++i) {
            eqShifts.emplace_back();
            evaluator->rotate_rows(eq, i, *gal_keys_server, eqShifts[i]);
        }

        evaluator->multiply_many(eqShifts, *rlk_server, comparisonResult);
        eqShifts[0] = comparisonResult;
        for(int i = 1; i < n; ++i) {
            evaluator->rotate_rows(comparisonResult, i, *gal_keys_server, eqShifts[i]);
        }

        evaluator->add_many(eqShifts, comparisonResult);
        evaluator->mod_switch_to_next_inplace(comparisonResult);

        Plaintext pt;
        decryptor->decrypt(comparisonResult, pt);
        vector<uint64_t> res;
        batch_encoder->decode(pt, res);
        this->correct = (res[n] == (encrypted_op <= plain_op));
    }

    void xcmp(Ciphertext enc_a, uint64_t b, Ciphertext& dest) {
        // SEAL does not support multiplying plain zero polynomial, hence the case split
        if(b < poly_modulus_degree - 1) {
            Plaintext T(poly_modulus_degree - b);
            for(size_t deg = 1; deg <= poly_modulus_degree - b - 1; deg++) {
                T.data()[deg] = plain_modulus - 1;
            }

            evaluator->multiply_plain(enc_a, T, dest);
            evaluator->relinearize_inplace(dest, *rlk_server);
        } else {
            encryptor->encrypt(string("0"), dest);
        }
    }

    void xcmp_eq(const Ciphertext &enc_a, uint64_t b, Ciphertext &dest) {
        Plaintext p(uint64_to_hex_string(inv));
        if(b > 0) {
            p = Plaintext(poly_modulus_degree - b + 1);
            p.data()[poly_modulus_degree - b] = plain_modulus - inv;
        }
        evaluator->multiply_plain(enc_a, p, dest);

        Ciphertext tmp;
        for(size_t i = 1; i < poly_modulus_degree; i *= 2) {
            evaluator->apply_galois(dest, poly_modulus_degree / i + 1, *gal_keys_server, tmp);
            evaluator->add_inplace(dest, tmp);
        }
    }

    void two_ctxt_xcmp(const Ciphertext &enc_a1, const Ciphertext &enc_a0, uint64_t b1, uint64_t b0, Ciphertext &dest) {
        xcmp(enc_a0, b0, dest);

        Ciphertext tmp;
        xcmp_eq(enc_a1, b1, tmp);
        evaluator->multiply_inplace(dest, tmp);
        evaluator->relinearize_inplace(dest, *rlk_server);
        xcmp(enc_a1, b1, tmp);
        evaluator->add_inplace(dest, tmp);
    }

    void three_ctxt_xcmp(const Ciphertext &enc_a2, const Ciphertext &enc_a1, Ciphertext &enc_a0, uint64_t b2, uint64_t b1, uint64_t b0, Ciphertext &dest) {
        two_ctxt_xcmp(enc_a1, enc_a0, b1, b0, dest);

        Ciphertext tmp;
        xcmp_eq(enc_a2, b2, tmp);
        evaluator->multiply_inplace(dest, tmp);
        evaluator->relinearize_inplace(dest, *rlk_server);
        xcmp(enc_a2, b2, tmp);
        evaluator->add_inplace(dest, tmp);
    }

    void xcmp_compare(){
        this->num_cmps = 1;
        
        uint64_t b0 = plain_op % poly_modulus_degree;
        uint64_t b1 = (plain_op / poly_modulus_degree) % poly_modulus_degree;
        uint64_t b2 = (plain_op / poly_modulus_degree / poly_modulus_degree) % poly_modulus_degree;

        auto encode = [](uint64_t n) {
            return (n == 0 ? string("1") : "1x^" + to_string(n));
        };

        Plaintext pt0(encode(encrypted_op % poly_modulus_degree));
        Plaintext pt1(encode((encrypted_op / poly_modulus_degree) % poly_modulus_degree));
        Plaintext pt2(encode((encrypted_op / poly_modulus_degree / poly_modulus_degree) % poly_modulus_degree));

        Ciphertext enc_a2, enc_a1, enc_a0;
        encryptor->encrypt(pt2, enc_a2);
        encryptor->encrypt(pt1, enc_a1);
        encryptor->encrypt(pt0, enc_a0);
        Ciphertext ans;
        if(n <= 12) {
            xcmp(enc_a0, b0, ans);
        } else if(n <= 24) {
            two_ctxt_xcmp(enc_a1, enc_a0, b1, b0, ans);
        } else if(n <= 36) {
            three_ctxt_xcmp(enc_a2, enc_a1, enc_a0, b2, b1, b0, ans);
        } else {
            throw "bitlength too large";
        }

        Plaintext pt;
        decryptor->decrypt(ans, pt);
        correct = (pt.data()[0] == (encrypted_op > plain_op));
        
    }

    void batched_comparison(){

        Ciphertext comparisonResult;

        Timer timer;
        timer.start();

        switch (comparison){
            case RANGE_COVER:
                range_cover_compare();
                break;
            case FOLKLORE:
                folklore_compare();
                break;
            case XCMP:
                xcmp_compare();
                break;
        }

        timer.end();

        // correct,bitlength,comparison,log_deg,hamming_weight,code_length,num_cmps,total_time,amortized_time
        cout 
             << correct << ","
             << n << ","
             << comparison << ","
             << poly_modulus_degree << ","
             << hamming_weight << ","
             << code_length << ","
             << num_cmps << ","
             << timer.get_time_in_microseconds() << ","
             << (float)timer.get_time_in_microseconds()/num_cmps << endl;
        // cout << endl;
        // cout << "\t    Total Time: " << timer.get_time_in_microseconds() << endl;
        // cout << "\tAmortized Time: " << (float)timer.get_time_in_microseconds()/num_cmps << endl;
        // return comparisonResult;

    }

};



int main(int argc, char* argv[]) {

    int runs = 10;
    if (argc>1) runs = atoi(argv[1]);

    omp_set_num_threads(32);
    
    ComparisonBenchmark* bench;
    for (int run=0;run<runs;run++){
        for (int n=4;n<=60;n+=2){

            if (n<=36){
                bench = new ComparisonBenchmark();
                bench->init(XCMP, n, 0);
                bench->batched_comparison();
                delete bench;
            }

            for (uint64_t hw: {2,4,8,16,32}){
                if ((float)n/hw>=9) continue;
                if (2*hw>=n+1) continue;
                bench = new ComparisonBenchmark();
                bench->init(RANGE_COVER, n, hw);
                bench->batched_comparison();
                delete bench;
            }

            bench = new ComparisonBenchmark();
            bench->init(FOLKLORE, n, 0);
            bench->batched_comparison();
            delete bench;

        }
    }


}
