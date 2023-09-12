#include "client.h"

string exec(const char* cmd)
{
    array<char, 128> buffer;
    string result;
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

Client::Client() { this->server = new TreeEvaluationServer; }

void Client::set_server(TreeEvaluationServer* _server)
{
    this->server = _server;
}

void Client::setup_crypto(size_t log_poly_modulus_degree,
    uint64_t prime_bitlength, bool _verbose)
{
    /* Parameter Selection */
    this->parms = new EncryptionParameters(scheme_type::bfv);
    this->parms->set_poly_modulus_degree(1 << log_poly_modulus_degree);
    this->log_poly_mod_degree = log_poly_modulus_degree;
    //this->parms->set_coeff_modulus(
    //CoeffModulus::BFVDefault(1 << log_poly_modulus_degree));
    uint64_t poly_mod_degree = 1 << log_poly_modulus_degree;
    this->parms->set_coeff_modulus(get_coeff_modulus(poly_mod_degree));

    uint64_t coeff_bitcount = 0;
    for (auto& prime : this->parms->coeff_modulus()) {
        coeff_bitcount += prime.bit_count();
    }

    this->parms->set_plain_modulus(
        PlainModulus::Batching(1 << log_poly_mod_degree, prime_bitlength));
    this->parms->save(parms_stream);

    if (_verbose) {
        cout << "-------------------------- Crypto Parameters "
                "--------------------------"
             << endl
             << "\tPoly Mod Degree: " << this->parms->poly_modulus_degree() << endl
             << "\tPlain Modulus: " << this->parms->plain_modulus().value() << " ("
             << this->parms->plain_modulus().bit_count() << " bits)" << endl
             << "\tCiphertext Modulus Bitcount: " << coeff_bitcount << endl;
    }

    this->context = new SEALContext(*(this->parms));
    this->keygen = new KeyGenerator(*(this->context));
    PublicKey pk;
    keygen->create_public_key(pk);
    sk = keygen->secret_key();
    this->evaluator = new Evaluator(*(this->context));
    this->encryptor = new Encryptor(*(this->context), (pk));
    this->batch_encoder = new BatchEncoder(*(this->context));
    this->sk = sk;
}

void Client::send_evaluation_keys(QueryParameters *query_parameters)
{
    Serializable<RelinKeys> rlk_client = keygen->create_relin_keys();
    this->comm_relin += rlk_client.save(data_stream);
    auto gal_keys_client = keygen->create_galois_keys();
    comm_gals += gal_keys_client.save(data_stream);
    this->sk.save(sk_stream);
}
// computes the range cover from [0, a]
void Client::compRangeSmallest(uint64_t a, uint64_t largest, uint64_t* range, int n, bool return_early)
{

    for (int i = 0; i < n + 1; i++) {
        range[i] = largest + 1;
    }
    if (return_early) {
        return;
    }
    // compute height of last common node in the path

    vector<int> exponents;
    int sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    int powers_sum = 0;
    for (auto& q : exponents)
        powers_sum += pow(2, q);

    while (a + 1 - powers_sum != 0) {
        int i = a - powers_sum;
        int j = 0;
        while (pow(2, j) - 1 <= i) {
            j++;
        }
        --j;
        int value = 0;

        for (auto& k : exponents) {
            value += pow(2, k - j);
        }

        range[j] = value;

        exponents.push_back(j);

        powers_sum = 0;
        for (auto& q : exponents)
            powers_sum += pow(2, q);
        int sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    }
}
// computes range cover from [a, 2^n - 1]
void Client::compRangeLargest(uint64_t a, uint64_t largest, uint64_t* range, int n)
{
    for (int i = 0; i < n + 1; i++) {
        range[i] = largest + 1;
    }
    // compute height of last common node in the path

    vector<int> exponents;
    uint64_t powers_sum = 0;
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
        powers_sum = 0;
        for (auto& q : exponents)
            powers_sum += pow(2, q);
        curr_largest = largest - powers_sum;
    }
}



/*
    Input: a vector of elements that need to be encoded into bits
    Output: The encode input which has also been encrypted. The underlying plaintext should have the following format:
            
            input:   [           A[0]                       A[1]                          A[2]                 ...  ]
            output:  [<encoding of A[0]> <BUFFER> <encoding of A[1]> <BUFFER> <encoding of A[2]> <BUFFER>      ...  ]

    The output can be multiple ciphertexts, if need be.
*/
vector<Ciphertext>
Client::attribute_vector_to_ciphertext(
    vector<uint64_t> attribute_vector,
    QueryParameters* query_parameters
){

    int n = query_parameters->n;
    int num_slots_per_element = query_parameters->num_slots_per_element;
    // int max_comps_per_ct = query_parameters->max_comps_per_ct;
    int row_count = query_parameters->row_count;

    // assert(row_count >= attribute_vector.size() * num_slots_per_element); // So that the slots don't go into the next row
    std::vector<Ciphertext> input_ciphertexts;

    if (query_parameters->comparison==RANGE_COVER){
        input_ciphertexts.resize(query_parameters->code_length);

        // middle vector in place of the array from old code
        vector<vector<vector<uint64_t>>> codes;

        // int k = query_parameters->k;

        // reduces serialization sizes
        this->encryptor->set_secret_key(this->sk);

        // size_t global_max = 0;

        uint64_t largest = pow(2, n) - 1;

        for (int i = 0; i < attribute_vector.size(); i++) {
            uint64_t range[n + 1];
            // if (leqSet.count(i) > 0) {
                compRangeLargest(attribute_vector[i], largest, range, n);
            // } else {
                // compRangeSmallest(attribute_vector[i] - 1, largest, range, n, (attribute_vector[i] == 0));
            // }

            std::vector<std::vector<uint64_t>> rangeCode;

            rangeCode.reserve(n + 1);

            uint64_t m = query_parameters->code_length;

            for (int i = n; i >= 0; --i) {
                if (range[i] == largest + 1) {
                    rangeCode.insert(rangeCode.begin(), std::vector<uint64_t>(m, 0));
                } else {
                    rangeCode.insert(rangeCode.begin(),
                        get_OP_CHaW_encoding(range[i], m, query_parameters->k, true));
                }
            }
            codes.push_back(rangeCode);

        }

        // Constructing plaintexts/ciphertexts from the bits of the range codes
        for (int i = 0; i < query_parameters->code_length; ++i) {
            vector<uint64_t> pod_matrix(query_parameters->slot_count, 0ULL);
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
            input_ciphertexts[i] = ct;
        }

    } else if (query_parameters->comparison==FOLKLORE) {
        Plaintext plaintext; 
        {
            vector<uint64_t> bit_reps(batch_encoder->slot_count(), parms->plain_modulus().value() - 1);
            for(size_t i = 0; i < attribute_vector.size(); i++) {
                int idx = i * query_parameters->num_slots_per_element + n;
                for(int j = 0; j < n; ++j) {
                    bit_reps[idx + j] = parms->plain_modulus().value() - ((attribute_vector[i] >> j) & 1);
                    bit_reps[idx + j + row_count] = bit_reps[idx + j];
                }
            }   
            batch_encoder->encode(bit_reps, plaintext);
        }
        input_ciphertexts.emplace_back();
        encryptor->encrypt(plaintext, input_ciphertexts.back());
    }

    query_parameters->num_input_ciphers = input_ciphertexts.size();
    return input_ciphertexts;
}

void Client::send(vector<Ciphertext>& range_bits_rangecode)
{
    for (int i = 0; i < range_bits_rangecode.size(); i++) {
        this->comm_query += range_bits_rangecode[i].save(data_stream);
    }
}

std::vector< Plaintext>
Client::load_and_decrypt(QueryParameters* query_parameters)
{
    Decryptor decryptor(*(this->context), sk);
    Ciphertext __temp_ct;
    Plaintext temp_1;
    Plaintext temp_2;
    std::vector<Plaintext> ret_vec;
    for (int i = 0; i < query_parameters->num_output_ciphers; i++) {
        this->comm_response += __temp_ct.load(*(this->context), data_stream);
        decryptor.decrypt(__temp_ct, temp_1);
        cout << "Final Noise Level: " << decryptor.invariant_noise_budget(__temp_ct) << endl;

        ret_vec.push_back(temp_1);
    }
    return ret_vec;
}

uint64_t
Client::extract_response_to_classification(std::vector<Plaintext>& response,
    QueryParameters* query_parameters)
{
    int slot_count = query_parameters->slot_count;
    int n = query_parameters->n;
    int limit = slot_count / query_parameters->num_slots_per_element;
    vector<uint64_t> indicator_plain;
    vector<uint64_t> classification_plain;
    for (int k = 0; k < response.size(); k++) {
        batch_encoder->decode(response[k], indicator_plain);
        for (int current_index = 0;current_index < query_parameters->poly_mod_degree; current_index++) {
            if (indicator_plain[current_index] == 0) {
                return indicator_plain[current_index+query_parameters->row_count];
            }
        }
    }
    std::cout << "Didn't find a response" << std::endl;
    return 0;
}

void Client::submit_classification_with_params(
    vector<uint64_t> attribute_vector,
    QueryParameters* query_parameters,
    bool _verbose
){

    Timer time_query;
    time_query.start();

    vector<uint64_t> classification_vec;
    // int max_reps = query_parameters->max_comps_per_ct / query_parameters->num_attr;
    // int max_reps = query_parameters->max_repetitions;
    // cout << max_reps << endl;
    for (int i = 0; i < query_parameters->num_attr; i++) {
        for (int j = 0; j < query_parameters->reps_in_ct_per_attr; j++) {
            classification_vec.push_back(attribute_vector[i]);
        }
    }
    vector<Ciphertext> __query_ciphertexts = attribute_vector_to_ciphertext(classification_vec, query_parameters);
    send(__query_ciphertexts);

    query_parameters->metrics_["time_query"] = time_query.end_and_get();
}

uint64_t find_true_respone(QueryParameters* query_parameters, vector<uint64_t> input, TreeEvaluationServer* server){
    DecisionTreeNode* current = server->tree_root;
    while(!current->is_leaf()){
        if (input[current->attribute_used] <= current->threshold_value){
            current = current->left_child;
        } else {
            current = current->right_child;
        }
    }
    return current->classification_value();
}

// End to End
bool Client::end_to_end_evaluation(QueryParameters* query_parameters,
    vector<uint64_t> input, bool _verbose)
{

    Timer time_setup_crypto, time_server_total, time_server_crypto, time_server_latency, time_extract_response;

    if (_verbose)
        cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                "<<<<<<<<<<<<<<<<<<<<<<"
             << endl;

    time_setup_crypto.start();
        this->setup_crypto(
            query_parameters->log_poly_mod_degree,
            query_parameters->prime_bitlength,
            _verbose
        );
    query_parameters->metrics_["time_setup_crypto"] = time_setup_crypto.end_and_get();

    // Setting up the server
    this->server->set_params(parms_stream, sk_stream, _verbose);
    this->send_evaluation_keys(query_parameters);
    this->server->set_parameters_used_for_debug(sk_stream);

    // Setting up input
    this->submit_classification_with_params(input, query_parameters, true);

    // Invoking server
    time_server_latency.start();
        this->server->respond_to_classification(data_stream, query_parameters, _verbose);
    query_parameters->metrics_["time_server_latency"] = time_server_latency.end_and_get();

    // Extracting Response
    time_extract_response.start();
        std::vector<Plaintext> _response_pts = this->load_and_decrypt(query_parameters);
        uint64_t response = this->extract_response_to_classification(_response_pts, query_parameters);
    query_parameters->metrics_["time_extract_response"] = time_extract_response.end_and_get();
    std::cout << "     Response: " << response << std::endl;

    query_parameters->metrics_["comm_relin"] = this->comm_relin;
    query_parameters->metrics_["comm_gals"] = this->comm_gals;
    query_parameters->metrics_["comm_query"] = this->comm_query;
    query_parameters->metrics_["comm_response"] = this->comm_response;
    query_parameters->add_parameters_to_metrics();

    // Checking Correctness
    uint64_t true_response = find_true_respone(query_parameters, input, server);
    cout << "True Response: " << true_response << endl;
    query_parameters->metrics_["correctness"] = (true_response == response);
    cout << ((true_response == response)? "CORRECT":"WRONG") << endl;

    // Printing timings
    if (_verbose) {
        cout << "------------------------ Timing "
                "----------------------------------------------"
             << endl;
        cout << "\tTotal Server    : " << setw(10)
             << query_parameters->metrics_["time_server_total"] << " ms" << endl;
        cout << "\tTotal Server Crypto    : " << setw(10)
             << query_parameters->metrics_["time_server_crypto"] << " ms" << endl;
        cout << "--------------------- Communication "
                "------------------------------------------"
             << endl
             << "\tData Independant: " << comm_relin / 1000 << " KB (Relin keys) + "
             << comm_gals / 1000 << " KB (Gal Keys)" << endl
             << "\tData Dependant: " << comm_query / 1000 << " KB (Query) + "
             << comm_response / 1000 << " KB (Reponse)" << endl
             << endl
             << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                "<<<<<<<<<<<<<<<<<<<<<<"
             << endl;
    }

    uint64_t filename = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count() % 100000000000;

    string filename_string=query_parameters->write_to_file + to_string(filename) + ".csv";
    cout << "Writing to file " << filename_string << endl;

    // rely on the user to ensure that results directory exists
    ofstream outFile;
    outFile.open(filename_string);
    for (pair<string, uint64_t> metric : query_parameters->metrics_) {
        outFile << metric.first << "," << metric.second << endl;
    }
    return true;
}
