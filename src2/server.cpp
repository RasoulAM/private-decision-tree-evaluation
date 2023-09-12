#include "server.h"
#include "parser.h"

#include <utility>
#include <stack>
#include <fstream>

using namespace std;
using namespace seal;

void Server::fetch_model(const string &fname) {
    Parser p;
    auto root = p.parseTree(fname);

    stack<pair<Node*, int>> unvisited;
    unvisited.emplace(root, 0);
    while(!unvisited.empty()) {
        auto p = unvisited.top();
        auto node = p.first;
        auto depth = p.second;
        unvisited.pop();

        if(!node->is_leaf()) {
            inner_nodes.push_back(node);
            unvisited.emplace(node->right, depth + 1);
            unvisited.emplace(node->left, depth + 1);
        } else {
            node->attr_idx = depth;
            leaves.push_back(node);
        }
    }
}

void Server::load_parameters_and_keys(stringstream& params_stream, bool debug_mode, bool _verbose) {
    enc_params = new EncryptionParameters();
    enc_params->load(params_stream);

    context = new SEALContext(*enc_params);

    relin_keys = new RelinKeys();
    galois_keys = new GaloisKeys();
    public_key = new PublicKey();
    relin_keys->load(*context, params_stream);
    galois_keys->load(*context, params_stream);
    public_key->load(*context, params_stream);

    batch_encoder = new BatchEncoder(*context);
    encryptor = new Encryptor(*context, *public_key);
    evaluator = new Evaluator(*context);

    poly_modulus_degree = enc_params->poly_modulus_degree();
    plain_modulus = enc_params->plain_modulus().value();

    inv = mod_inverse(poly_modulus_degree, plain_modulus);

    if(debug_mode) {
        SecretKey secret_key;
        secret_key.load(*context, params_stream);
        decryptor = new Decryptor(*context, secret_key);
    }

    {
        string hex_poly = uint64_to_hex_string(leaves[0]->threshold + 1);
        for(size_t i = 1; i < leaves.size(); i++) {
            hex_poly = uint64_to_hex_string(leaves[i]->threshold + 1) + "x^" + to_string(i) + " + " + hex_poly;
        }
        encryptor->encrypt(hex_poly, labels);
    }
}

vector<Ciphertext> Server::load_inputs(stringstream &data_stream, size_t input_size) {
    vector<Ciphertext> ret;
    Ciphertext T;
    for(size_t i = 0; i < input_size; i++) {
        T.load(*context, data_stream);
        ret.push_back(T);
    }
    data_stream.str("");
    return ret;
}

void Server::send_results(const vector<Ciphertext>& results, stringstream &data_stream) {
    for(auto &ctxt: results) {
        ctxt.save(data_stream);
    }
}

void Server::xcmp(const Ciphertext &enc_a, uint64_t b, Ciphertext &dest) {
    // SEAL does not support multiplying plain zero polynomial, hence the case split
    if(b < poly_modulus_degree - 1) {
        Plaintext T(poly_modulus_degree - b);
        for(size_t deg = 1; deg <= poly_modulus_degree - b - 1; deg++) {
            T.data()[deg] = plain_modulus - 1;
        }

        evaluator->multiply_plain(enc_a, T, dest);
        evaluator->relinearize_inplace(dest, *relin_keys);
    } else {
        encryptor->encrypt(string("0"), dest);
    }
}

void Server::xcmp_eq(const Ciphertext &enc_a, uint64_t b, Ciphertext &dest) {
    Plaintext p(uint64_to_hex_string(inv));
    if(b > 0) {
        p = Plaintext(poly_modulus_degree - b + 1);
        p.data()[poly_modulus_degree - b] = plain_modulus - inv;
    }
    evaluator->multiply_plain(enc_a, p, dest);

    Ciphertext tmp;
    for(size_t i = 1; i < poly_modulus_degree; i *= 2) {
        evaluator->apply_galois(dest, poly_modulus_degree / i + 1, *galois_keys, tmp);
        evaluator->add_inplace(dest, tmp);
    }
}

void Server::two_ctxt_xcmp(const Ciphertext &enc_a1, const Ciphertext &enc_a2, uint64_t b1, uint64_t b2, Ciphertext &dest) {
    xcmp(enc_a2, b2, dest);

    Ciphertext tmp;
    xcmp_eq(enc_a1, b1, tmp);
    evaluator->multiply_inplace(dest, tmp);
    evaluator->relinearize_inplace(dest, *relin_keys);
    xcmp(enc_a1, b1, tmp);
    evaluator->add_inplace(dest, tmp);
}

void Server::make_room_for_packing(Ciphertext &ctxt) {
    Ciphertext tmp;
    for(size_t i = 1; i < leaves.size(); i *= 2) {
        evaluator->apply_galois(ctxt, poly_modulus_degree / i + 1, *galois_keys, tmp);
        evaluator->add_inplace(ctxt, tmp);
    }
}

uint64_t Server::evaluate_plain(const vector<uint64_t> &attr_vec) {
    auto cur = inner_nodes[0];
    while(!cur->is_leaf()) {
        if(attr_vec[cur->attr_idx] <= cur->threshold) {
            cur = cur->left;
        } else {
            cur = cur->right;
        }
    }
    return cur->threshold;
}

void Server::do_server_computation(stringstream& data_stream, size_t input_size, Metrics* metrics, size_t num_ctxt_per_value, bool mult_path, bool _verbose) {
    const vector<Ciphertext> attr_vec = load_inputs(data_stream, input_size);
    metrics->metrics_["inner_nodes"]=inner_nodes.size();
    
    Timer server_time;
    server_time.start();

    Timer gen_comp_bits_time;
    gen_comp_bits_time.start();

    #pragma omp parallel for
    for(auto &node: inner_nodes) {
        if(num_ctxt_per_value == 1) {
            xcmp(attr_vec[node->attr_idx], node->threshold, node->right->ctxt);
        } else if(num_ctxt_per_value == 2) {
            const auto &enc_a1 = attr_vec[2 * node->attr_idx];
            const auto &enc_a2 = attr_vec[2 * node->attr_idx + 1];
            auto b1 = node->threshold / DomainExtensionModulo, b2 = node->threshold % DomainExtensionModulo;
            
            two_ctxt_xcmp(enc_a1, enc_a2, b1, b2, node->right->ctxt);

        } else if(num_ctxt_per_value == 3) {
            const auto &enc_a1 = attr_vec[3 * node->attr_idx];
            const auto &enc_a2 = attr_vec[3 * node->attr_idx + 1];
            const auto &enc_a3 = attr_vec[3 * node->attr_idx + 2];

            auto b1 = (node->threshold / DomainExtensionModulo) / DomainExtensionModulo;
            auto b2 = (node->threshold / DomainExtensionModulo) % DomainExtensionModulo;
            auto b3 = node->threshold % DomainExtensionModulo;

            auto &dest = node->right->ctxt;
            two_ctxt_xcmp(enc_a2, enc_a3, b2, b3, dest);

            Ciphertext tmp;
            xcmp_eq(enc_a1, b1, tmp);
            evaluator->multiply_inplace(dest, tmp);
            evaluator->relinearize_inplace(dest, *relin_keys);
            xcmp(enc_a1, b1, tmp);
            evaluator->add_inplace(dest, tmp);
        }
    }

    for(auto &node: inner_nodes) {
        evaluator->negate(node->right->ctxt, node->left->ctxt);
        evaluator->add_plain_inplace(node->left->ctxt, string("1"));
    }

    metrics->metrics_["time_generate_comparison_bits"] = gen_comp_bits_time.end_and_get();

    // aggregate comparison bits
    for(size_t i = 1; i < inner_nodes.size(); i++) {
        auto &node = inner_nodes[i];
        evaluator->add_inplace(node->left->ctxt, node->ctxt);
        evaluator->add_inplace(node->right->ctxt, node->ctxt);
    }

    // extract results at leaves
    for(auto &node: leaves) {
        evaluator->sub_plain_inplace(node->ctxt, uint64_to_hex_string(node->attr_idx));
        make_room_for_packing(node->ctxt);
    }

    // pack indicators
    Ciphertext indicators;
    indicators = leaves[0]->ctxt;
    for(size_t i = 1; i < leaves.size(); i++) {
        auto &node = leaves[i];
        uint64_t mask = rand_uint64() % (poly_modulus_degree - 1) + 1;
        evaluator->multiply_plain_inplace(node->ctxt, uint64_to_hex_string(mask) + "x^" + to_string(i));
        evaluator->relinearize_inplace(node->ctxt, *relin_keys);
        evaluator->add_inplace(indicators, node->ctxt);
    }

    metrics->metrics_["time_server"] = server_time.end_and_get();

    send_results({labels, indicators}, data_stream);
}
