#ifndef __CLIENT__H
#define __CLIENT__H

#include "seal/seal.h"
#include "server.h"
#include "tree_utils.h"
#include "utils.h"
#include <cassert>
#include <set>
#include <stdlib.h>
#include <thread>
#include <unistd.h>
using namespace std;
using namespace seal;

class Client {
public:
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    EncryptionParameters* parms;
    uint64_t log_poly_mod_degree;
    uint64_t coefficient_hex_length;
    SEALContext* context;
    KeyGenerator* keygen;
    SecretKey sk;
    Evaluator* evaluator;
    Encryptor* encryptor;
    BatchEncoder* batch_encoder;
    stringstream parms_stream;
    stringstream data_stream;
    stringstream sk_stream;

    TreeEvaluationServer* server;

    long comm_relin = 0, comm_gals = 0, comm_query = 0, comm_response = 0;

    // Preparation steps

    Client();

    void set_server(TreeEvaluationServer* server);
    void setup_crypto(size_t log_poly_modulus_degree = 13, uint64_t prime_bitlength = 20, bool _verbose = true);
    void send_evaluation_keys(QueryParameters* query_parameters);

    void compRangeSmallest(uint64_t a, uint64_t largest, uint64_t* range, int n, bool return_early);
    void compRangeLargest(uint64_t a, uint64_t largest, uint64_t* range, int n);

    vector<Ciphertext> attribute_vector_to_ciphertext(vector<uint64_t> attribute_vector, QueryParameters* query_parameters);
    void send(vector<Ciphertext>& range_bits_rangecode);
    std::vector<Plaintext> load_and_decrypt(QueryParameters* query_parameters);
    uint64_t extract_response_to_classification(std::vector<Plaintext>& response, QueryParameters* query_parameters);

    void submit_classification_with_params(vector<uint64_t> attribute_vector, QueryParameters* query_parameters, bool _verbose);

    bool end_to_end_evaluation(QueryParameters* query_parameters, vector<uint64_t> input, bool _verbose = true);
};

#endif
