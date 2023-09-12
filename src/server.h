#ifndef __SERVER__H
#define __SERVER__H

#include "seal/seal.h"
#include "tree_utils.h"
#include "utils.h"
#include <queue>

using namespace std;
using namespace seal;

class TreeEvaluationServer {
public:
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    seal::EncryptionParameters* parms;
    seal::SEALContext* context;

    // Tree Root (for unbalanced case)
    DecisionTreeNode* tree_root;

    vector<DecisionTreeNode> decision_tree;
    vector<DecisionTreeNode> transformed_decision_tree;

    Evaluator* evaluator;
    RelinKeys* rlk_server;
    GaloisKeys* gal_keys_server;
    Decryptor* noise_calculator;
    Encryptor* enc;

    BatchEncoder* batch_encoder;

    TreeEvaluationServer();

    void initialize(uint64_t poly_modulus_degree = 8192);
    void initialize_model(std::string, QueryParameters*);
    // void initialize_model_preorder(std::vector<std::pair<int, int>> decision_tree, QueryParameters* q_params);
    void initialize_params_with_input(stringstream& parms_stream);
    void load_keys(stringstream& data_stream);
    void set_parameters_used_for_debug(stringstream& sk_stream);
    void set_params(stringstream& parms_stream, stringstream& sk_stream, bool _verbose = true);
    vector<Ciphertext> constant_weight_arith(std::vector<std::vector<Plaintext>>& cts_x_bits_,
        std::vector<std::vector<Ciphertext>>& cts_y_bits_, QueryParameters* query_parameters);
    Ciphertext folklore_comparison(
        std::vector<Plaintext>& cts_x_bits_,
        std::vector<Ciphertext>& cts_y_bits_,
        int p, bool parallel);
    vector<Plaintext> prepare_server_comparison_encoding(vector<uint64_t>& server_values, QueryParameters* q_params);
    vector<Ciphertext> batched_comparison(vector<Ciphertext>& client_input, vector<vector<uint64_t>>& server_comp_values, QueryParameters* q_params);
    void respond_to_classification(stringstream& data_stream, QueryParameters* q_params, bool _verbose = true);
    vector<Ciphertext> bfs_path_summing(vector<Ciphertext>& reference_ciphertext, QueryParameters* query_parameters);
    Ciphertext prepare_client_result(vector<Ciphertext>& batched_result, QueryParameters* q);
    // Ciphertext rotate_sum(Ciphertext& sum_paths, QueryParameters* q_params);
    // std::vector<std::pair<Ciphertext, Ciphertext>> compress_return_vec(std::vector<std::pair<Ciphertext, Ciphertext>>& uncompressed_ret, QueryParameters* q_params);
    Plaintext create_random_mask(vector<int>& indices_to_keep, vector<int>& values_for_kept_indices, QueryParameters* q_params);
    Plaintext create_near_zero_mask(vector<int>& indices_to_keep, vector<int>& values_for_kept_indices, QueryParameters* q_params);

    vector<vector<uint64_t>> generate_server_comp_values(QueryParameters* query_parameters);
};

#endif
