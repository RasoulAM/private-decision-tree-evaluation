#include "server.h"
#include "parser.h"
#include <cassert>

using namespace std;
using namespace seal;

TreeEvaluationServer::TreeEvaluationServer()
{
}

/*std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}
*/
template <typename T>
std::vector<T>
transform_for_rotation_summing(std::vector<T>& full_tree)
{
    int
        last_layer_size
        = (full_tree.size() + 1) / 2;
    int
        depth
        = log2(full_tree.size() + 1) - 1;

    std::vector<T> rotating_tree;
    rotating_tree.reserve(last_layer_size * depth);

    int current_layer_length = 1;
    int added = 0;
    int total_to_add = last_layer_size * depth;
    int current_layer_start_index = 0;
    while (added != total_to_add) {
        // number of repetitions of each item in a row
        for (int j = 0; j < current_layer_length; ++j) {
            for (int i = 0; i < last_layer_size / current_layer_length; i++) {
                rotating_tree.push_back(
                    full_tree
                        [current_layer_start_index + j]);
                added += 1;
            }
        }
        current_layer_start_index += current_layer_length;
        current_layer_length *= 2;
    }
    return rotating_tree;
}

void TreeEvaluationServer::initialize(uint64_t poly_modulus_degree)
{
    this->parms = new EncryptionParameters(scheme_type::bfv);
    this->parms->set_poly_modulus_degree(poly_modulus_degree);
    this->parms->set_coeff_modulus(get_coeff_modulus(poly_modulus_degree));

    this->parms->set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));

    // this->parms = new EncryptionParameters(scheme_type::bfv);
    // this->parms->load(parms_stream);

    this->context = new SEALContext(*(this->parms));
    KeyGenerator keygen(*context);
    PublicKey pk;
    keygen.create_public_key(pk);

    evaluator = new Evaluator(*context);
    this->rlk_server = new RelinKeys();
    keygen.create_relin_keys(*(this->rlk_server));

    vector<uint32_t> elts;
    for (int i = 1; i < poly_modulus_degree; i *= 2) {
        elts.push_back(2 * i + 1);
    }
    this->gal_keys_server = new GaloisKeys();
    keygen.create_galois_keys(elts, *(this->gal_keys_server));

    this->enc = new Encryptor(*context, pk);

    this->batch_encoder = new BatchEncoder(*context);
}

void TreeEvaluationServer::initialize_model(
    string model_filename,
    QueryParameters* query_parameters
){
    Parser p;
    this->tree_root = p.parseTree(model_filename);
    query_parameters->set_tree_params(this->tree_root);
}

void TreeEvaluationServer::
    initialize_params_with_input(stringstream& parms_stream)
{
    this->parms = new EncryptionParameters();
    this->parms->load(parms_stream);

    parms_stream.seekg(0, parms_stream.beg);
    context = new SEALContext(*(this->parms));

    this->evaluator = new Evaluator(*context);
    this->batch_encoder = new BatchEncoder(*context);

}

void TreeEvaluationServer::load_keys(stringstream& data_stream)
{
    this->rlk_server = new RelinKeys();
    this->gal_keys_server = new GaloisKeys();
    this->rlk_server->load(*context, data_stream);
    this->gal_keys_server->load(*context, data_stream);
}

void TreeEvaluationServer::set_parameters_used_for_debug(stringstream& sk_stream)
{
    SecretKey sk;
    sk.load(*context, sk_stream);
    this->noise_calculator = new Decryptor(*context, sk);
    KeyGenerator keygen(*context, sk);

    PublicKey pk;
    keygen.create_public_key(pk);
    this->enc = new Encryptor(*context, pk);
}

void TreeEvaluationServer::set_params(stringstream& parms_stream,
    stringstream& sk_stream, bool _verbose)
{
    this->initialize_params_with_input(parms_stream);
    bool DEBUG = false;
    if (DEBUG)
        this->set_parameters_used_for_debug(sk_stream);
}

// Constant-weight Equality Operator
vector<Ciphertext> TreeEvaluationServer::constant_weight_arith(
    vector<vector<Plaintext>>& cts_x_bits_,
    vector<vector<Ciphertext>>& cts_y_bits_,
    QueryParameters* query_parameters)
{
    assert(cts_x_bits_.size() == cts_y_bits_.size() && "Size mismatch");
    int num_cts = cts_x_bits_.size();
    int code_length = query_parameters->code_length;
    int hamming_weight = query_parameters->k;

    vector<vector<Ciphertext>> operands(num_cts, vector<Ciphertext>(0));

    #pragma omp parallel for collapse(2)
    for (int ct_ind=0; ct_ind < num_cts; ct_ind++){
        for (int i = 0; i < code_length; i++) {
            if (!cts_x_bits_[ct_ind][i].is_zero()){
                Ciphertext temp;
                this->evaluator->multiply_plain(cts_y_bits_[ct_ind][i], cts_x_bits_[ct_ind][i], temp);
                #pragma omp critical
                {
                    operands[ct_ind].push_back(temp);
                }
            }
        }
    }

    vector<vector<Ciphertext>> cts_ops(num_cts, vector<Ciphertext>(hamming_weight));
    #pragma omp parallel for
    for (int ct_ind=0; ct_ind < num_cts; ct_ind++){
        this->evaluator->add_many(operands[ct_ind], cts_ops[ct_ind][0]);
    }

    #pragma omp parallel for collapse(2)
    for (int ct_ind=0; ct_ind < num_cts; ct_ind++){
        for (int i = 1; i < hamming_weight; i++) {
            Plaintext plain_matrix;
            vector<uint64_t> pod_matrix(batch_encoder->slot_count(), i);
            this->batch_encoder->encode(pod_matrix, plain_matrix);
            this->evaluator->sub_plain(cts_ops[ct_ind][0], plain_matrix, cts_ops[ct_ind][i]);
        }
    }

    vector<Ciphertext> result_cts(num_cts);

    if (hamming_weight > 1) {
        uint64_t ceil_log_k = ceil(log2(hamming_weight));

        for (int i = 0; i < ceil_log_k; i++) {
            #pragma omp parallel for collapse(2)
            for (int ct_ind=0; ct_ind < num_cts; ct_ind++){
                for (int j = 0; j < 1 << (ceil_log_k - 1 - i); j++) {

                    if (j + (1 << (ceil_log_k - 1 - i)) < cts_ops[ct_ind].size()) {
                        this->evaluator->multiply_inplace(cts_ops[ct_ind][j],
                            cts_ops[ct_ind][j + (1 << (ceil_log_k - 1 - i))]);
                        this->evaluator->relinearize_inplace(cts_ops[ct_ind][j],
                            *(rlk_server));
                    }
                    
                }
            }
        }

        // find multiplicative inverse of k! mod p
        uint64_t inv = prime_mod_inverse(iter_factorial(hamming_weight), *(this->parms->plain_modulus().data()));

        Plaintext pt;
        vector<uint64_t> pod_matrix(this->batch_encoder->slot_count(),
            inv);
        this->batch_encoder->encode(pod_matrix, pt);
        
        #pragma omp parallel for
        for (int ct_ind=0; ct_ind < num_cts; ct_ind++){
            this->evaluator->multiply_plain(cts_ops[ct_ind][0], pt, result_cts[ct_ind]);
        }
    } else {
        #pragma omp parallel for
        for (int ct_ind=0; ct_ind < num_cts; ct_ind++){
            result_cts[ct_ind] = cts_ops[ct_ind][0];
        }
    }

    return result_cts;
}

vector<Plaintext> TreeEvaluationServer::prepare_server_comparison_encoding(
    vector<uint64_t>& server_values,
    QueryParameters* query_parameters
){

    int n = query_parameters->n;
    int k = query_parameters->k;

    // inner is the array from previous code
    vector<vector<uint64_t>> paths(server_values.size(), vector<uint64_t>());

    for (int i=0; i < server_values.size(); i++) {
        uint64_t cur = server_values[i];
        for (int j = 0; j < n; j++) {
            paths[i].push_back(cur);
            cur /= 2;
        }
        paths[i].push_back(0);
    }

    uint64_t m = query_parameters->code_length;

    // middle is for array in old code
    std::vector<std::vector<std::vector<uint64_t>>> pathcodes;

    for (int j = 0; j < paths.size(); ++j) {
        vector<vector<uint64_t>> pathcode;
        pathcode.reserve(n + 1);
        for (int i = n; i >= 0; --i) {
            pathcode.insert(pathcode.begin(),
                get_OP_CHaW_encoding(paths[j][i], m, k, true));
        }
        pathcodes.push_back(pathcode);
    }

    std::vector<Plaintext> cts_range_bits_pathcode;

    Plaintext plain_matrix;
    Ciphertext ct;
    for (int i = 0; i < m; ++i) {
        vector<uint64_t> pod_matrix(query_parameters->slot_count, 0ULL);
        for (int k = 0; k < pathcodes.size(); k++) {
            for (int j = 0; j < n + 1; ++j) {
                if (pathcodes[k][j].size() <= i) {
                    continue;
                } else {
                    pod_matrix[k * query_parameters->num_slots_per_element + j] = pathcodes[k][j][i];
                    pod_matrix[k * query_parameters->num_slots_per_element + j + query_parameters->row_count] = pathcodes[k][j][i];
                }
            }
        }
        batch_encoder->encode(pod_matrix, plain_matrix);
        cts_range_bits_pathcode.push_back(plain_matrix);
    }
    return cts_range_bits_pathcode;
}


/*
    INPUT: the client input is mentioned in the client file
           the server input is each element that needs to be compared with the corresponding (encrypted) client element
        
    this function: encode each server value and constuct a plaintext, compare that with the client input.

    Output: If each client input take up (n + BUFFER) bits then the result of the comparison should be in bit (n+1)
*/


vector<Ciphertext> TreeEvaluationServer::batched_comparison(
    vector<Ciphertext>& client_input,
    vector<vector<uint64_t>>& server_comp_values,
    QueryParameters* q_params
){

    int num_cts = ceil((float)q_params->max_repetitions/q_params->reps_in_ct_per_attr);
    vector<Ciphertext> comparisonResult(num_cts);
    int n = q_params->n;

    vector<vector<Ciphertext>> client_input_vec(num_cts, client_input);

    if (q_params->comparison==RANGE_COVER) { // Constant-weight Comparison

        int k = q_params->k;

        vector<vector<Plaintext>> cts_range_bits_pathcode_vec;

        for (int i=0;i<num_cts;i++)
            cts_range_bits_pathcode_vec.push_back(prepare_server_comparison_encoding(server_comp_values[i], q_params));

        vector<Ciphertext> batchedComparison = constant_weight_arith(cts_range_bits_pathcode_vec, client_input_vec,
            q_params);

        vector<vector<Ciphertext>> ops(num_cts, vector<Ciphertext>(n));

        #pragma omp parallel for collapse(2)
        for (int ct_ind=0;ct_ind<num_cts;ct_ind++){
            for (int i = 1; i <= n; ++i) {
                // Ciphertext rotations;
                evaluator->rotate_rows(batchedComparison[ct_ind], -i,
                    *(this->gal_keys_server), ops[ct_ind][i-1]);
            }
        }

        #pragma omp parallel for
        for (int ct_ind=0;ct_ind<num_cts;ct_ind++){
            ops[ct_ind].push_back(batchedComparison[ct_ind]);
            evaluator->add_many(ops[ct_ind], comparisonResult[ct_ind]);
        }

    } else if (q_params->comparison==FOLKLORE) {
        for (int ct_ind=0;ct_ind<num_cts;ct_ind++){
            Plaintext a; 
            {
                vector<uint64_t> bit_reps(batch_encoder->slot_count(), 1);
                for(size_t i = 0; i < server_comp_values[ct_ind].size(); i++) {
                    int idx = i * q_params->num_slots_per_element + n;
                    for(int j = 0; j < n; ++j) {
                        bit_reps[idx + j] = ((server_comp_values[ct_ind][i] + 1) >> j) & 1;
                        bit_reps[idx + j + q_params->row_count] = bit_reps[idx + j];
                    }
                }   
                batch_encoder->encode(bit_reps, a);
            }

            const auto &neg_b = client_input_vec[ct_ind][0];
            Plaintext one("1");
            Ciphertext eq, gt;
            evaluator->add_plain(neg_b, a, eq);
            evaluator->square_inplace(eq);
            evaluator->relinearize_inplace(eq, *rlk_server);
            evaluator->negate_inplace(eq);
            evaluator->add_plain_inplace(eq, one);
            // eq = 1 - (a - b)^2
            evaluator->add_plain(neg_b, one, gt);
            evaluator->multiply_plain_inplace(gt, a);
            // gt = a * (1 - b)
            evaluator->relinearize_inplace(gt, *rlk_server);

            vector<Ciphertext> eqShifts;
            eqShifts.push_back(gt);
            for(int i = 1; i < n; ++i) {
                eqShifts.emplace_back();
                evaluator->rotate_rows(eq, i, *gal_keys_server, eqShifts[i]);
                // eqShifts[i] = the i+1 most significant bits are equal
            }

            evaluator->multiply_many(eqShifts, *rlk_server, comparisonResult[ct_ind]);
            eqShifts[0] = comparisonResult[ct_ind];
            for(int i = 1; i < n; ++i) {
                evaluator->rotate_rows(comparisonResult[ct_ind], i, *gal_keys_server, eqShifts[i]);
            }

            evaluator->add_many(eqShifts, comparisonResult[ct_ind]);
            evaluator->mod_switch_to_next_inplace(comparisonResult[ct_ind]);
        }
    }
    
    return comparisonResult;
}

vector<Ciphertext> TreeEvaluationServer::bfs_path_summing(vector<Ciphertext>& reference_ciphertext, QueryParameters* query_parameters)
{
    int n = query_parameters->n;

    vector<uint64_t> pod_matrix(query_parameters->slot_count, 0);

    Plaintext pt;

    batch_encoder->encode(pod_matrix, pt);

    Ciphertext zero_ciphertext;
    this->enc->encrypt(pt, zero_ciphertext);
    if(query_parameters->comparison == FOLKLORE) {
        evaluator->mod_switch_to_next_inplace(zero_ciphertext);
    }
    vector<Ciphertext> finished_ciphers;

    vector<DecisionTreeNode*> all_queue_elements_rotation_amount;

    std::queue<DecisionTreeNode*> traversal_queue;
    traversal_queue.push(this->tree_root);
    while (!traversal_queue.empty()) {
        DecisionTreeNode* current_element = traversal_queue.front();
        traversal_queue.pop();

        if (!current_element->is_leaf()) {
            all_queue_elements_rotation_amount.push_back(current_element);
            traversal_queue.push(current_element->left_child);
            traversal_queue.push(current_element->right_child);
        }
    }

    vector<Ciphertext> rotated_ciphs(2*all_queue_elements_rotation_amount.size());

    #pragma omp parallel
    {
        Plaintext pt;
        batch_encoder->encode(vector<uint64_t>(query_parameters->poly_mod_degree,1), pt);
        #pragma omp for
        for (int i = 0; i < all_queue_elements_rotation_amount.size(); ++i) {
            
            int ct_ind = all_queue_elements_rotation_amount[i]->second_index / query_parameters->reps_in_ct_per_attr;
            int rotation_amount = (all_queue_elements_rotation_amount[i]->second_index % query_parameters->reps_in_ct_per_attr
                                + all_queue_elements_rotation_amount[i]->first_index * query_parameters->reps_in_ct_per_attr)
                                * query_parameters->num_slots_per_element
                                + n;

            evaluator->rotate_rows(
                reference_ciphertext[ct_ind], rotation_amount,
                *gal_keys_server, rotated_ciphs[2*i]
            );
            evaluator->negate(rotated_ciphs[2*i], rotated_ciphs[2*i+1]);
            evaluator->add_plain_inplace(rotated_ciphs[2*i+1], pt);
        }
    }

    int current_index = 0;
    std::queue<BFSQueueElement> bfs_queue;
    bfs_queue.push(BFSQueueElement(this->tree_root, zero_ciphertext, 0));
    while (!bfs_queue.empty()) {
        BFSQueueElement current_element = bfs_queue.front();
        bfs_queue.pop();
        if (current_element.tree_node->is_leaf()) {
            // both children are NULL

            vector<uint64_t> pod_matrix(query_parameters->slot_count, current_element.path_length);

            Plaintext pt;
            batch_encoder->encode(pod_matrix, pt);
            this->evaluator->sub_plain_inplace(current_element.current_node_ciphertext, pt);
            finished_ciphers.push_back(current_element.current_node_ciphertext);
        } else {

            // accumulate in the n-th slot
            Ciphertext left_child_ciph;
            evaluator->add(rotated_ciphs[current_index], current_element.current_node_ciphertext, left_child_ciph);
            BFSQueueElement left_child_element(current_element.tree_node->left_child, left_child_ciph, current_element.path_length + 1);

            current_index += 1;
            Ciphertext right_child_ciph;
            evaluator->add(rotated_ciphs[current_index], current_element.current_node_ciphertext, right_child_ciph);

            BFSQueueElement right_child_element(current_element.tree_node->right_child, right_child_ciph, current_element.path_length + 1);
            current_index += 1;

            bfs_queue.push(left_child_element);
            bfs_queue.push(right_child_element);
        }
    }
    return finished_ciphers;
}

Ciphertext TreeEvaluationServer::prepare_client_result(
    std::vector<Ciphertext>& batched_result,
    QueryParameters* query_parameters
){

    Ciphertext ret_vec;
    
    int batch_p_num = *(this->parms->plain_modulus().data());

        std::vector<uint64_t> leaf_classification_values = get_leaf_classification_values(this->tree_root);
        
        vector<uint64_t> values_vec(query_parameters->poly_mod_degree,0);

        #pragma omp parallel
        {
            Plaintext pt_mask; 
            #pragma omp for
            for (int i=0;i<batched_result.size();i++){
                vector<int> to_keep({i,i+ query_parameters->row_count});
                vector<int> to_keep_vals({(rand()%(batch_p_num-1))+1,(rand()%(batch_p_num-1))+1});
                pt_mask=create_near_zero_mask(
                    to_keep,
                    to_keep_vals,
                    query_parameters
                );
                evaluator->rotate_rows_inplace(batched_result[i], -i, *gal_keys_server);
                evaluator->multiply_plain_inplace(batched_result[i], pt_mask);

                values_vec[query_parameters->row_count+i] = leaf_classification_values[i];
            }
        }

        evaluator->add_many(batched_result, ret_vec);

        Plaintext values_pt;
        batch_encoder->encode(values_vec, values_pt);

        evaluator->add_plain_inplace(ret_vec, values_pt);

    this->evaluator->mod_switch_to_inplace(ret_vec, context->last_context_data()->parms_id());
    
    return ret_vec;
}

Plaintext TreeEvaluationServer::create_random_mask(vector<int>& indices_to_keep, vector<int>& values_for_kept_indices, QueryParameters* q_params)
{
    Plaintext mask;
    vector<uint64_t> mask_matrix(batch_encoder->slot_count(), 0);

    int batch_p_num = *(this->parms->plain_modulus().data());
    for (int j = 0; j < q_params->slot_count; ++j) {
        mask_matrix[j] = (std::rand() % (batch_p_num - 1)) + 1; // ensures we are multiplying by a non-zero number
    }

    for (int i = 0; i < indices_to_keep.size(); i++) {
        mask_matrix[indices_to_keep[i]] = values_for_kept_indices[i];
    }
    this->batch_encoder->encode(mask_matrix, mask);
    return mask;
}

Plaintext TreeEvaluationServer::create_near_zero_mask(vector<int>& indices_to_keep, vector<int>& values_for_kept_indices, QueryParameters* q_params)
{
    Plaintext mask;
    vector<uint64_t> mask_matrix(batch_encoder->slot_count(), 0);

    for (int i = 0; i < indices_to_keep.size(); i++) {
        mask_matrix[indices_to_keep[i]] = values_for_kept_indices[i];
    }
    this->batch_encoder->encode(mask_matrix, mask);
    return mask;
}

// Annotate the nodes of the tree with their position in the plaintext vector
void traverse_and_fill(DecisionTreeNode* root, vector<vector<uint64_t>>& arrs, map<int, int>& first_unused_index, int max_reps)
{
    if (!root->is_leaf()) {
        int first_index = root->attribute_used;
        int second_index = first_unused_index[root->attribute_used];

        arrs[second_index/max_reps][first_index*max_reps+(second_index%max_reps)] = root->threshold_value;
        first_unused_index[root->attribute_used] += 1;

        // root->set_position(first_index * max_reps + second_index);
        root->set_position(first_index, second_index);

        traverse_and_fill(root->left_child, arrs, first_unused_index, max_reps);
        traverse_and_fill(root->right_child, arrs, first_unused_index, max_reps);
    }
}

vector<vector<uint64_t>> TreeEvaluationServer::generate_server_comp_values(QueryParameters* query_parameters){

    uint64_t number_of_attributes = query_parameters->num_attr;
    int num_cts = ceil((float)query_parameters->max_repetitions/query_parameters->reps_in_ct_per_attr);
    vector<vector<uint64_t>> server_comp_values(num_cts, vector<uint64_t>(number_of_attributes*query_parameters->reps_in_ct_per_attr,0));

    map<int, int> first_unused_indices;
    for (int i = 0; i < number_of_attributes; i++) {
        first_unused_indices[i] = 0;
    }

    traverse_and_fill(this->tree_root, server_comp_values, first_unused_indices, query_parameters->reps_in_ct_per_attr);

    return server_comp_values;
}

void TreeEvaluationServer::respond_to_classification(
    stringstream& data_stream,
    QueryParameters* query_parameters,
    bool _verbose
){

    // Setting threads
    if (query_parameters->num_threads == 0) query_parameters->num_threads = omp_get_max_threads();
    omp_set_num_threads(query_parameters->num_threads);
    if (_verbose)
        cout << "\tNumber of Threads: " << query_parameters->num_threads << endl
             << "--------------------------------------------------------------------"
             << endl;

    Timer time_server_crypto, time_server_batch_compare, time_server_total;
    time_server_total.start();
    
    Timer time_load_data;
    time_load_data.start();
    ////////////////////////////////////////////////////////
        // Loading Keys, then data
        this->load_keys(data_stream);
        vector<Ciphertext> encrypted_query(query_parameters->num_input_ciphers);
        for (int i = 0; i < query_parameters->num_input_ciphers; i++) {
            encrypted_query[i].load(*context, data_stream);
        }
    ////////////////////////////////////////////////////////
    query_parameters->metrics_["time_load_data"] = time_load_data.end_and_get();
    query_parameters->metrics_["noise_budget_after_comparison"] = this->noise_calculator->invariant_noise_budget(encrypted_query[0]);

    time_server_crypto.start();

    // Generate the server values to compare with the clients values
    vector<vector<uint64_t>> server_comp_values = generate_server_comp_values(query_parameters);

    Timer time_batched_comparison;
    time_batched_comparison.start();
    ////////////////////////////////////////////////////////
        vector<Ciphertext> comparisonResult = batched_comparison(encrypted_query, server_comp_values, query_parameters);
    ////////////////////////////////////////////////////////
    query_parameters->metrics_["time_batched_comparison"] = time_batched_comparison.end_and_get();
    query_parameters->metrics_["noise_budget_after_comparison"] = this->noise_calculator->invariant_noise_budget(comparisonResult[0]);
    /*for (int i = 0; i < 5; i++) {
		this->evaluator->mod_switch_to_next_inplace(comparisonResult);
	}*/
    
    // Sum up the path
    Timer time_summing;
    time_summing.start();
    ////////////////////////////////////////////////////////
        vector<Ciphertext> comparisonResultList;
        if (query_parameters->path_eval == SUM) {
            comparisonResultList = bfs_path_summing(comparisonResult, query_parameters);
        } else if (query_parameters->path_eval == PROD) {
            // TODO
        }
    ////////////////////////////////////////////////////////
    query_parameters->metrics_["time_summing"] = time_summing.end_and_get();
    query_parameters->metrics_["noise_budget_after_summing"] = this->noise_calculator->invariant_noise_budget(comparisonResultList[0]);

    // Sum up the path
    Timer time_prepare_result;
    time_prepare_result.start();
    ////////////////////////////////////////////////////////
    // Mask everything except required result
        Ciphertext client_result = prepare_client_result(comparisonResultList, query_parameters);
    ////////////////////////////////////////////////////////
    query_parameters->metrics_["time_prepare_result"] = time_prepare_result.end_and_get();
    query_parameters->metrics_["noise_budget_after_prepare_result"] = this->noise_calculator->invariant_noise_budget(client_result);

    query_parameters->metrics_["time_server_crypto"] = time_server_crypto.end_and_get();

    if (_verbose)
        cout << "--- End of process ---" << endl;

    data_stream.str("");

    auto size_encrypted_answer = 0;
        size_encrypted_answer += client_result.save(data_stream);

    query_parameters->metrics_["time_server_total"] = time_server_total.end_and_get();
}
