#ifndef __TREEUTILS_H
#define __TREEUTILS_H

#include "seal/seal.h"
#include "utils.h"
#include <cassert>
#include <algorithm>
#include <map>
#include <queue>
#include <random>

using namespace seal;

enum PathEvaluation {
    SUM=0,
    PROD=2
};

enum ComparisonType {
    RANGE_COVER=0,
    FOLKLORE=1,
    XCMP=2
};


class DecisionTreeNode {
public:
    uint64_t attribute_used;
    uint64_t threshold_value;

    uint64_t position_in_comparison_vector;
    uint64_t ct_ind;
    uint64_t rotation_amount;

    uint64_t first_index, second_index;

    DecisionTreeNode* left_child;
    DecisionTreeNode* right_child;

    DecisionTreeNode(uint64_t tv, uint64_t au)
    {
        attribute_used = au;
        threshold_value = tv;
        left_child = NULL;
        right_child = NULL;
    }

    void set_position(uint64_t pos)
    {
        this->position_in_comparison_vector = pos;
    }

    void set_position(uint64_t first_index, uint64_t second_index)
    {
        this->first_index=first_index;
        this->second_index=second_index;
    }


    // void set_position(uint64_t pos, uint64_t max_cmp_per_ct){
    //     this->position_in_comparison_vector = pos;
    //     this->ct_ind = pos / max_cmp_per_ct;
    //     this->rotation_amount = (pos % max_cmp_per_ct);
    // }

    bool is_leaf(){
        return (this->left_child==NULL) && (this->right_child==NULL);
    }

    uint64_t classification_value(){
        assert(this->is_leaf() && "Non-leaf node does not have classification value");
        return this->threshold_value;
    }
};

class QueryParameters {
public:
    // Encryption Parameters
    uint64_t log_poly_mod_degree; // Set by Hamming weight
    uint64_t poly_mod_degree;
    int row_count;
    uint64_t prime_bitlength = 20;
    uint64_t num_output_ciphers = 1;
    uint64_t num_input_ciphers;

    size_t slot_count = 0;

    // Secure mode only
    int max_repetitions = 0;
    int num_interal_nodes = 0;

    int max_comps_per_ct = 0;
    int reps_in_ct_per_attr = 0;

    // Common Params
    int k; // Hamming weight
    int n; // Bitlength of thresholds
    int num_slots_per_element;
    int code_length;
    PathEvaluation path_eval; // Balanced or non-balanced version (TODO: Change this variable name!)
    ComparisonType comparison;
    // bool compress_return; // in the non-balanced version, should we compress the return (slower but better space efficiency)
    string write_to_file; // write results to a file for further analysis

    int num_attr;
    int num_leaves; // when compressing, need this to check the correct number of slots

    // Performance params
    map<string, uint64_t> metrics_;
    uint64_t num_threads = 0;

    QueryParameters(int n, int depth, int k, uint64_t num_threads = 0, PathEvaluation path_eval = SUM, ComparisonType comparison=RANGE_COVER, string write_to_file = "results/");
    void set_tree_params(DecisionTreeNode* root);
    void add_parameters_to_metrics();
    
};

class BFSQueueElement {
public:
    DecisionTreeNode* tree_node;
    Ciphertext current_node_ciphertext;
    int path_length;
    BFSQueueElement(DecisionTreeNode* tn, Ciphertext cnc, int pl)
    {
        tree_node = tn;
        current_node_ciphertext = cnc;
        path_length = pl;
    }
};

void permute_ret_vec(std::vector<std::pair<Ciphertext, Ciphertext>>& ret_vec);
std::vector<uint64_t> get_leaf_classification_values(DecisionTreeNode* root);
#endif
