#ifndef __UTILS__H
#define __UTILS__H

#include "seal/seal.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <time.h>
#include <vector>
#include <chrono>

using namespace std;

class Timer
{
    public:
     std::chrono::steady_clock::time_point start_;
     std::chrono::steady_clock::time_point end_;
      Timer();
      void start();
      void end();
      long double end_and_get();
      void reset();
      long double get_time_in_milliseconds();
      long double get_time_in_microseconds();
};


template <typename T>
std::string int_to_hex(T i);

std::string vector_to_poly_string(vector<uint64_t> coeffs);

uint64_t hex_string_to_uint(string str);

pair<uint64_t, uint64_t> monomial_to_integer(string monomial);

vector<uint64_t> poly_string_to_vector(std::string poly_string, uint64_t max_size);

/*
    return a ^ e mod n
*/
uint64_t mod_exp(uint64_t a, uint64_t e, uint64_t n);

uint64_t prime_mod_inverse(uint64_t a, uint64_t n);

/*
    return n!
*/
uint64_t iter_factorial(uint64_t n);
uint64_t mod_iter_factorial(uint64_t n, uint64_t mod);

/*
    return n choose k
*/
uint64_t choose(uint64_t n, uint64_t k);

uint64_t find_choose(uint64_t k, uint64_t n);

float log2_choose(uint64_t n, uint64_t k);

uint64_t find_log2_choose(uint64_t k, uint64_t log_n);

/*
    return smallest power of 2 bigger than x
*/
uint64_t bigger_power2(uint64_t x);

std::vector<uint64_t> get_OP_CHaW_encoding(uint64_t __number, uint64_t encoding_size, uint64_t hamming_weight, bool __verbose = true);

vector<pair<int, int>> read_balanced_tree_from_file(string model_filename);
vector<uint64_t> read_input_from_file(string attribute_vec_filename);
vector<seal::Modulus> get_coeff_modulus(uint64_t poly_modulus_degree);

uint64_t find_code_length(uint64_t n, uint64_t k);

#endif
