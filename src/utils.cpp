#include "utils.h"
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
#include <cassert>

using namespace std;

Timer::Timer(){
    start_ = chrono::steady_clock::now();
    end_ = chrono::steady_clock::now();
}
void Timer::start(){
    start_ = chrono::steady_clock::now();
}

void Timer::end(){
    end_ = chrono::steady_clock::now();
}

long double Timer::end_and_get(){
    end_ = chrono::steady_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end_ - start_);
    return elapsed.count();
}

long double Timer::get_time_in_milliseconds(){
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end_ - start_);
    return elapsed.count();
}

long double Timer::get_time_in_microseconds(){
    auto elapsed = chrono::duration_cast<chrono::microseconds>(end_ - start_);
    return elapsed.count();
}

template <typename T>
string int_to_hex(T i)
{
    stringstream stream;
    stream << hex << i;
    return stream.str();
}

// string xx[] = {"0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","10","11","12","13","14","15","16","17","18","19","1a","1b","1c","1d","1e","1f","20","21","22","23","24","25","26","27","28","29","2a","2b","2c","2d","2e","2f","30","31","32","33","34","35","36","37","38","39","3a","3b","3c","3d","3e","3f","40","41","42","43","44","45","46","47","48","49","4a","4b","4c","4d","4e","4f","50","51","52","53","54","55","56","57","58","59","5a","5b","5c","5d","5e","5f","60","61","62","63","64","65","66","67","68","69","6a","6b","6c","6d","6e","6f","70","71","72","73","74","75","76","77","78","79","7a","7b","7c","7d","7e","7f","80","81","82","83","84","85","86","87","88","89","8a","8b","8c","8d","8e","8f","90","91","92","93","94","95","96","97","98","99","9a","9b","9c","9d","9e","9f","a0","a1","a2","a3","a4","a5","a6","a7","a8","a9","aa","ab","ac","ad","ae","af","b0","b1","b2","b3","b4","b5","b6","b7","b8","b9","ba","bb","bc","bd","be","bf","c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","ca","cb","cc","cd","ce","cf","d0","d1","d2","d3","d4","d5","d6","d7","d8","d9","da","db","dc","dd","de","df","e0","e1","e2","e3","e4","e5","e6","e7","e8","e9","ea","eb","ec","ed","ee","ef","f0","f1","f2","f3","f4","f5","f6","f7","f8","f9","fa","fb","fc","fd","fe","ff"};

// string int_to_hex_2( int num ){
//     string _ans="";
//     while (num > 0){
//         _ans = xx[num%256] + _ans;
//         num = num >> 8;
//     }
//     if (_ans.size() == 0)
//         _ans = "0";
//     return _ans;
// }

string vector_to_poly_string(vector<uint64_t> coeffs)
{
    chrono::high_resolution_clock::time_point time1 = chrono::high_resolution_clock::now();
    string ans = int_to_hex((int)(coeffs[0]));
    // cout << int_to_hex((int)(coeffs[0])) << " " << int_to_hex_2((int)(coeffs[0])) << "|" << endl;
    // cout << "---> " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - time1).count() << endl;
    for (int i = 1; i < coeffs.size(); i++) {
        // cout << i << ": " << int_to_hex((int)(coeffs[i])) << " " << ans.length() << endl;
        ans = int_to_hex((int)(coeffs[i])) + "x^" + to_string(i) + " + " + ans;
        // cout << "---> " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - time1).count() << endl;
    }
    // cout << "---> " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - time1).count() << endl;
    return ans;
}

uint64_t hex_string_to_uint(string str)
{
    uint64_t x;
    stringstream ss;
    ss << hex << str;
    ss >> x;
    return x;
}

pair<uint64_t, uint64_t> monomial_to_integer(string monomial)
{
    uint64_t index, coeff;
    if (monomial.find("x^") != string::npos) {
        size_t hit = monomial.find("x^");
        index = stoi(monomial.substr(hit + 2));
        coeff = hex_string_to_uint(monomial.substr(0, hit));
    } else if (monomial.find("x") != string::npos) {
        size_t hit = monomial.find("x");
        index = 1;
        coeff = hex_string_to_uint(monomial.substr(0, hit));
    } else {
        index = 0;
        coeff = hex_string_to_uint(monomial);
    }
    return pair<uint64_t, uint64_t>(index, coeff);
}

vector<uint64_t> poly_string_to_vector(string poly_string, uint64_t max_size)
{
    string __temp = poly_string;
    vector<uint64_t> ans(max_size, 0);
    pair<uint64_t, uint64_t> mono;
    while (__temp.find("+") != string::npos) {
        mono = monomial_to_integer(__temp.substr(0, __temp.find("+")));
        ans[mono.first] = mono.second;
        __temp = __temp.substr(__temp.find("+") + 1);
    }
    mono = monomial_to_integer(__temp);
    ans[mono.first] = mono.second;
    return ans;
}

/*
    return a ^ e mod n
*/
uint64_t mod_exp(uint64_t a, uint64_t e, uint64_t n)
{
    if (e == 0)
        return 1;

    long k = ceil(log2(e));
    uint64_t res;
    res = 1;

    for (long i = k - 1; i >= 0; i--) {
        res = (res * res) % n;
        if ((e / (uint64_t)(pow(2, i))) % 2 == 1)
            res = (res * a) % n;
    }
    return res;
}

uint64_t prime_mod_inverse(uint64_t a, uint64_t n)
{
    return mod_exp(a, n - 2, n);
}

/*
    return n!
*/
uint64_t iter_factorial(uint64_t n)
{
    uint64_t ret = 1;
    for (uint64_t i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}

uint64_t mod_iter_factorial(uint64_t n, uint64_t mod)
{
    uint64_t ret = 1;
    for (uint64_t i = 1; i <= n; ++i){
        ret *= i;
        ret %= mod;
    }
    return ret;
}


/*
    return n choose k
*/
uint64_t choose(uint64_t n, uint64_t k)
{
    if (k > n) {
        return 0;
    }
    uint64_t r = 1;
    for (uint64_t d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}

uint64_t find_choose(uint64_t k, uint64_t n)
{
    uint64_t _ans = k;
    while (choose(_ans, k) < n) {
        _ans++;
    }
    return _ans;
}

float log2_choose(uint64_t n, uint64_t k)
{
    if (k > n) {
        cout << "Can't have k > n" << endl;
        return -1; // Error case
    }
    float log2_r = 0;
    for (uint64_t d = 1; d <= k; ++d) {
        log2_r += log2(n) - log2(d);
        n--;
    }
    return log2_r;
}

uint64_t find_log2_choose(uint64_t k, uint64_t log_n)
{
    uint64_t _ans = k;
    while (log2_choose(_ans, k) < log_n) {
        _ans++;
    }
    return _ans;
}

/*
    return smallest power of 2 bigger than x
*/
uint64_t bigger_power2(uint64_t x)
{
    uint64_t ans = 1;
    while (ans < x)
        ans *= 2;
    return ans;
}

vector<uint64_t> get_OP_CHaW_encoding(uint64_t __number, uint64_t encoding_size, uint64_t hamming_weight, bool __verbose)
{
    vector<uint64_t> ans(encoding_size, 0ULL);
    uint64_t log2_mod_size = log2_choose(encoding_size, hamming_weight);
    assert(log2(__number) < log2_mod_size && "Overflow occurred, everything okay?");
    long remainder = __number, k_prime = hamming_weight;
    for (long pointer = encoding_size - 1; pointer >= 0; pointer--) {
        if (remainder >= choose(pointer, k_prime)) {
            ans[pointer] = 1ULL;
            remainder -= choose(pointer, k_prime);
            k_prime -= 1;
        }
    }
    return ans;
}

/*
    Hash the string __s 
*/
uint64_t hash_(string __s)
{
    hash<string> ptr_hash;
    return ptr_hash(__s);
}

vector<pair<int, int>> read_balanced_tree_from_file(string model_filename)
{

    vector<pair<int, int>> tree;
    ifstream model_infile(model_filename);

    string line;
    while (getline(model_infile, line)) {
        istringstream iss(line);
        int a, b;
        if (!(iss >> a >> b)) {
            break;
        } // error
        tree.push_back(make_pair(a, b));
    }

    return tree;
}

vector<uint64_t> read_input_from_file(string attribute_vec_filename)
{
    vector<uint64_t> input;
    ifstream client_infile(attribute_vec_filename);

    string line;
    while (getline(client_infile, line)) {
        istringstream iss(line);
        uint64_t a;
        if (!(iss >> a)) {
            break;
        } // error

        input.push_back(a);
    }
    return input;
}

vector<seal::Modulus> get_coeff_modulus(uint64_t poly_modulus_degree)
{
    if (poly_modulus_degree == 16384) {
        return seal::CoeffModulus::Create(poly_modulus_degree, { 48, 48, 48, 48, 48, 48, 48 });
    } else {
        return seal::CoeffModulus::BFVDefault(poly_modulus_degree);
    }
}

uint64_t find_code_length(uint64_t n, uint64_t k)
{
    if (k==0) return 0;
    // TODO: this will fail for large values of k cause k! doesn't fit in uint64_t
    uint64_t cur_n = std::pow((uint64_t)2, n);
    uint64_t double_m = std::pow(cur_n, 1 / (double)k);
    // m = std::ceil(std::pow(cur_n * iter_factorial(k), 1 / (double)k) + k);
    for (int i=2;i<=k;i++){
        double_m *= std::pow(i, 1 / (double)k);
    }
    uint64_t m = std::ceil(double_m)+k;
    while(log2_choose(m, k) <= n){
        m++;
    }
    while (log2_choose(m, k) > n) {
        --m;
    }
    ++m;

    return m;
}