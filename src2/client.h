#ifndef __CLIENT__H
#define __CLIENT__H

#include "seal/seal.h"
#include "utils.h"
#include "server.h"

#include <vector>

class Client {
    public:
    Client(Server *_server);
    bool run_protocol(const vector<uint64_t> &attr_vec, const struct Config &config);

    private:
    seal::EncryptionParameters *enc_params = nullptr;
    seal::SEALContext *context = nullptr;
    seal::KeyGenerator *keygen = nullptr;
    seal::SecretKey secret_key;
    seal::Evaluator *evaluator = nullptr;
    seal::Encryptor *encryptor = nullptr;
    seal::BatchEncoder *batch_encoder = nullptr;
    
    // Stringstreams
    std::stringstream params_stream;
    std::stringstream data_stream;

    Server *server = nullptr;
    Metrics *metrics = nullptr;

    // Preparation steps
    void set_server(Server *server);
    void setup_crypto(uint64_t log_poly_modulus_degree, uint64_t prime_bitlength, bool _verbose);
    void send_parameters_and_keys(bool debug_mode);
    vector<seal::Plaintext> encode_attr_vec(const vector<uint64_t> &attr_vec, size_t num_ctxt_per_value);
    void encrypt_and_send(const vector<seal::Plaintext> &plain_input, bool _verbose);

    vector<seal::Plaintext> load_and_decrypt(size_t output_size, bool _verbose);
    uint64_t extract_response(const vector<seal::Plaintext> &results, bool mult_path);
};

#endif
