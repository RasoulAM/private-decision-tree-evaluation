#include "client.h"

#include <iomanip>
#include <cassert>
#include <fstream>
#include <thread>

using namespace std;
using namespace seal;

Client::Client(Server *_server): server{_server} {}

// Ideal for batch encoding
void Client::setup_crypto(uint64_t log_poly_modulus_degree, uint64_t prime_bitlength, bool _verbose) {
    /* Parameter Selection */
    enc_params = new EncryptionParameters(scheme_type::bfv);
    enc_params->set_poly_modulus_degree(1 << log_poly_modulus_degree);

    // Setting the coeffcient modulus
    enc_params->set_coeff_modulus(CoeffModulus::BFVDefault(1 << log_poly_modulus_degree));
    enc_params->set_plain_modulus(PlainModulus::Batching(1 << log_poly_modulus_degree, prime_bitlength));

    context = new SEALContext(*enc_params);
    keygen = new KeyGenerator(*context);
    secret_key = keygen->secret_key();
    encryptor = new Encryptor(*context, secret_key);
    batch_encoder = new BatchEncoder(*context);

    if(_verbose) {
        uint64_t coeff_bitcount = 0;
        for (auto& prime : enc_params->coeff_modulus()) {
            coeff_bitcount += prime.bit_count();
        }
        cout << "-------------------------- Crypto Parameters "
                "--------------------------"
             << endl
             << "\tPoly Mod Degree: " << enc_params->poly_modulus_degree() << endl
             << "\tPlain Modulus: " << enc_params->plain_modulus().value() << " ("
             << enc_params->plain_modulus().bit_count() << " bits)" << endl
             << "\tCiphertext Modulus Bitcount: " << coeff_bitcount << endl;
    }
}

void Client::send_parameters_and_keys(bool debug_mode) {
    metrics->metrics_["comm_params"] = enc_params->save(params_stream);

    // If require custom rotation keys, uncomment below
    vector<uint32_t> elts;
    for(size_t i = 2; i <= enc_params->poly_modulus_degree(); i *= 2) {
        elts.push_back(i + 1);
    }

    Serializable<RelinKeys> rlk_client = keygen->create_relin_keys();
    metrics->metrics_["comm_relin_keys"] = rlk_client.save(params_stream);
    Serializable<GaloisKeys> gal_keys_client = keygen->create_galois_keys(elts);
    metrics->metrics_["comm_gal_keys"] = gal_keys_client.save(params_stream);
    Serializable<PublicKey> public_key = keygen->create_public_key();
    metrics->metrics_["comm_pk"] = public_key.save(params_stream);
    
    // If debugging
    if(debug_mode) secret_key.save(params_stream);
}

void Client::encrypt_and_send(const vector<Plaintext> &plain_input, bool _verbose) {
    data_stream.str("");
    metrics->metrics_["comm_request"] = 0;
    for(auto &p: plain_input) {
        metrics->metrics_["comm_request"] += encryptor->encrypt_symmetric(p).save(data_stream);
    }
}

vector<Plaintext> Client::encode_attr_vec(const vector<uint64_t> &attr_vec, size_t num_ctxt_per_value) {
    vector<Plaintext> ret;

    auto encode = [](uint64_t n) {
        return (n == 0 ? string("1") : "1x^" + to_string(n));
    };

    if(num_ctxt_per_value == 1) {
        for(auto &n: attr_vec) {
            ret.emplace_back(encode(n));
        }
    } else if(num_ctxt_per_value == 2) {
        for(auto &n: attr_vec) {
            ret.emplace_back(encode(n / DomainExtensionModulo));
            ret.emplace_back(encode(n % DomainExtensionModulo));
        }
    } else if(num_ctxt_per_value == 3) {
        for(auto &n: attr_vec) {
            ret.emplace_back(encode((n / DomainExtensionModulo) / DomainExtensionModulo));
            ret.emplace_back(encode((n / DomainExtensionModulo) % DomainExtensionModulo));
            ret.emplace_back(encode(n % DomainExtensionModulo));
        }   
    }

    return ret;
}

vector<Plaintext> Client::load_and_decrypt(size_t output_size, bool _verbose) {
    std::vector<Plaintext> result;
    Decryptor decryptor(*context, secret_key);
    Ciphertext __temp_ct;
    Plaintext __temp_pt(enc_params->poly_modulus_degree());
    metrics->metrics_["comm_response"] = 0;
    for(size_t i = 0; i < output_size; i++) {
        metrics->metrics_["comm_response"] += __temp_ct.load(*(context), data_stream);
        decryptor.decrypt(__temp_ct, __temp_pt);
        result.push_back(__temp_pt);
    }
    return result;
}

uint64_t Client::extract_response(const vector<Plaintext> &response, bool mult_path) {
    if(mult_path) {
        return response[0].data()[0] - 1;
    }

    auto labels = response[0], indicators = response[1];

    if(labels.coeff_count() > indicators.coeff_count()) {
        return labels.data()[indicators.coeff_count()] - 1;
    }

    for(size_t i = 0;;i++) {
        if(indicators.data()[i] == 0) {
            return labels.data()[i] - 1;
        }
    }
}

// End to End
bool Client::run_protocol(const vector<uint64_t> &attr_vec, const Config &config) {
    Timer time_setup_crypto, time_send_input, time_server_latency, time_extract_response;
    metrics = new Metrics();

    bool _verbose = config.verbose;
    bool debug_mode = config.debug_mode;
    bool mult_path = config.mult_path;
    string write_to_file = config.write_to_file;
    
    size_t num_ctxt_per_value = 0;
    if(config.bitlength <= 12) {
        num_ctxt_per_value = 1;
    } else if(config.bitlength <= 24) {
        num_ctxt_per_value = 2;
    } else if(config.bitlength <= 36) {
        num_ctxt_per_value = 3;
    } else {
        throw "bitlength too large";
    }

    metrics->metrics_["bitlength"] = config.bitlength;

    if(_verbose) {
        cout << string(88, '<') << endl;
    }

    time_setup_crypto.start();
    if(mult_path) {
        setup_crypto(14, 20, _verbose);
    } else if(num_ctxt_per_value == 1) {
        setup_crypto(12, 16, _verbose);
    } else {
        setup_crypto(13, 18, _verbose);
    }
    metrics->metrics_["time_setup_crypto"] = time_setup_crypto.end_and_get();

    // Sending Parameters
    send_parameters_and_keys(debug_mode);
    server->load_parameters_and_keys(params_stream, debug_mode, _verbose);

    // Sending Input
    // GENERATE INPUT HERE
    time_send_input.start();
    metrics->metrics_["num_attr"]=attr_vec.size();
    vector<Plaintext> encoded = encode_attr_vec(attr_vec, num_ctxt_per_value);
    encrypt_and_send(encoded, _verbose);
    metrics->metrics_["time_send_input"] = time_send_input.end_and_get();

    // Invoking server
    time_server_latency.start();
    server->do_server_computation(data_stream, encoded.size(), metrics, num_ctxt_per_value, mult_path, _verbose);
    metrics->metrics_["time_server_latency"] = time_server_latency.end_and_get();

    // Extracting Response
    time_extract_response.start();
    {
        size_t output_size = (mult_path ? 1 : 2); // Specify how many cts are sent over the network
        std::vector<Plaintext> _response_pts = load_and_decrypt(output_size, _verbose);
        uint64_t result = extract_response(_response_pts, mult_path);
        uint64_t correct_answer = server->evaluate_plain(attr_vec);
        if(result == correct_answer) {
            cout << "-- PASS --" << endl;
        } else {
            cout << "-- FAIL --" << endl;
            cout << "Expected: " << correct_answer << " , got: " << result << endl;
        }
        metrics->metrics_["correctness"] = (result == correct_answer ? 1 : 0);
    }
    metrics->metrics_["time_extract_response"] = time_extract_response.end_and_get();

    // add some more parameters
    metrics->metrics_["num_of_threads"] = thread::hardware_concurrency();
    metrics->metrics_["poly_modulus_degree"] = enc_params->poly_modulus_degree();
    metrics->metrics_["plain_modulus"] = enc_params->plain_modulus().value();
    metrics->metrics_["mult_path"] = config.mult_path;

    // Printing timings
    if(_verbose) {
        cout << "----------------------------------- Timing "
                "-----------------------------------"
             << endl;
        cout << "\tTotal Server    : " << setw(10)
             << metrics->metrics_["time_server"] << " ms" << endl;
             
        cout << "------------------------------- Communication "
                "--------------------------------"
             << endl << "\tData Independant: "
             << metrics->metrics_["comm_relin_keys"] / 1000000 << " MB (Relin keys) + "
             << metrics->metrics_["comm_gal_keys"] / 1000000 << " MB (Gal Keys) + " 
             << metrics->metrics_["comm_pk"] / 1000 << " KB (Public Keys)" 
             << endl << "\tData Dependant: "
             << metrics->metrics_["comm_request"] / 1000 << " KB (Query) + "
             << metrics->metrics_["comm_response"] / 1000 << " KB (Reponse)" << endl
             << endl
             << string(88, '<')
             << endl;
    }

    uint64_t random_num = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count() % 100000000000;

    string fname = write_to_file + to_string(random_num) + ".csv";
    cout << "Writing to file " << fname << endl;

    // rely on the user to ensure that results directory exists
    ofstream ofs;
    ofs.open(fname);
    for(auto &metric: metrics->metrics_) {
        ofs << metric.first << ", " << metric.second << "\n";
    }
    ofs.close();

    return true;
}

