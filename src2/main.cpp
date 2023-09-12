#include "utils.h"
#include "client.h"
#include "server.h"

#include <iostream>
#include <unistd.h>
#include <omp.h>

using namespace std;

// example: ./main -v -m tree.in -a attr_vec.in -n 14
int main(int argc, char *argv[]) {

    omp_set_num_threads(32);

    string model_fname;
    vector<uint64_t> attr_vec;
    Config config;
    string err_msg;

    int opt;
    while((opt = getopt(argc, argv, "dvs:m:a:n:")) != -1) {
        switch(opt) {
            case 'd': {
                config.debug_mode = true;
                break;
            }
            case 'v': {
                config.verbose = true;
                break;
            }
            case 's': {
                config.write_to_file = string(optarg);
                break;
            }
            case 'm': {
                model_fname = string(optarg);
                break;
            }
            case 'a': {
                read_attr_vec(optarg, attr_vec);
                break;
            }
            case 'n': {
                config.bitlength = stoi(string(optarg));
                break;
            }
        }
    }

    if(model_fname.empty()) {
        err_msg = "Empty classification model.";
    } else if(attr_vec.empty()) {
        err_msg = "Empty attribute vector.";
    } else if(config.bitlength < 1 || config.bitlength > 36) {
        err_msg = "Invalid bitlength, or you didn't specify one.";
    }

    if(!err_msg.empty()) {
        cerr << err_msg << endl;
        return 0;
    }

    auto server = Server();
    server.fetch_model(model_fname);

    Client client(&server);
    client.run_protocol(attr_vec, config);
}
