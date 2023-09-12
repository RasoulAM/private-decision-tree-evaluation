#include "client.h"
#include "seal/seal.h"
#include "server.h"
#include "tree_utils.h"
#include "utils.h"
#include "parser.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <vector>

using namespace seal;
using namespace std;

int main(int argc, char* argv[])
{

    if (argc < 3) {
        cout << "Need at least a model and a vector" << endl;
        return -1;
    }

    string model_filename;
    string attribute_vec_filename;

    int opt;
    int bitlength = 32;
    int hamming_weight = 8;
    bool flexible_mode = false;
    PathEvaluation path_eval;
    ComparisonType comparison=RANGE_COVER;
    string write_to_file = "results/"; // The directory to write to (include forward slash)
    while ((opt = getopt(argc, argv, "fs:t:v:w:n:pe:")) != -1) {
        switch (opt) {
        case 's':
            write_to_file = string(optarg);
            break;
        case 't':
            model_filename = string(optarg);
            break;
        case 'v':
            attribute_vec_filename = string(optarg);
            break;
        case 'p':
            path_eval=SUM;
            break;
        case 'w':
            hamming_weight = atoi(optarg);
            break;
        case 'n':
            bitlength = atoi(optarg);
            break;
        case 'e':
            switch (stoi(optarg)){
                case 0:
                    comparison=RANGE_COVER;
                    break;
                case 1:
                    comparison=FOLKLORE;
                    break;
                default:
                    comparison=RANGE_COVER;
            }
            break;
        case ':':
            std::cout << "option needs a value: " << string(optarg) << endl;
            break;
        case '?':
            std::cout << "unknown option: " << optopt << "\n";
            break;
        }
    }
    vector<uint64_t> input = read_input_from_file(attribute_vec_filename);

    QueryParameters* query_parameters = new QueryParameters(bitlength, input.size(), hamming_weight, 32, path_eval, comparison, write_to_file);
    Client* client = new Client();
    client->server->initialize_model(model_filename, query_parameters);
    client->end_to_end_evaluation(query_parameters, input, true);
    delete query_parameters;
}
