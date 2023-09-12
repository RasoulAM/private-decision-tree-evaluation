# Private Decision Tree Evaluation

This is the artifact for the paper: "Level Up: Private Non-Interactive Decision Tree Evaluation using Levelled Homomorphic Encryption"

## Requirements
To build the project, you need the following:
- CMake 3.10 or higher
- C++ compiler with C++17 support
- Microsoft SEAL 

## Building the project
To build the project
```
    cd src
    mkdir build
    cd build
    cmake ..
    make
    cmake ../../
    cd src2
    mkdir build
    cd build
    cmake ..
    make
```

## Example to run

Run in the build directory
```
    ./main -t <path-to-tree> -v <path-to-input> -n <bitlength>
```

## How to run experiments
Assumes the binary is available in ```src/build```.

In ```experiments/``` directory, run
```
    cd experiments
    ./prepare_data # Prepare data
    python3 run_all_experiments.py # Run all experiments

```


