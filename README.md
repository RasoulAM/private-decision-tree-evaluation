# private-decision-tree-evaluation


## Format of Tree
Traverse the tree in preorder format. Each line is one node.

If it's a non-leaf node:
```
<threshold-value> <attribute-index>
```

If it's a leaf node:
```
<classification-value> -1
```

## Format of the input
Attribute $i$ in on line $i$

## Building the project
To build the project
```
    cd src
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