import argparse
from ast import parse
import string
import random
from math import comb

CLASSIFICATION_VALUE_BITLENGTH=8


class Leaf(object):
    def __init__(self, v):
        self.type = "leaf"
        self.value = int(v)

class Internal(object):
    def __init__(self, threshold, feature, left, right):
        self.type = "split"
        self.threshold = int(threshold)
        self.feature = int(feature)
        self.left = left.__dict__
        self.right = right.__dict__


def generate_and_write_tree(preorder, args):
    bitlength=args.bitlength
    num_attributes = args.num_attributes
    path = args.path

    with open(path, 'w') as output_file:
        for is_leaf in preorder:
            if is_leaf:
                threshold = random.randint(0, 2**CLASSIFICATION_VALUE_BITLENGTH-1)
                attribute_index = -1
                output_file.write('{} {}\n'.format(threshold, attribute_index))
            else:
                threshold = random.randint(0, 2**bitlength-1)
                attribute_index = random.randint(0, num_attributes-1)
                output_file.write('{} {}\n'.format(threshold, attribute_index))

def generate_balanced_tree_rec(max_depth, bitlength, num_attributes, seed=None):

    if max_depth>0:
        left = generate_balanced_tree_rec(max_depth-1, bitlength, num_attributes, seed)
        right = generate_balanced_tree_rec(max_depth-1, bitlength, num_attributes, seed)
        threshold = random.randint(0, 2**bitlength-1)
        feature = random.randint(0, num_attributes-1)
        t = Internal(threshold, feature, left, right)
    else:
        t = Leaf(random.randint(0, 2**CLASSIFICATION_VALUE_BITLENGTH-1))

    return t    


def generate_balanced_tree(max_depth, path, bitlength, num_attributes, seed=None):
    if seed is not None:
        random.seed(seed)
    with open(path, 'w') as output_file:
        stack=[0]
        while len(stack) > 0:
            current_depth = stack[0]
            if current_depth == max_depth:
                threshold = random.randint(0, 2**CLASSIFICATION_VALUE_BITLENGTH-1)
                attribute_index = -1
                output_file.write(f'{threshold} {attribute_index}\n')
                
                stack = stack[1:]
            elif current_depth < max_depth:
                threshold = random.randint(0, 2**bitlength-1)
                attribute_index = random.randint(0, num_attributes-1)
                output_file.write(f'{threshold} {attribute_index}\n')
                stack = [current_depth+1, current_depth+1] + stack[1:]
            else:
                print("This should not happen!")

def generate_input(path, bitlength, num_attributes, seed=None):
    if seed is not None:
        random.seed(seed)
    with open(path, 'w') as output_file:
        for _ in range(num_attributes):
            attribute_value = random.randint(0, 2**bitlength-1)
            output_file.write(f'{attribute_value}\n')
    

def generate_balanced_tree_from_args(args):
    max_depth = args.max_depth
    path = args.path
    bitlength=args.bitlength
    num_attributes = args.num_attributes
    generate_balanced_tree(max_depth=max_depth,path=path,bitlength=bitlength,num_attributes=num_attributes)


def catalan(n):
    return comb(2*n, n) / (n+1)

def recursive_uniform_tree_generator(remaining_inner_nodes):
    if remaining_inner_nodes == 0:
        return [True]

    remaining_inner_nodes -= 1
    weights = [ catalan(i)*catalan(remaining_inner_nodes-i) for i in range(0,remaining_inner_nodes+1) ]
    
    left_inner_nodes=random.choices(range(0,remaining_inner_nodes+1), weights=weights)[0]
    right_inner_nodes=remaining_inner_nodes-left_inner_nodes

    return [False] + recursive_uniform_tree_generator(left_inner_nodes) + recursive_uniform_tree_generator(right_inner_nodes)

def generate_random_unbalanced_tree(args):
    number_of_nodes=args.number_of_nodes
    preorder_representation=recursive_uniform_tree_generator(number_of_nodes)
    generate_and_write_tree(preorder_representation, args)

# if __name__ == '__main__':

#     parser = argparse.ArgumentParser()

#     # common for both cases    
#     parser.add_argument('--bitlength', type=int, default=32)
#     parser.add_argument('--num_attributes', type=int, default=4)
#     parser.add_argument('--balanced', const=True, default=False, nargs='?')

#     # balanced case
#     parser.add_argument('--max_depth', type=int, default=4)
#     parser.add_argument('--path', type=str, required=True)

#     # unbalanced case
#     parser.add_argument('--number_of_nodes', type=int, default=31)

#     args = parser.parse_args()

#     if args.balanced:
#         print("Balanced!")
#         generate_balanced_tree_from_args(args)
#     else:
#         print("Not Balanced!")
#         generate_random_unbalanced_tree(args)
        

import os 
import json

WORKSPACE_DIR='/home/r5akhava/private-decision-tree-evaluation/experiments'

if __name__ == '__main__':
    write_path = 'datasets_synthetic'
    for max_depth in range(2,21):
        for bitlength in [8, 12, 16, 24, 32]:
            for num_attributes in range(2, 200, 10):
                t = generate_balanced_tree_rec(max_depth, bitlength, num_attributes)
                with open(os.path.join(WORKSPACE_DIR, write_path, f'tree_depth_{max_depth}_n_{bitlength}_attr_{num_attributes}.json'), 'w+') as f:
                    f.write(json.dumps(t.__dict__))

