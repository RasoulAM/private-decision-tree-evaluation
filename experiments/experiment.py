from asyncio import subprocess
from re import L, sub
from sys import stdout
import yaml
import subprocess
import os
from generate_random_tree import *
import time

filename='/root/private-decision-tree-evaluation/experiments/experiment.yaml'
ROOT_PATH='/root/private-decision-tree-evaluation'

with open(filename) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)

    bitlength_list=data['main_experiment']['bitlength']
    max_depth_list=data['main_experiment']['max_depth']
    num_attributes=4

    for max_depth in max_depth_list:
        for bitlength in bitlength_list:
            tree_path=os.path.join(ROOT_PATH, "experiments/data/temp_tree.txt")
            input_path=os.path.join(ROOT_PATH, "experiments/data/input.txt")
            generate_balanced_tree(max_depth=max_depth, path=tree_path, bitlength=bitlength, num_attributes=num_attributes, seed=42)
            generate_input(path=input_path, bitlength=bitlength, num_attributes=num_attributes, seed=41)

            cmd=f'./main '\
                ' -f' \
                f' -t {tree_path}' \
                f' -v {input_path}' \
                ' -w 16 -c -s' \
                f' -n {bitlength}'
            
            res = subprocess.run(
                cmd, shell=True,
                cwd='/root/private-decision-tree-evaluation/build/',
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            if "CORRECT" in res.stdout.decode("utf-8"):
                print("PASS")
            else:
                print("FAIL", cmd)
            with open("log.txt", "a") as log_file:
                log_file.write(cmd+"\n")
                log_file.write(res.stdout.decode("utf-8"))
                log_file.write(res.stderr.decode("utf-8"))