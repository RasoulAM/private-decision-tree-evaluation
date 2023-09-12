from asyncio import subprocess
from re import L, sub
from sys import stdout
import yaml
import subprocess
import os
from generate_random_tree import *
import time

from commons import *

experiment_name = 'dataset_experiment_xcmp'
data_path = 'experiments/datasets_quantized'

with open(os.path.join(PROJECT_ROOT,'experiments/experiment.yaml')) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)

dataset_name_list = data[experiment_name]['dataset']
bitlength_list = data[experiment_name]['bitlength']
mult_path_list = data[experiment_name]['mult_path']
repeat = data[experiment_name]['repeat']

def is_valid_config(dataset, bitlength, mult_path):
    if mult_path and bitlength > 16:
        return False
    if mult_path and dataset=="spam" and bitlength > 8:
        return False
    return True

for _ in range(repeat):
    for mult_path in mult_path_list:
        for bitlength in bitlength_list:
            print(f"Running tests with bit length {bitlength}...")
            for dataset_name in dataset_name_list:
                if not is_valid_config(dataset_name, bitlength, mult_path):
                    continue
                data_path = f'experiments/datasets_quantized/{dataset_name}'
                write_path=os.path.join(PROJECT_ROOT, f'experiments/results-{version_tag}/src2', dataset_name)
                if not os.path.exists(write_path):
                    os.makedirs(write_path)
                # else:
                #     print(f"Directory {write_path} exists")

                with open(os.path.join(PROJECT_ROOT, data_path, f"{dataset_name}.config"), "r") as f:
                    num_attributes=int(f.read())

                tree_path=os.path.join(PROJECT_ROOT, os.path.join(data_path, f"tree_{dataset_name}_n_{bitlength}.json"))
                input_path=os.path.join(PROJECT_ROOT, os.path.join(data_path, f"input_{dataset_name}_n_{bitlength}.txt"))
                generate_input(path=input_path, bitlength=bitlength, num_attributes=num_attributes)

                cmd =   f'./main ' \
                        f' -m {tree_path}' \
                        f' -a {input_path}' \
                        f' -n {bitlength}' \
                        f' -s {write_path}/'

                if mult_path:
                    cmd += ' -x'
                
                res = subprocess.run(
                    cmd, shell=True,
                    cwd=os.path.join(PROJECT_ROOT,'src2/build/'),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )

                if "PASS" in res.stdout.decode("utf-8"):
                    print("PASS")
                else:
                    print("FAIL", cmd)
                    exit(0)

                with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/log-src2.txt"), "a") as log_file:
                    log_file.write(cmd+"\n")
                    log_file.write(res.stdout.decode("utf-8"))
                    log_file.write(res.stderr.decode("utf-8"))
