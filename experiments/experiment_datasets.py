from asyncio import subprocess
from re import L, sub
from sys import stdout
import yaml
import subprocess
import os
from generate_random_tree import *
import time

from commons import *

experiment_name = 'dataset_experiment'

def valid_setting(dataset, bitlength, hw):
    if hw != 0 and bitlength/hw > 10:
        return False
    if hw != 0 and bitlength/hw < 2:
        return False
    return True

with open(os.path.join(PROJECT_ROOT,'experiments/experiment.yaml')) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)

    dataset_name_list=data[experiment_name]['dataset']
    bitlength_list=data[experiment_name]['bitlength']
    hamming_weight_list=data[experiment_name]['hamming_weight']
    repeat=data[experiment_name]['repeat']

    for rep in range(repeat):
        print(rep)
        for dataset_name in dataset_name_list:
            data_path = f'experiments/datasets_quantized/{dataset_name}'
            write_path=os.path.join(PROJECT_ROOT, f'experiments/results-{version_tag}/src1', dataset_name)
            if not os.path.exists(write_path):
                os.makedirs(write_path)
            else:
                print(f"Directory {write_path} exists")

            with open(os.path.join(PROJECT_ROOT, data_path, f"{dataset_name}.config"), "r") as f:
                num_attributes=int(f.read())
                
            for bitlength in bitlength_list:
                for hamming_weight in hamming_weight_list:

                    if not valid_setting(dataset_name, bitlength, hamming_weight):
                        continue

                    tree_path=os.path.join(PROJECT_ROOT, os.path.join(data_path, f"tree_{dataset_name}_n_{bitlength}.json"))
                    input_path=os.path.join(PROJECT_ROOT, os.path.join(data_path, f"input_{dataset_name}_n_{bitlength}.txt"))
                    generate_input(path=input_path, bitlength=bitlength, num_attributes=num_attributes)

                    cmd=f'./main '\
                        f' -t {tree_path}' \
                        f' -v {input_path}' \
                        f' -s {write_path}/' \
                        f' -n {bitlength}' \
                        f' -w {hamming_weight}' \
                        f' -e {1 if hamming_weight==0 else 0}'
                    
                    res = subprocess.run(
                        cmd, shell=True,
                        cwd=os.path.join(PROJECT_ROOT,'src/build/'),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    if "CORRECT" in res.stdout.decode("utf-8"):
                        pass
                        # print("PASS")
                    else:
                        # print("FAIL", cmd)
                        with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/fail-src1.txt"), "a") as log_file:
                            log_file.write(cmd+"\n")
                    with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/log-src1.txt"), "a") as log_file:
                        log_file.write(cmd+"\n")
                        log_file.write(res.stdout.decode("utf-8"))
                        log_file.write(res.stderr.decode("utf-8"))