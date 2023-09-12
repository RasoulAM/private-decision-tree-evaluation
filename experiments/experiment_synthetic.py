from asyncio import subprocess
from re import L, sub
from sys import stdout
import yaml
import subprocess
import os
from generate_random_tree import *
import time
from tqdm import tqdm
import random

from commons import *


def valid_setting(max_depth, num_attributes, bitlength, hw):
    if hw != 0 and bitlength/hw > 10:
        return False
    return True

def get_hamming_weight(bitlength):
    if bitlength <= 24:
        return 4
    if bitlength <= 30:
        return 8
    else:
        return 16

with open(os.path.join(PROJECT_ROOT,'experiments/experiment.yaml')) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)

    for experiment_name in ['synthetic_experiment_num_nodes','synthetic_experiment_num_attr']:
        bitlength_list=data[experiment_name]['bitlength']
        max_depth_list=data[experiment_name]['max_depth']
        num_attributes_list=data[experiment_name]['num_attributes']
        repeat=data[experiment_name]['repeat']

        for rep in tqdm(range(repeat)):
            for max_depth in max_depth_list:

                data_path = f'experiments/datasets_synthetic'
                write_path=os.path.join(PROJECT_ROOT, f'experiments/results-{version_tag}/src1/synthetic')
                if not os.path.exists(write_path):
                    os.makedirs(write_path)

                write_path_xcmp=os.path.join(PROJECT_ROOT, f'experiments/results-{version_tag}/src2/synthetic')
                if not os.path.exists(write_path_xcmp):
                    os.makedirs(write_path_xcmp)
                    
                write_path_sortinghats=os.path.join(PROJECT_ROOT, f'experiments/results-{version_tag}/sortinghats/synthetic')
                if not os.path.exists(write_path_sortinghats):
                    os.makedirs(write_path_sortinghats)

                bitlength = 11
                for num_attributes in num_attributes_list:
                    
                    path_to_model = os.path.join(PROJECT_ROOT, data_path, f'tree_depth_{max_depth}_attr_{num_attributes}')

                    cmd =   '/home/r5akhava/SortingHat/src/rust_pdte/target/release/homdte'\
                            f' --dir={path_to_model}'\
                            ' --input-size=1'\
                            ' --artificial'\
                            ' --parallel'
                    
                    res = subprocess.run(
                        cmd, shell=True,
                        cwd=os.path.join(PROJECT_ROOT, data_path),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    
                    # get the time in millisconds as a 8 digit number
                    id_num = int(time.time())

                    with open(os.path.join(write_path_sortinghats, f"{id_num}.csv"), "w") as out_file:
                        out_file.write("num_attributes,dataset,depth,leaf_count,internal_count,class_count,duration\n")
                        out_file.write(str(num_attributes)+","+res.stdout.decode("utf-8"))

                    with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/log-synthetic-sortinghats.txt"), "a") as log_file:
                        log_file.write(cmd+"\n")
                        log_file.write(res.stdout.decode("utf-8"))
                        log_file.write(res.stderr.decode("utf-8"))


                for bitlength in bitlength_list:
                    for num_attributes in num_attributes_list:
                            
                        hamming_weight=get_hamming_weight(bitlength)

                        if not valid_setting(max_depth, num_attributes, bitlength, hamming_weight):
                            continue

                        tree_path=os.path.join(PROJECT_ROOT, data_path, f"tree_depth_{max_depth}_n_{bitlength}_attr_{num_attributes}.json")
                        input_path=os.path.join(PROJECT_ROOT, data_path, f"input_n_{bitlength}_attr_{num_attributes}.txt")
                        generate_input(path=input_path, bitlength=bitlength, num_attributes=num_attributes)

                        # src1

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
                            with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/fail-synthetic-src1.txt"), "a") as log_file:
                                log_file.write(cmd+'\n')
                        with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/log-synthetic-src1.txt"), "a") as log_file:
                            log_file.write(cmd+"\n")
                            log_file.write(res.stdout.decode("utf-8"))
                            log_file.write(res.stderr.decode("utf-8"))

                        # src2

                        cmd =   f'./main ' \
                                f' -m {tree_path}' \
                                f' -a {input_path}' \
                                f' -n {bitlength}' \
                                f' -s {write_path_xcmp}/'
                        
                        res = subprocess.run(
                            cmd, shell=True,
                            cwd=os.path.join(PROJECT_ROOT,'src2/build/'),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE
                        )

                        if "PASS" in res.stdout.decode("utf-8"):
                            pass
                            # print("PASS")
                        else:
                            # print("FAIL", cmd)
                            # exit(0)
                            with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/fail-synthetic-src2.txt"), "a") as log_file:
                                log_file.write(cmd+'\n')
                        with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/log-synthetic-src2.txt"), "a") as log_file:
                            log_file.write(cmd+"\n")
                            log_file.write(res.stdout.decode("utf-8"))
                            log_file.write(res.stderr.decode("utf-8"))
