from asyncio import subprocess
from re import L, sub
from sys import stdout
import yaml
import subprocess
import os
import time

from commons import *

experiment_name = 'dataset_experiment_sortinghats'
data_path = 'experiments/datasets_quantized'

with open(os.path.join(PROJECT_ROOT,'experiments/experiment.yaml')) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)

dataset_name_list = data[experiment_name]['dataset']
repeat = data[experiment_name]['repeat']


for i in range(repeat):
    for dataset_name in dataset_name_list:   
        write_path=os.path.join(PROJECT_ROOT, f'experiments/results-{version_tag}/sortinghats', dataset_name)
        if not os.path.exists(write_path):
            os.makedirs(write_path)
            
        with open(os.path.join(PROJECT_ROOT, data_path, dataset_name, f"{dataset_name}.config"), "r") as f:
            num_attributes=int(f.read())

        cmd =   '/home/r5akhava/SortingHat/src/rust_pdte/target/release/homdte'\
                f' --dir={dataset_name}'\
                ' --input-size=1'\
                ' --artificial'\
                ' --parallel'
        
        res = subprocess.run(
            cmd, shell=True,
            cwd=os.path.join(PROJECT_ROOT, data_path),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        with open(os.path.join(write_path, f"{i}.csv"), "w") as out_file:
            out_file.write("dataset,depth,leaf_count,internal_count,class_count,duration\n")
            out_file.write(res.stdout.decode("utf-8"))

        with open(os.path.join(PROJECT_ROOT, f"experiments/results-{version_tag}/log-sortinghats.txt"), "a") as log_file:
            log_file.write(cmd+"\n")
            log_file.write(res.stdout.decode("utf-8"))
            log_file.write(res.stderr.decode("utf-8"))
