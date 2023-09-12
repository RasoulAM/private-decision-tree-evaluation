import yaml
from commons import *
import subprocess

def run_command(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        print(f"Command failed with error: {stderr.decode()}")
        return False
    else:
        return True

def main():
    commands1 = [
        "python3 experiment_datasets.py",
        "python3 experiment_datasets_xcmp.py",
        "python3 experiment_datasets_sortinghats.py",
        "python3 experiment_synthetic.py",
        f"mkdir results-{version_tag}/cmp/ && " \
         "../src/build/cmp_bench 10 > results-{version_tag}/cmp/cmp_bench.csv",
    ]

    commands2 = [
        "python3 visualization.py",
        "python3 visualization_synthetic.py",
        "python3 visualization_cmp.py",
    ]

    for command in commands1:
        if not run_command(command):
            return
    print("All experiments done!")
    
    for command in commands2:
        if not run_command(command):
            return
    print("All visualizations done!")

if __name__ == "__main__":
    main()
