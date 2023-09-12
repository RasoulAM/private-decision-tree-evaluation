import os

import matplotlib.pyplot as plt
import pandas as pd
import yaml
from commons import *

num_attributes={
    "breast": 30,
    "steel": 33,
    "heart": 13,
    "spam": 57,
}

def get_points_list(x, y):
    assert len(x.values) == len(y.values)
    return [(x.values[i], y.values[i]) for i in range(len(x.values))]

def read_directory(dir_path):
    res=[]
    for path in os.listdir(dir_path):
        # check if current path is a file
        if os.path.isfile(os.path.join(dir_path, path)):
            x = pd.read_csv(os.path.join(dir_path, path), names=['metric', 'value'], index_col=0).T
            res.append(x)
    if len(res) > 0:
        return pd.concat(res)
    else:
        return None    
        
def read_sortinghats(dir_path):
    res=[]
    for path in os.listdir(dir_path):
        # check if current path is a file
        if os.path.isfile(os.path.join(dir_path, path)):
            x = pd.read_csv(os.path.join(dir_path, path))
            res.append(x)
    if len(res) > 0:
        return pd.concat(res)
    else:
        return None   

with open(os.path.join(PROJECT_ROOT,'experiments/experiment.yaml')) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)

    dataset_name_list=data['plotting']['dataset']
    results_path=f'experiments/results-{version_tag}'

    # Create a dictionary with a list of upper and lower bounds for each dataset
    time_bounds = {
        'breast': [-10, 700],
        'steel': [-10, 250],
        'heart': None,
        'spam': None,
    }
    
    comm_bounds = {
        'breast': [800, 20000],
        'steel': [900, 20000],
        'heart': [300, 20000],
        'spam': [1000, 26000],        
    }
    
    
    for dataset_name in dataset_name_list:
        print(dataset_name)

        f1, plt_time=plt.subplots()
        f2, plt_comm=plt.subplots()
        
        # Get the current figure size in inches
        current_fig_size = plt.gcf().get_size_inches()

        # Set the desired aspect ratio: height = 0.5 * width
        desired_width = current_fig_size[0]
        desired_height = desired_width / 1.618

        # Update the figure size with the new dimensions
        f1.set_size_inches(desired_width, desired_height)   
        f2.set_size_inches(desired_width, desired_height)   
             
        # src1

        raw_data=read_directory(os.path.join(PROJECT_ROOT, results_path, 'src1', dataset_name))
        
        raw_data=raw_data[(2 <= raw_data["bitlength"]) & (raw_data["bitlength"]<=16)]
        
        if raw_data is not None:
            # Keep correct results, average them
            raw_data=raw_data[raw_data["correctness"]==1]
            data_avg=raw_data.groupby(['bitlength', 'hamming_weight', 'comparison'], as_index=False).mean()
            data_std=raw_data.groupby(['bitlength', 'hamming_weight', 'comparison'], as_index=False).std()

            # Folklore
            folklore_data_avg=data_avg[data_avg["comparison"]==1]
            folklore_data_std=data_std[data_std["comparison"]==1]
            
            plt_time.fill_between(
                folklore_data_avg['bitlength'], 
                folklore_data_avg['time_server_crypto'] - folklore_data_std['time_server_crypto'],
                folklore_data_avg['time_server_crypto'] + folklore_data_std['time_server_crypto'],
                alpha=0.2
            )
            plt_time.plot(folklore_data_avg['bitlength'], folklore_data_avg['time_server_crypto'], label=f'Folklore-PDTE')
            plt_comm.plot(folklore_data_avg['bitlength'], folklore_data_avg['comm_query']/1000, label=f'Folklore-PDTE')

            # Range cover
            for hamming_weight in [2,4]:
                filtered_avg=data_avg[(data_avg["hamming_weight"]==hamming_weight) & (data_avg["comparison"]==0)]
                filtered_std=data_std[(data_std["hamming_weight"]==hamming_weight) & (data_std["comparison"]==0)]
                
                plt_time.fill_between(
                    filtered_avg['bitlength'],
                    filtered_avg['time_server_crypto'] - filtered_std['time_server_crypto'], 
                    filtered_avg['time_server_crypto'] + filtered_std['time_server_crypto'],
                    alpha=0.2
                )
                plt_time.plot(filtered_avg['bitlength'], filtered_avg['time_server_crypto'], label=f'RCC-PDTE (h={hamming_weight})')
                plt_comm.plot(filtered_avg['bitlength'], filtered_avg['comm_query']/1000, label=f'RCC-PDTE (h={hamming_weight})')


        #src2

        raw_data_xcmp=read_directory(os.path.join(PROJECT_ROOT, results_path, 'src2', dataset_name))
        raw_data_xcmp=raw_data_xcmp[(2 <= raw_data_xcmp["bitlength"] )&( raw_data_xcmp["bitlength"]<=16)]
        
        
        if raw_data_xcmp is not None:

            # Only valid entries
            raw_data_xcmp = raw_data_xcmp[raw_data_xcmp['correctness']==1]
            data_xcmp_avg=raw_data_xcmp.groupby(['bitlength', 'mult_path'], as_index=False).mean()
            data_xcmp_std=raw_data_xcmp.groupby(['bitlength', 'mult_path'], as_index=False).std()

            # Sum Path
            filtered=data_xcmp_avg[data_xcmp_avg['mult_path']==0]
            filtered_std=data_xcmp_std[data_xcmp_std['mult_path']==0]
            plt_time.fill_between(
                filtered['bitlength'],
                filtered['time_server'] - filtered_std['time_server'],
                filtered['time_server'] + filtered_std['time_server'],
                alpha=0.2
            )
            plt_time.plot(filtered['bitlength'], filtered['time_server'], label=f'XXCMP-PDTE')
            plt_comm.plot(filtered['bitlength'], filtered['comm_request']/1000, label=f'XXCMP-PDTE')

            # Mult Path
            # filtered=data_xcmp_avg[data_xcmp_avg['mult_path']==1]
            # plt_time.fill_between(
            #     filtered['bitlength'],
            #     filtered['time_server'] - filtered_std['time_server'],
            #     filtered['time_server'] + filtered_std['time_server'],
            #     alpha=0.2
            # )
            # plt_time.plot(filtered['bitlength'], filtered['time_server'], label=f'XCMP (mult)')
            # plt_comm.plot(filtered['bitlength'], filtered['comm_request']/1000, label=f'XCMP (mult)')

        # SortingHats

        raw_data_sortinghats=read_sortinghats(os.path.join(PROJECT_ROOT, results_path, 'sortinghats', dataset_name))
        
        if raw_data_sortinghats is not None:
            # Add approximated communication
            raw_data_sortinghats['comm'] = raw_data_sortinghats['dataset'].apply(lambda x: num_attributes[x])*64*2048*2/8000
            data_sortinghats_avg=raw_data_sortinghats.groupby(['dataset'], as_index=False).mean()
            data_sortinghats_std=raw_data_sortinghats.groupby(['dataset'], as_index=False).std()
            
            # Rewrite the line below, but make the dot black
            plt_time.plot(11, data_sortinghats_avg['duration'], 'o', color='black', label=f'SortingHats')
            plt_comm.plot(11, data_sortinghats_avg['comm'], 'o', color='black', label=f'SortingHats (Approx.)')  

        
        figures_dir=os.path.join(PROJECT_ROOT, results_path, 'figures')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)

        plt_time.set_ylabel("Milliseconds")
        plt_time.set_xlabel("Bit Precision")
        plt_time.legend()
        if time_bounds[dataset_name] is not None:
            # set the yaxis range
            plt_time.set_ylim(time_bounds[dataset_name][0], time_bounds[dataset_name][1])
        plt_time.set_title(f"Inference Time for {dataset_name.capitalize()} dataset ({num_attributes[dataset_name]} attributes)")
        f1.savefig(os.path.join(figures_dir, f"time-{dataset_name}.pdf"), bbox_inches = 'tight')
        plt.close()
        
        plt_comm.set_ylabel("KBytes")
        plt_comm.set_xlabel("Bit Precision")
        plt_comm.legend()
        plt_comm.set_yscale('log')
        if comm_bounds[dataset_name] is not None:
            # set the yaxis range
            plt_comm.set_ylim(comm_bounds[dataset_name][0], comm_bounds[dataset_name][1])
        plt_comm.set_title(f"Communication for {dataset_name.capitalize()} dataset ({num_attributes[dataset_name]} attributes)")
        f2.savefig(os.path.join(figures_dir, f"comm-{dataset_name}.pdf"), bbox_inches = 'tight')
        plt.close()

        # # Comm-Comp Graphs
        # for bitlength in [8, 11, 12, 16, 20 , 24, 36]:
        #     f3, plt_time_comm = plt.subplots()
        #     points = []

        #     selected_data = data_avg[(data_avg['bitlength']==bitlength) & (data_avg['comparison']==1)]
        #     plt_time_comm.plot(selected_data['comm_query'], selected_data['time_server_crypto'], 'o', label="Folklore")
        #     points += get_points_list(selected_data['comm_query'], selected_data['time_server_crypto'])

        #     selected_data = data_avg[(data_avg['bitlength']==bitlength) & (data_avg['comparison']==0)]
        #     for hamming_weight in [4,6,8,16]:
        #         filtered=selected_data[selected_data["hamming_weight"]==hamming_weight]
        #         plt_time_comm.plot(filtered['comm_query'], filtered['time_server_crypto'], 'o', label=f"k={hamming_weight}")
        #         points += get_points_list(filtered['comm_query'], filtered['time_server_crypto'])

        #     selected_data_xcmp = data_xcmp_avg[(data_xcmp_avg['bitlength']==bitlength) & (data_xcmp_avg['mult_path']==0)]
        #     plt_time_comm.plot(selected_data_xcmp['comm_request'], selected_data_xcmp['time_server'], 'o', label="XCMP")
        #     points += get_points_list(selected_data_xcmp['comm_request'], selected_data_xcmp['time_server'])

        #     if bitlength <= 11:
        #         plt_time_comm.plot(data_sortinghats_avg['comm'], data_sortinghats_avg['duration'], 'o', label='SortingHats')
        #         points += get_points_list(data_sortinghats_avg['comm'], data_sortinghats_avg['duration'])

        #     # Throughput lines
        #     # # MB/s -> KB/ms
        #     # for throughput in [100, 1000, 10000, 100000]:
        #     #     pareto_point=min(points, key=lambda x: x[0] + x[1]*throughput)
        #     #     plt_time_comm.axline(pareto_point, slope=-1/throughput, linestyle='--', label=f'Throughput={throughput}')

        #     plt_time_comm.legend()
        #     plt_time_comm.set_title(f"bitlength={bitlength}, Throughputs in MB/s")
        #     # plt_time_comm.set_xscale('log')
        #     plt_time_comm.set_ylabel("Milliseconds")
        #     plt_time_comm.set_xlabel("KBytes")
        #     f3.savefig(os.path.join(figures_dir, f"time_comm/{dataset_name}-n={bitlength}.pdf"))
        
        #     plt.close()
        

