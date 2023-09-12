import os

import matplotlib.pyplot as plt
import pandas as pd
import yaml
from commons import *
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

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

    results_path=f'experiments/results-{version_tag}'

    bitlength_color_list = [
        (8,'r'),
        # (12, 'c'),
        (16,'b'),
        (26,'g'),
        # (32,'g'),
    ]

    raw_data=read_directory(os.path.join(PROJECT_ROOT, results_path, f'src1/synthetic'))
    comm_div=10**3
    if raw_data is not None:
        # Keep correct results, average them
        raw_data=raw_data[raw_data["correctness"]==1]
        raw_data['comm_query'] = raw_data['comm_query'] / comm_div
        raw_data['comm_response'] = raw_data['comm_response'] / comm_div
        data_avg=raw_data.groupby(['bitlength', 'num_attr', 'num_internal_nodes'], as_index=False).mean()
        data_std=raw_data.groupby(['bitlength', 'num_attr', 'num_internal_nodes'], as_index=False).std()

    raw_data_xcmp=read_directory(os.path.join(PROJECT_ROOT, results_path, f'src2/synthetic'))
    if raw_data_xcmp is not None:
        # Keep correct results, average them
        raw_data_xcmp = raw_data_xcmp[raw_data_xcmp['correctness']==1]
        raw_data_xcmp['comm_request'] = raw_data_xcmp['comm_request'] / comm_div
        raw_data_xcmp['comm_response'] = raw_data_xcmp['comm_response'] / comm_div
        data_xcmp_avg=raw_data_xcmp.groupby(['bitlength', 'mult_path', 'num_attr', 'inner_nodes'], as_index=False).mean()
        data_xcmp_std=raw_data_xcmp.groupby(['bitlength', 'mult_path', 'num_attr', 'inner_nodes'], as_index=False).std()
    
    raw_data_sortinghats = read_sortinghats(os.path.join(PROJECT_ROOT, results_path, f'sortinghats/synthetic'))
    # Formula for communication: num_attributes * (RLWE CT size = (q=64) * (N=2048) * 2) / 8000
    raw_data_sortinghats['comm'] = raw_data_sortinghats['num_attributes']*64*2048*2/8000
    if raw_data_sortinghats is not None:
        data_sortinghats_avg = raw_data_sortinghats.groupby(['num_attributes', 'depth'], as_index=False).mean()
        data_sortinghats_std = raw_data_sortinghats.groupby(['num_attributes', 'depth'], as_index=False).std()

    f1, plt_time=plt.subplots()
    f2, plt_comm=plt.subplots()
    
    
    ##################################################################################################################################
    
    # create custom legend labels
    rcc_label = mlines.Line2D([], [], color='black', linestyle='-', label='RCC-PDTE')
    xcmp_label = mlines.Line2D([], [], color='black', linestyle='--', label='XXCMP-PDTE')
    sortinghats_label = mlines.Line2D([], [], color='black', linestyle='-.', label='SortingHats')
    n8_handle = mpatches.Rectangle((0,0),1,1, color='red',    label='n=8')
    n16_handle = mpatches.Rectangle((0,0),1,1, color='blue',  label='n=16')
    n26_handle = mpatches.Rectangle((0,0),1,1, color='green', label='n=26')
    
    plt_time.add_artist(plt_time.legend(handles=[rcc_label, xcmp_label, sortinghats_label], loc='upper left'))
    plt_time.add_artist(plt_time.legend(handles=[n8_handle, n16_handle, n26_handle], loc='upper right'))
    
    plt_comm.add_artist(plt_comm.legend(handles=[rcc_label, xcmp_label, sortinghats_label], loc='lower right'))
    plt_comm.add_artist(plt_comm.legend(
        handles=[n8_handle, n16_handle, n26_handle],
        bbox_to_anchor=(0.62, 0),
        loc='lower center'
    ))

    depth=6
    print(r"\toprule")
    print(r"Precision (bits) & RCC Approx. Runtime & XXCMP Approx. Runtime & SortingHats Approx. Runtime \\")
    print(r"\midrule")
    
    for bitlength, color in bitlength_color_list:
        print("{:2.0f} & ".format(bitlength), end="")

        selected=data_avg[(data_avg['bitlength']==bitlength) & (data_avg['num_internal_nodes']==2**depth-1)]
        plt_time.plot(selected['num_attr'], selected['time_server_crypto'], color, label=f"RCC n={bitlength}")
        plt_comm.plot(selected['num_attr'], selected['comm_query'], color, label=f"RCC n={bitlength}")
        
        print("{:0.0f} - {:0.0f}".format(
            selected['time_server_crypto'].mean()-selected['time_server_crypto'].std(),
            selected['time_server_crypto'].mean()+selected['time_server_crypto'].std()
        ), end=" & ")

        selected = data_xcmp_avg[(data_xcmp_avg['bitlength']==bitlength) & (data_xcmp_avg['inner_nodes']==2**depth-1)]
        plt_time.plot(selected['num_attr'], selected['time_server'], color+'--', label=f"XCMP n={bitlength}")
        plt_comm.plot(selected['num_attr'], selected['comm_request'], color+'--', label=f"XCMP n={bitlength}")
        
        print("{:0.0f} - {:0.0f}".format(
            selected['time_server'].mean()-selected['time_server'].std(),
            selected['time_server'].mean()+selected['time_server'].std()
        ), end=" & ")
        
        if bitlength <= 11:
            selected = data_sortinghats_avg[data_sortinghats_avg['depth']==depth]
            plt_time.plot(selected['num_attributes'], selected['duration'], 'r-.', label=f"Sortinghats")
            plt_comm.plot(selected['num_attributes'], selected['comm'], 'r-.', label=f"Sortinghats")

            print("{:0.0f} - {:0.0f}".format(
                selected['duration'].mean()-selected['duration'].std(),
                selected['duration'].mean()+selected['duration'].std()
            ), end="")
            
        else:
            print(" - ", end="")
        
        print(r"\\")
            
    print(r"\bottomrule")
    

    figures_dir=os.path.join(PROJECT_ROOT, results_path, 'figures/sythetic')
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    plt_time.set_ylabel("Milliseconds")
    plt_time.set_xlabel("Number of Attributes")
    # plt_time.legend()
    plt_time.set_title("Inference Time versus Number of Attributes")
    f1.savefig(os.path.join(figures_dir, f"time-depth-{depth}.pdf"), bbox_inches = 'tight')
    plt.close()

    plt_comm.set_ylabel("KBytes")
    plt_comm.set_xlabel("Number of Attributes")
    # plt_comm.legend()
    plt_comm.set_title("Commmunication versus Number of Attributes")
    plt_comm.set_yscale('log')
    f2.savefig(os.path.join(figures_dir, f"comm-depth-{depth}.pdf"), bbox_inches = 'tight')
    plt.close()
    
    ################################################################################################################################## 

    # create custom legend labels
    rcc_label = mlines.Line2D([], [], color='black', linestyle='-', label='RCC-PDTE')
    xcmp_label = mlines.Line2D([], [], color='black', linestyle='--', label='XXCMP-PDTE')
    sortinghats_label = mlines.Line2D([], [], color='black', linestyle='-.', label='SortingHats')
    n8_handle = mpatches.Rectangle((0,0),1,1, color='red', label='n=8')
    n16_handle = mpatches.Rectangle((0,0),1,1, color='blue', label='n=16')
    n26_handle = mpatches.Rectangle((0,0),1,1, color='green', label='n=26')
    
    f1, plt_time=plt.subplots()
    # set the location of the legend manually between 'center left' and 'upper left'
    
    plt_time.add_artist(plt_time.legend(
        handles=[rcc_label, xcmp_label, sortinghats_label],loc='upper left'
    ))
    plt_time.add_artist(plt_time.legend(
        handles=[n8_handle, n16_handle, n26_handle],
        bbox_to_anchor = (0, 0.7),
        loc='center left'
    ))
    
    
    f2, plt_comm=plt.subplots()
    plt_comm.add_artist(plt_comm.legend(handles=[rcc_label, xcmp_label, sortinghats_label], loc='upper left'))
    plt_comm.add_artist(plt_comm.legend(handles=[n8_handle, n16_handle, n26_handle], loc='upper right'))


    num_attr = 32
    print()
    print(r"\toprule")
    print(r"Precision (bits) & RCC-PDTE & XXCMP-PDTE & SortingHats \\")
    print(r"\midrule")
    for bitlength, color in bitlength_color_list:
        print("{:2.0f} & ".format(bitlength), end="")

        selected_data=raw_data[(raw_data['bitlength']==bitlength) & (raw_data['num_attr']==num_attr) & (raw_data['num_internal_nodes']<=2000)]

        data_avg=selected_data.groupby(['num_internal_nodes', 'bitlength', 'num_attr'], as_index=False).mean()
        data_std=selected_data.groupby(['num_internal_nodes', 'bitlength', 'num_attr'], as_index=False).std()

        plt_time.fill_between(
            data_avg['num_internal_nodes'],
            data_avg['time_server_crypto']-data_std['time_server_crypto'],
            data_avg['time_server_crypto']+data_std['time_server_crypto'],
            color=color,
            alpha=0.1
        )
        plt_time.plot(data_avg['num_internal_nodes'], data_avg['time_server_crypto'], color, label=f"RCC n={bitlength}")
        plt_comm.plot(data_avg['num_internal_nodes'], data_avg['comm_query'], color, label=f"RCC n={bitlength}")
        
        print("{:0.0f}".format(        
            data_avg['comm_query'].mean()
            # data_avg['comm_query'].mean()-data_avg['comm_query'].std(),
            # data_avg['comm_query'].mean()+data_avg['comm_query'].std()
        ), end=" & ")

        #### 
        selected_data=raw_data_xcmp[(raw_data_xcmp['bitlength']==bitlength) & (raw_data_xcmp['num_attr']==num_attr) & (raw_data_xcmp['inner_nodes']<=2000)]

        data_xcmp_avg=selected_data.groupby(['inner_nodes', 'bitlength', 'mult_path', 'num_attr'], as_index=False).mean()
        data_xcmp_std=selected_data.groupby(['inner_nodes', 'bitlength', 'mult_path', 'num_attr'], as_index=False).std()

        plt_time.fill_between(
            data_xcmp_avg['inner_nodes'],
            data_xcmp_avg['time_server']-data_xcmp_std['time_server'],
            data_xcmp_avg['time_server']+data_xcmp_std['time_server'],
            color=color,
            alpha=0.1
        )
        plt_time.plot(data_xcmp_avg['inner_nodes'], data_xcmp_avg['time_server'], color+'--', label=f"XCMP n={bitlength}")
        plt_comm.plot(data_xcmp_avg['inner_nodes'], data_xcmp_avg['comm_request'], color+'--', label=f"XCMP n={bitlength}")

        print("{:0.0f}".format(        
            data_xcmp_avg['comm_request'].mean()
            # data_xcmp_avg['comm_request'].mean()-data_xcmp_avg['comm_request'].std(),
            # data_xcmp_avg['comm_request'].mean()+data_xcmp_avg['comm_request'].std()
        ), end=" & ")

        if bitlength <= 11:
            selected = data_sortinghats_avg[(data_sortinghats_avg['num_attributes']==num_attr) & (data_sortinghats_avg['internal_count']<=2000)]
            
            plt_time.plot(selected['internal_count'], selected['duration'], 'r-.', label=f"Sortinghats")
            plt_comm.plot(selected['internal_count'], selected['comm'], 'r-.', label=f"Sortinghats")

            print("{:0.0f}".format(        
                selected['duration'].mean()
            ), end="")
        else:
            print("-", end="")
        
        print(r"\\")
            

    figures_dir=os.path.join(PROJECT_ROOT, results_path, 'figures/sythetic')
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    plt_time.set_ylabel("Milliseconds")
    plt_time.set_xlabel("Number of Decision Nodes")
    plt_time.set_xscale('log')
    plt_time.set_xticks([4,8,16,32,64,128,256,512,1024], [4,8,16,32,64,128,256,512,1024])
    plt_time.set_ylim(-50,10000)
    # plt_time.legend()
    plt_time.set_title("Inference Time versus Number of Decision Nodes")
    f1.savefig(os.path.join(figures_dir, f"time-num-attr-{num_attr}.pdf"), bbox_inches = 'tight')
    plt.close()

    plt_comm.set_ylabel("KBytes")
    plt_comm.set_xlabel("Number of Decision Nodes")
    # plt_comm.legend()
    plt_comm.set_title("Communication versus Number of Decision Nodes")
    plt_comm.set_yscale('log')
    f2.savefig(os.path.join(figures_dir, f"comm-num-attr-{num_attr}.pdf"), bbox_inches = 'tight')
    plt.close()
