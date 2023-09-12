import os

import matplotlib.pyplot as plt
import pandas as pd
import yaml
from matplotlib.patches import ConnectionPatch, Rectangle
from commons import *
# PROJECT_ROOT='/home/r5akhava/private-decision-tree-evaluation/'

data = pd.read_csv(
        os.path.join(PROJECT_ROOT,f"experiments/results-{version_tag}/cmp/cmp_bench.csv"),
        # correct,bitlength,comparison,log_deg,hamming_weight,code_length,num_cmps,total_time,amortized_time
        names=["correct","bitlength","comparison","log_deg","hamming_weight","code_length","num_cmps","total_time","amortized_time"],
    )
data['amortized_time'] =data['amortized_time']/1000


fig, ax = plt.subplots()
plt.xlabel("Bit Precision")
plt.ylabel("Milliseconds")
plt.title("Amortized Runtime for Comparison")

# # Adjust figure size and main plot position
fig.subplots_adjust(left=0.1, bottom=0.1, right=1, top=0.9)

# Create separate Axes object for the magnified part and place it outside the main plot
axins = fig.add_axes([0.68, 0.5, 0.3, 0.3])  # left, bottom, width, height

# Folklore
data1=data[data['comparison']==1]
data1_avg = data1.groupby(['bitlength'], as_index=False).mean()
data1_std = data1.groupby(['bitlength'], as_index=False).std()
ax.plot(data1_avg['bitlength'], data1_avg['amortized_time'], label="Folklore")
ax.fill_between(
    data1_avg['bitlength'],
    data1_avg['amortized_time']-data1_std['amortized_time'],
    data1_avg['amortized_time']+data1_std['amortized_time'],
    alpha=0.1
)
axins.plot(data1_avg['bitlength'], data1_avg['amortized_time'])
axins.fill_between(
    data1_avg['bitlength'],
    data1_avg['amortized_time']-data1_std['amortized_time'],
    data1_avg['amortized_time']+data1_std['amortized_time'],
    alpha=0.1
)

# Range Cover
for hw, x_low, x_high in [
        (2,None, None),
        (4,None, 22),
        (8,None, 34),
        (16,None, None),
        (32,None, None),
    ]:
    data1=data[(data['comparison']==0) & (data['hamming_weight']==hw)]
    if x_low is not None:
        data1 = data1[data1['bitlength']>=x_low]
    if x_high is not None:
        data1 = data1[data1['bitlength']<=x_high]
    
    data1_avg = data1.groupby(['bitlength'], as_index=False).mean()
    data1_std = data1.groupby(['bitlength'], as_index=False).std()
    ax.plot(data1_avg['bitlength'], data1_avg['amortized_time'], label=f'RCC (h={hw})')
    ax.fill_between(
        data1_avg['bitlength'],
        data1_avg['amortized_time']-data1_std['amortized_time'],
        data1_avg['amortized_time']+data1_std['amortized_time'],
        alpha=0.1
    )
    axins.plot(data1_avg['bitlength'], data1_avg['amortized_time'], label=f'RCC (h={hw})')
    axins.fill_between(
        data1_avg['bitlength'],
        data1_avg['amortized_time']-data1_std['amortized_time'],
        data1_avg['amortized_time']+data1_std['amortized_time'],
        alpha=0.1
    )

# XCMP
data1=data[data['comparison']==2]
data1_avg = data1.groupby(['bitlength'], as_index=False).mean()
data1_std = data1.groupby(['bitlength'], as_index=False).std()
ax.plot(data1_avg['bitlength'], data1_avg['amortized_time'], label="XXCMP")
ax.fill_between(
    data1_avg['bitlength'],
    data1_avg['amortized_time']-data1_std['amortized_time'],
    data1_avg['amortized_time']+data1_std['amortized_time'],
    alpha=0.1
)
axins.plot(data1_avg['bitlength'], data1_avg['amortized_time'])
axins.fill_between(
    data1_avg['bitlength'],
    data1_avg['amortized_time']-data1_std['amortized_time'],
    data1_avg['amortized_time']+data1_std['amortized_time'],
    alpha=0.1
)

# SortingHats (base on their paper)
ax.plot(11, 0.4, "o", label="SortingHats")
axins.plot(11, 0.4, "o")


# Ilia et al.
avg_times = []

# Open and read the file
filename='results-v9/ilia/case1.txt'
with open(filename, 'r') as file:
    lines = file.readlines()

    # Iterate over the lines
    for line in lines:
        # If the line contains 'Avg. time per integer'
        if "Avg. time per integer" in line:
            # Extract the value and append to the list
            value = float(line.split(':')[1].split('ms')[0].strip())
            avg_times.append(value)

# Calculate and return the average
ilia = sum(avg_times) / len(avg_times)
ax.plot(64, ilia, "o", label="Iliashenko et al.")
axins.plot(64, ilia, "o")

#put the legend top left
ax.legend(loc='upper left')

# axins.set_yticks([0,2,4,6,8])
# axins.set_yticklabels([0,2,4,6,8])


# Set the limits for the magnifying glass
x1, x2, y1, y2 = 4, 37, -0.5, 5
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# Add a title to the inset plot
# axins.set_title("Zoomed In")

# Add a box (rectangle) around the magnified part of the main plot
rect = Rectangle((x1, y1), x2 - x1, y2 - y1, linewidth=1, edgecolor='grey', facecolor='none', linestyle="--")
ax.add_patch(rect)

# Create lines connecting the magnified part of the main plot to the magnified part's box
xyA = (x2, y1)
xyB = (x2, y2)
coordsA = "data"
coordsB = "data"
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA=coordsA, coordsB=coordsB,
                      axesA=axins, axesB=ax, arrowstyle="-", linestyle="--", color="grey", linewidth=1)
axins.add_artist(con)

xyA = (x1, y1)
xyB = (x1, y2)
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA=coordsA, coordsB=coordsB,
                      axesA=axins, axesB=ax, arrowstyle="-", linestyle="--", color="grey", linewidth=1)
axins.add_artist(con)

# Display the plot
# plt.show()
plt.savefig(f"results-{version_tag}/figures/cmp_compare_all.pdf", bbox_inches='tight')