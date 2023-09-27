#!/usr/bin/env python

import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#########################################################################################################################

### parse arguments

parser = argparse.ArgumentParser()

#input and output files
parser.add_argument("-m", "--matrix", help="Input matrix", action = "store")
parser.add_argument("-o", "--outFile", help="Output heatmap file", action = "store")
parser.add_argument("-w", "--windowFile", help="Input genomic windows file", action = "store")

args = parser.parse_args()


#########################################################################################################################

#open files

matrix_file = pd.read_csv(args.matrix, sep=',', lineterminator='\n', header=None)
window_file = pd.read_csv(args.windowFile, sep='\t', lineterminator='\n', header=None)
output = args.outFile


#########################################################################################################################

#rename axis based on genomic window positions
matrix_file.rename(index=window_file[1], columns=window_file[1], inplace=True)


#########################################################################################################################

#plot and save output

#plot and save output
plt.rcParams['figure.figsize'] = [30, 30]
sns.set(font_scale=3)

data = np.log(matrix_file + 0.0001) * 100 #transform the data to log scale and multiply by 100 to get a percentage

heatmap_plot = sns.heatmap(data, cmap="YlGnBu", square=True, cbar_kws={'label': 'Barcode sharing %', 'shrink': 0.5})


# Calculate appropriate tick locations and labels for the y-axis (rows)
num_ticks_y = round(len(data.index) / 60)  # Adjust the number of desired ticks
tick_locs_y = np.around(np.linspace(0, len(data.index) - 1, num_ticks_y)).astype(int)
tick_labels_y = data.index[tick_locs_y] # get the labels from the index (chromosome positions)

heatmap_plot.set_yticks(tick_locs_y)
heatmap_plot.set_yticklabels(tick_labels_y)

# Calculate appropriate tick locations and labels for the x-axis (columns)
num_ticks_x = round(len(data.columns) / 60)  # Adjust the number of desired ticks
tick_locs_x = np.around(np.linspace(0, len(data.columns) - 1, num_ticks_x)).astype(int)
tick_labels_x = data.columns[tick_locs_x] # get the labels from the index (chromosome positions)

heatmap_plot.set_xticks(tick_locs_x) 
heatmap_plot.set_xticklabels(tick_labels_x)

# Calculate appropriate tick locations and labels for the colorbar
cbar = heatmap_plot.collections[0].colorbar
num_ticks = 2  # Adjust the number of desired ticks

# Ensure that tick_locs is a 1D array
tick_locs = np.linspace(np.nanmin(data.to_numpy()), np.nanmax(data.to_numpy()), num_ticks)

# Transform tick_locs back to the original scale

cbar.set_ticks(tick_locs)
cbar.set_ticklabels(['0','100'])

# Plot the figure
fig = heatmap_plot.get_figure()
fig.figure.axes[-1].yaxis.label.set_size(50)
fig.savefig(output)
