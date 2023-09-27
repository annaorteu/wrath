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
parser.add_argument("-m1", "--matrix1", help="Input first matrix", action = "store")
parser.add_argument("-m2", "--matrix2", help="Input second matrix", action = "store")
parser.add_argument("-o", "--outFile", help="Output heatmap file", action = "store")
parser.add_argument("-w", "--windowFile", help="Input genomic windows file", action = "store")

args = parser.parse_args()


#########################################################################################################################

#open files

matrix_file1 = pd.read_csv(args.matrix1, sep=',', lineterminator='\n', header=None)
matrix_file2 = pd.read_csv(args.matrix2, sep=',', lineterminator='\n', header=None)
window_file = pd.read_csv(args.windowFile, sep='\t', lineterminator='\n', header=None)
output = args.outFile


#########################################################################################################################

#Transpose one of the matrices and join the two triangles
matrix_file3 = matrix_file1.transpose().add(matrix_file2)

#Drop the last column and last row as its full of NaNs
matrix_file3 = matrix_file3.iloc[:-1, :-1]

#rename axis based on genomic window positions
matrix_file3.rename(index=window_file[1], columns=window_file[1], inplace=True)


#########################################################################################################################

#plot and save output
plt.rcParams['figure.figsize'] = [30, 30]
sns.set(font_scale=3)

data = np.log(matrix_file3 + 0.0001) * 100

heatmap_plot = sns.heatmap(data, cmap="YlGnBu", square=True, cbar_kws={'label': 'Barcode sharing %', 'shrink': 0.5})


# Calculate appropriate tick locations and labels for the y-axis (rows)
num_ticks_y = round(len(data.index) / 60)  # Adjust the number of desired ticks
tick_locs_y = np.around(np.linspace(0, len(data.index) - 1, num_ticks_y)).astype(int)
tick_labels_y = data.index[tick_locs_y]

heatmap_plot.set_yticks(tick_locs_y)
heatmap_plot.set_yticklabels(tick_labels_y)

# Calculate appropriate tick locations and labels for the x-axis (columns)
num_ticks_x = round(len(data.columns) / 60)  # Adjust the number of desired ticks
tick_locs_x = np.around(np.linspace(0, len(data.columns) - 1, num_ticks_x)).astype(int)
tick_labels_x = data.columns[tick_locs_x]

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
fig.figure.axes[-1].tick_params(labelsize=40)
fig.savefig(output)
