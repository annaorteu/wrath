#!/usr/bin/env python
# Description: This script takes a matrix file and a list of outliers and outputs a list of SVs and a heatmap
# Usage: python sv_detection_and_heatmap.py -m matrix_file -o outliers_file -s output_file -f window_size -c chromosome -w window_file -p plot_file
# Input: matrix_file = matrix file with genomic windows as row and column names
#        outliers_file = list of outliers with row and column numbers
#        window_size = size of genomic windows
#        chromosome = chromosome name
#        window_file = file with genomic window positions
# Output: output_file = list of SVs with start and end positions and length in genomic windows
#         plot_file = heatmap plot
# Modules required: argparse, pandas, numpy, matplotlib, seaborn, sklearn
# Date: 27 September 2023
# Author: Anna Orteu
#########################################################################################################################

import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering

#########################################################################################################################

### parse arguments

parser = argparse.ArgumentParser()

#input and output files
parser.add_argument("-m", "--matrix", help="Input matrix", action = "store")
parser.add_argument("-o", "--outliers", help="Input detected outliers", action = "store")
parser.add_argument("-p", "--plot", help="Output heatmap plot", action = "store")
parser.add_argument("-s", "--outFile", help="Output SVs", action = "store")
parser.add_argument("-f", "--winSize", help="Window size", type=int, action = "store", default = 1)
parser.add_argument("-c", "--chromosome", help="Chromosome name", type=str, action = "store")
parser.add_argument("-w", "--windowFile", help="Input genomic windows file", action = "store")

args = parser.parse_args()


#########################################################################################################################

#open files

matrix_file = pd.read_csv(args.matrix, sep=',', lineterminator='\n', header=None)
outliers_file = pd.read_csv(args.outliers,sep=',', lineterminator='\n')
window_file = pd.read_csv(args.windowFile, sep='\t', lineterminator='\n', header=None)
outplot = args.plot
output = args.outFile
window_size = args.winSize
chrom = args.chromosome


#########################################################################################################################

#rename axis based on genomic window positions
matrix_file.rename(index=window_file[1], columns=window_file[1], inplace=True)

#########################################################################################################################

#plot settings
plt.rcParams['figure.figsize'] = [30, 30] #set figure size

#clustering of outliers by proximity. Max distance between pairs of points set to 3
points=outliers_file[["ncol","nrow"]].to_numpy()
clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=3, linkage='single').fit(points)
outliers_file.loc[:, 'group']=clustering.labels_

# Identification of breakpoints
groups=outliers_file.groupby(['group'])
d={'mincol':groups.min(['ncol'])['ncol'], 'maxcol':groups.max(['ncol'])['ncol'], 'minrow':groups.min(['nrow'])['nrow'], 'maxrow':groups.max(['nrow'])['nrow']}
breakPoints = pd.DataFrame(data=d)

#plot heatmap in half a triangle and the detected outliers in the other
sns.set(font_scale=3)

data = np.log(matrix_file + 0.0001) * 100 #transform the data to log scale and multiply by 100 to get a percentage

heatmap_plot = sns.heatmap(data, cmap="YlGnBu", square=True, cbar_kws={'label': 'Barcode sharing %', 'shrink': 0.5})
heatmap_plot.scatter(x=breakPoints['minrow'], y=breakPoints['mincol'], color='k')
heatmap_plot.scatter(x=breakPoints['maxrow']+1, y=breakPoints['maxcol']+1, color='m')

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
fig.savefig(outplot)

#calculate lengths of svs and sort them by length
breakPoints['length']=breakPoints['maxcol']-breakPoints['minrow']
breakPoints.sort_values(by=['length'], ascending=False, inplace=True)

output_table={'SV_id':breakPoints.index.values, 'chromosome':chrom, 'start':breakPoints['minrow']*window_size, 'end':breakPoints['maxcol']*window_size, 'length':breakPoints['length']*window_size}
output_df=pd.DataFrame(output_table)
output_df.to_csv(output,index=False)
