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

args = parser.parse_args()


#########################################################################################################################

#open files

matrix_file = pd.read_csv(args.matrix, sep=',', lineterminator='\n', header=None)
output = args.outFile


#########################################################################################################################

#plot and save output 
plt.rcParams['figure.figsize'] = [30, 30] #set figure size

heatmap_plot = sns.heatmap(np.log(matrix_file+0.0001)*100, cmap="YlGnBu", square=True)
fig = heatmap_plot.get_figure()
fig.savefig(output) 