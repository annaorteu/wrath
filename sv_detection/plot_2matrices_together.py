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
plt.rcParams['figure.figsize'] = [30, 30] #set figure size

heatmap_plot = sns.heatmap(np.log(matrix_file3+0.0001)*100, cmap="YlGnBu", square=True)
fig = heatmap_plot.get_figure()
fig.savefig(output)
