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

matrix_file.index.name = 'nrow'
matrix_file.reset_index(inplace=True)

matrix_file=matrix_file.drop(columns=[50]) #drop last column, full of nans
df_file = pd.melt(matrix_file, id_vars="nrow", var_name='ncol') #melt dataframe 
df_file['newcol'] = df_file['ncol']+df_file['nrow'] #shift columns by row index to get positions around diagonal
df_file = df_file.drop(columns=['ncol']) #drop old column names
df_file = df_file.pivot(index='nrow', columns='newcol') #spread the dataframe
df_file = df_file.fillna(0) #substitute nans for 0s

heatmap_plot = sns.heatmap(np.log(df_file+0.0001)*100, cmap="YlGnBu", square=True)
fig = heatmap_plot.get_figure()
fig.savefig(output) 


