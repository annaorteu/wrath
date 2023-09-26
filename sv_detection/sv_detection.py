#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering

#########################################################################################################################

### parse arguments

parser = argparse.ArgumentParser()

#input and output files
parser.add_argument("-m", "--matrix", help="Input matrix", action = "store")
parser.add_argument("-o", "--outliers", help="Input detected outliers", action = "store")
parser.add_argument("-s", "--outFile", help="Output SVs", action = "store")
parser.add_argument("-w", "--winSize", help="Window size", type=int, action = "store", default = 1)

args = parser.parse_args()


#########################################################################################################################

#open files

matrix_file = pd.read_csv(args.matrix, sep=',', lineterminator='\n', header=None)
outliers_file = pd.read_csv(args.outliers,sep=',', lineterminator='\n')
output = args.outFile
window_size = args.winSize


#########################################################################################################################

#clustering of outliers by proximity. Max distance between pairs of points set to 3
points=outliers_file[["ncol","nrow"]].to_numpy()
clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=3, linkage='single').fit(points)
outliers_file.loc[:, 'group']=clustering.labels_

# Identification of breakpoints
groups=outliers_file.groupby(['group'])
d={'mincol':groups.min(['ncol'])['ncol'], 'maxcol':groups.max(['ncol'])['ncol'], 'minrow':groups.min(['nrow'])['nrow'], 'maxrow':groups.max(['nrow'])['nrow']}
breakPoints = pd.DataFrame(data=d)

#calculate lengths of svs and sort them by length
breakPoints['length']=breakPoints['maxcol']-breakPoints['minrow']
breakPoints.sort_values(by=['length'], ascending=False, inplace=True)

output_table={'SV_id':breakPoints.index.values, 'start':breakPoints['minrow']*window_size, 'end':breakPoints['maxcol']*window_size, 'length':breakPoints['length']*window_size}
output_df=pd.DataFrame(output_table)
output_df.to_csv(output)
