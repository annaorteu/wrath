#!/usr/bin/env python
# Description: This script takes a barcode file (bed) and a list of windows (bed) and outputs a jaccard matrix of barcode sharing between windows
# Usage: python jaccard_matrix.py -w window_file -b barcode_file -o output_file -t threads
# Input: window_file = file with genomic window positions
#        barcode_file = file with barcodes and positions
# Output: output_file = jaccard matrix
# Modules required: argparse, sys, gzip, random, pysam, math, numpy, pandas
# Date: 27 September 2023
# Author: Anna Orteu
#########################################################################################################################

import argparse, sys, gzip, random, pysam, math
import numpy as np
import pandas as pd

from threading import Thread

from multiprocessing import Process

if sys.version_info>=(3,0):
    from multiprocessing import SimpleQueue
else:
    from multiprocessing.queues import SimpleQueue

from time import sleep

import time
start_time = time.time()


#########################################################################################################################


### parse arguments

parser = argparse.ArgumentParser()

#input and output files
parser.add_argument("-w", "--winFile", help="Input window file", action = "store")
parser.add_argument("-b", "--barcodeFile", help="Input barcode file", action = "store")
parser.add_argument("-o", "--outFile", help="Output jaccard matrix file", action = "store")

#other
parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("--test", help="Test - runs 10 windows", action='store_true')
parser.add_argument("--verbose", help="Verbose output", action = "store_true")

args = parser.parse_args()


#########################################################################################################################

#open files

if args.outFile:
    outFile = open(args.outFile, "wt")
else: outFile = sys.stdout

#read windows
windowFile = pd.read_csv(args.winFile, sep='\t', lineterminator='\n', header=None)
num_win = windowFile.shape[0]

#create a matrix of n x n, n = number of windows to compare
windowFile=pd.DataFrame(windowFile)
windowFile.index.name = 'index'
windowFile.reset_index(inplace=True)

#read barcodes
tbx = pysam.TabixFile(args.barcodeFile)


#########################################################################################################################

#functions

'''A function that reads from the input queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def freqs_wrapper(inQueue, resultQueue, number_win, inFile):
    while True:
        windowNumber,windowLine = inQueue.get() # retrieve window
        if windowNumber == -1:
            resultQueue.put((-1,None,)) # this is the way of telling everything we're done
            break
        array_i = np.zeros((number_win))
        array_u = np.ones((number_win))
        bedfile1 = inFile.fetch(windowLine[0], windowLine[1], windowLine[2],  parser=pysam.asBed(), multiple_iterators=True)
        barcodes1 = [rowbed1.name for rowbed1 in bedfile1]
        for index2, row2 in windowFile.iloc[windowNumber:,:].iterrows():
            bedfile2 = inFile.fetch(row2[0], row2[1], row2[2], parser=pysam.asBed(), multiple_iterators=True)
            barcodes2 = [rowbed2.name for rowbed2 in bedfile2]
            intersect = np.intersect1d(barcodes1, barcodes2)
            union = np.union1d(barcodes1, barcodes2)
            array_i[index2] = intersect.size
            if union.size>1:
                array_u[index2] = union.size
            elif union.size<1:
                array_u[index2] = 1
                array_i[index2] = 0
            array_u[index2] = union.size
        outArray = np.divide(array_i, array_u)
        resultQueue.put((windowNumber, outArray,))


'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the line number, and this increases consecutively.'''
def sorter(doneQueue, writeQueue, verbose, nWorkerThreads):
    global resultsReceived
    sortBuffer = {}
    expect = 0
    threadsComplete = 0 #this will keep track of the worker threads and once they're all done this thread will break
    while True:
        windowNumber, results = doneQueue.get()
        #check if we're done
        if windowNumber == -1: threadsComplete += 1
        if threadsComplete == nWorkerThreads:
            writeQueue.put((-1,None,))
            break #this is the way of telling everything we're done
        resultsReceived += 1
        if verbose:
            sys.stderr.write("Sorter received window {}\n".format(windowNumber))
        if windowNumber == expect:
            writeQueue.put((windowNumber,results))
            if verbose:
                sys.stderr.write("window {} sent to writer\n".format(windowNumber))
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    results = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,results))
                    if verbose:
                        sys.stderr.write("window {} sent to writer\n".format(expect))
                    expect +=1
                except:
                    break
        else:
            #otherwise this line is ahead of us, so add to buffer dictionary
            sortBuffer[str(windowNumber)] = results



'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, out, verbose):
    global resultsWritten
    while True:
        windowNumber, results = writeQueue.get()
        #check if we're done
        if windowNumber == -1: break
        if verbose:
            sys.stderr.write("Writer received window {}\n".format(windowNumber))
        np.savetxt(out, results, fmt='%.10f', newline=',')
        out.write("\n")
        resultsWritten += 1

'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        sys.stderr.write("{} windows queued | {} windows analysed | {} windows written\n".format(windowQueued,resultsReceived,resultsWritten))





#########################################################################################################################

#counting stat that will let keep track of how far we are
windowQueued = 0
resultsReceived = 0
resultsWritten = 0
linesWritten = 0

'''Create queues to hold the data one will hold the line info to be passed to the analysis'''
inQueue = SimpleQueue()
#one will hold the results (in the order they come)
resultQueue = SimpleQueue()
#one will hold the sorted results to be written
writeQueue = SimpleQueue()


'''start worker Processes for analysis. The command should be tailored for the analysis wrapper function
of course these will only start doing anything after we put data into the line queue
the function we call is actually a wrapper for another function.(s) This one reads from the line queue, passes to some analysis function(s), gets the results and sends to the result queue'''
workerThreads = []
sys.stderr.write("\nStarting {} worker threads\n".format(args.threads))
for x in range(args.threads):
  workerThread = Process(target=freqs_wrapper, args = (inQueue, resultQueue, num_win, tbx,))
  workerThread.daemon = True
  workerThread.start()
  workerThreads.append(workerThread)


'''thread for sorting results'''
sorterThread = Thread(target=sorter, args=(resultQueue, writeQueue, args.verbose, args.threads,))
sorterThread.daemon = True
sorterThread.start()

'''start thread for writing the results'''
writerThread = Thread(target=writer, args=(writeQueue, outFile, args.verbose,))
writerThread.daemon = True
writerThread.start()

'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
checkerThread = Thread(target=checkStats)
checkerThread.daemon = True
checkerThread.start()

#########################################################################################################################


if not args.test:
    for windowIdx, windowLine in windowFile.iterrows():
        inQueue.put((windowQueued,windowLine))
        windowQueued += 1
else:
    for windowIdx, windowLine in windowFile.iterrows():
        inQueue.put((windowQueued,windowLine))
        windowQueued += 1
        if windowQueued == 10: break


#########################################################################################################################

#Now we send completion signals to all worker threads
for x in range(args.threads):
    inQueue.put((-1,None,)) # -1 tells the threads to break

sys.stderr.write("\nClosing worker threads\n".format(args.threads))
for x in range(len(workerThreads)):
    workerThreads[x].join()

sorterThread.join()
writerThread.join()

sys.stderr.write("\nDone\n")

sys.stderr.write("My program took {} to run\n".format(time.time() - start_time))

outFile.close()

sys.exit()
