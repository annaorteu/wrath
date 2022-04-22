#!/usr/bin/env python

# ---------------------------------------------------------
#                  Simon Martin, 2022
#                simon.martin@ed.ac.uk
#     Script to parse haplotag parcodes in index reads
# ---------------------------------------------------------

#example_command
# python parse_haptag_barcodes.py -R read1.fq.gz read2.fq.gz -I index1.fq.gz index2.fq.gz \
# --output_label my_experiment_label --demult_file demult_file.txt --output_dir /path/to/ouput/ --count_barcodes

import argparse
import gzip
from collections import defaultdict
from multiprocessing import Process, SimpleQueue

#a simple read object
class Read:
    def __init__(self,name,seq,qual):
        self.name = name
        self.seq = seq
        self.qual = qual
    
    def as_text(self):
        return "\n".join([self.name, self.seq, "+", self.qual])

#this class is a bit like a defaultdict, but it doesn not store the missing thing 
class missing_dict(dict):
    def __init__(self, missing):
        self.missing=missing
    def __missing__(self, key):
        return self.missing

#another dictionary class that just returns the key if it is missing
class mirror_dict(dict):
    def __missing__(self, key):
        return key

#function to read in a sequence read from a fastq file (4 lines)
def get_read(readFile, reverse_complement=False):
    name = readFile.readline().strip()
    seq = readFile.readline().strip()
    plus = readFile.readline()
    qual = readFile.readline().strip()
    if reverse_complement:
        seq = revComp(seq)
        qual = qual[::-1]
    return Read(name,seq,qual)

#reverse complement function (currently not using this, but it's here if needed.
complementTrans = str.maketrans("ACGT", "TGCA")

def revComp(seq):
    return seq.translate(complementTrans)[::-1]

#argsort function (borrowed) for writing observed barcodes in order
def argsort(seq):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__)

def fastq_writer(queue, outfile_name):
    outfile = gzip.open(outfile_name, "wt")
    while True:
        read = queue.get()
        if read == None: break #for ending process
        outfile.write(read.as_text() + "\n")
    outfile.close()

#function to make a dictionary for matching a barcode (exactly or with one mismatch)
def make_barcode_dict(barcode_file, missing):
    d = missing_dict(missing)
    with open(barcode_file, "rt") as BCfile:
        for line in BCfile:
            tag,bc = line.split()
            if bc in d:
                print("WARNING: Barcode {} is linked with both {} and {}.".format(bc, d[bc], tag))
                d.pop(bc)
                continue
            d[bc] = tag
    return(d)

#make a dictionary to return the correct barcode in the case of a mismatch
def make_barcode_correction_dict(barcodes, missing):
    #a default dictionary that returns the requested key if the item is not in the dictionary
    d = mirror_dict()
    redundant_alt_bcs = set()
    for bc in barcodes:
        #add actual barcide to dictionary
        d[bc] = bc
        for i in range(len(bc)):
            #for each position, try all possible alternative bases
            for base in "ACGT":
                if base != bc[i]:
                    alt_bc = bc[:i] + base + bc[i+1:]
                    if alt_bc in barcodes:
                        print("WARNING: ignoring alternative barcode {} for {} because it matches an existing barcode.".format(alt_bc, bc))
                        continue
                    if alt_bc in d:
                        print("WARNING: ignoring alternative barcode {} because it matches both {} and {}.".format(alt_bc, bc, d[alt_bc]))
                        redundant_alt_bcs.add(alt_bc)
                        continue
                    d[alt_bc] = bc
    #remove problem alternative barcodes that were matched by subsequent ones
    for alt_bc in redundant_alt_bcs:
        d.pop(alt_bc)
    return(d)


###############################################################################

#arguments
parser = argparse.ArgumentParser()
parser.add_argument("-R", "--read_files", help="Read files", nargs=2, required = True)
parser.add_argument("-I", "--index_read_files", help="Index read files", nargs=2, required = True)

parser.add_argument("--output_label", help="Output file label", default="")
parser.add_argument("--output_dir", help="Output file directory", default=".")

parser.add_argument("--barcode_files", help="Barcode files for A B C D (separated by spaces)",
                    nargs=4, required = False, default=["BC_A.txt", "BC_B.txt", "BC_C.txt", "BC_D.txt"])

parser.add_argument("--demult_file", help="File giving C tag and sample name, if you want to demultiplex")

parser.add_argument("--exact_match_only", help="Only allow exact barcode matches", action="store_true")

parser.add_argument("--count_barcodes", help="Output counts for all barcodes seen", action="store_true")

args = parser.parse_args()

###############################################################################

#a dictionary of what to return for missing barcodes
missing = dict((x, x+"00",) for x in "ABCD")

#read the barcode files and make the dictionaries
barcode_files = dict(zip(["A","B","C","D"], args.barcode_files))

barcode_dicts = dict([(x, make_barcode_dict(barcode_files[x], missing = missing[x]),) for x in "ABCD"])

if not args.exact_match_only:
    barcode_correction_dicts = dict([(x, make_barcode_correction_dict(list(barcode_dicts[x].keys()),
                                                                      missing="N"),) for x in "ABCD"])

#if demultiplexing by individual, parse that file
if args.demult_file:
    with open(args.demult_file, "rt") as df:
        demult_dict = dict([line.split() for line in df])
    #check that all barcode tags are associated with an individual
    for barcode,tag in barcode_dicts["C"].items():
        if tag not in demult_dict:
            print("Warning: tag {} is not in the demultiplex file, so will be sent to unassigned output.".format(tag))
    
    samples=list(set(demult_dict.values()))
    print("Demultiplexing {} samples.".format(len(samples)))


#open input files
index_files = [gzip.open(f, "rt") for f in args.index_read_files]
read_files = [gzip.open(f, "rt") for f in args.read_files]

#start writer threads for output files
writer_procs = []
out_queues = {}

#output for unassigned reads
out_queues["unassignedR1"] = SimpleQueue()
writer = Process(target=fastq_writer, args = (out_queues["unassignedR1"],
                                              args.output_dir + "/" + ".".join(["unassigned", args.output_label, "R1.fastq.gz"])))
writer.daemon = True
writer.start()

out_queues["unassignedR2"] = SimpleQueue()
writer = Process(target=fastq_writer, args = (out_queues["unassignedR2"],
                                              args.output_dir + "/" + ".".join(["unassigned", args.output_label, "R2.fastq.gz"])))
writer.daemon = True
writer.start()


#if demultiplexing, set up output files for each individual
if args.demult_file:
    for sample in samples:
        out_queues[sample] = {}
        out_queues[sample]["R1"] = SimpleQueue()
        writer = Process(target=fastq_writer, args = (out_queues[sample]["R1"],
                                                      args.output_dir + "/" + ".".join([sample, args.output_label, "R1.fastq.gz"])))
        writer.daemon = True
        writer.start()
        writer_procs.append(writer)
        
        out_queues[sample]["R2"] = SimpleQueue()
        writer = Process(target=fastq_writer, args = (out_queues[sample]["R2"],
                                                      args.output_dir + "/" + ".".join([sample, args.output_label, "R2.fastq.gz"])))
        writer.daemon = True
        writer.start()
        writer_procs.append(writer)
        
        #and reads that can be assigned to an individual, but not to a molecule
        out_queues[sample]["unassignedR1"] = SimpleQueue()
        writer = Process(target=fastq_writer, args = (out_queues[sample]["unassignedR1"],
                                                args.output_dir + "/" + ".".join([sample, "unassigned", args.output_label, "R1.fastq.gz"])))
        writer.daemon = True
        writer.start()
        writer_procs.append(writer)

        out_queues[sample]["unassignedR2"] = SimpleQueue()
        writer = Process(target=fastq_writer, args = (out_queues[sample]["unassignedR2"],
                                                args.output_dir + "/" + ".".join([sample, "unassigned", args.output_label, "R2.fastq.gz"])))
        writer.daemon = True
        writer.start()
        writer_procs.append(writer)
else:
    #if not demultiplexing start writers for all read 1s and read 2s
    out_queues["R1"] = SimpleQueue()
    writer = Process(target=fastq_writer, args = (out_queues["R1"],
                                                  args.output_dir + "/" + ".".join([args.output_label, "R1.fastq.gz"])))
    writer.daemon = True
    writer.start()
    writer_procs.append(writer)
    
    out_queues["R2"] = SimpleQueue()
    writer = Process(target=fastq_writer, args = (out_queues["R2"],
                                                  args.output_dir + "/" + ".".join([args.output_label, "R2.fastq.gz"])))
    writer.daemon = True
    writer.start()
    writer_procs.append(writer)


#if counting barcodes
if args.count_barcodes:
    #dictionary counter that defualts to zero for new barcodes
    barcodeCounts = dict([(x, defaultdict(int),) for x in "ABCD"])
    #file names for writing
    barcodeCounts_files = dict([(x, args.output_dir + "/" + args.output_label + "." + x + ".counts.txt",) for x in "ABCD"])

###############################################################################


#A loop that reads one read at a time (well, one each from R1, R2, I1 and I2
#and determines the BX tag based on the barcodes
#and then writes to either the ouput files or the unassigned files (if no match was found)
while True:
    I1 = get_read(index_files[0], reverse_complement=False)
    I2 = get_read(index_files[1], reverse_complement=False)
    R1 = get_read(read_files[0])
    R2 = get_read(read_files[1])
    
    #if we're at the end of the file, just stop
    if I1.name=="": break
    
    #extract barcodes A, B, C, D from indices
    barcodes = {"C":I1.seq[:6], "A":I1.seq[7:], "D": I2.seq[:6], "B":I2.seq[7:]}
    
    #if allowing mismatches, correct barcodes now    
    if not args.exact_match_only:
        for x in "ABCD":
            barcodes[x] = barcode_correction_dicts[x][barcodes[x]]
        #and actually reconstruct index sequence from corrected barcodes
        I1.seq = barcodes["C"] + I1.seq[6] + barcodes["A"]
        I2.seq = barcodes["D"] + I2.seq[6] + barcodes["B"]
    
    #if counting barcodes, add counts
    if args.count_barcodes:
        for x in "ABCD": barcodeCounts[x][barcodes[x]] += 1
    
    #Look up the codes based on the barcode sequences
    codes = dict([(x, barcode_dicts[x][barcodes[x]],) for x in "ABCD"])
    
    #check that we got a match for each code, otherwise mark this read as unassigned
    assigned = True
    for x in "ABCD":
        if codes[x] == missing[x]:
            assigned = False
            break
    
    #construct the BX, RX and QX tags
    BXtag = codes["A"] + codes["C"] + codes["B"] + codes["D"]
    
    RXtag = I1.seq + "+" + I2.seq
    
    QXtag = I1.qual + "+" + I2.qual
    
    name_extension = "_" + BXtag + "_" + I1.seq + "_" + revComp(I2.seq)
    
    #add extension and tags to read names
    R1.name += name_extension + "\tBX:Z:"+BXtag + "\t" + "RX:Z:"+RXtag + "\t" "QX:Z:"+QXtag
    R2.name += name_extension + "\tBX:Z:"+BXtag + "\t" + "RX:Z:"+RXtag + "\t" "QX:Z:"+QXtag
    
    #write to outputs
    if args.demult_file:
        #if demultiplexing, check individual code
        try:
            sample = demult_dict[codes["C"]]
            if assigned:
                out_queues[sample]["R1"].put(R1)
                out_queues[sample]["R2"].put(R2)
            else:
                #these are assignable to an sample (C tag), but not to a molecule (A, B and D)
                out_queues[sample]["unassignedR1"].put(R1)
                out_queues[sample]["unassignedR2"].put(R2)
        except:
            out_queues["unassignedR1"].put(R1)
            out_queues["unassignedR2"].put(R2)
    else:
        #otherwise write to assigned or unassigned output files
        if assigned:
            out_queues["R1"].put(R1)
            out_queues["R2"].put(R2)
        else:
            out_queues["unassignedR1"].put(R1)
            out_queues["unassignedR2"].put(R2)

#end writer processes
out_queues["unassignedR1"].put(None)
out_queues["unassignedR2"].put(None)

if args.demult_file:
    for sample in samples:
        out_queues[sample]["R1"].put(None)
        out_queues[sample]["R2"].put(None)
        out_queues[sample]["unassignedR1"].put(None)
        out_queues[sample]["unassignedR2"].put(None)
else:
    out_queues["R1"].put(None)
    out_queues["R2"].put(None)


#close all writers to close files and clear write buffers
for proc in writer_procs: proc.join()

#if writing barcode counts
if args.count_barcodes:
    for x in "ABCD":
        barcodeList = list(barcodeCounts[x].keys())
        barcodeCountList = [barcodeCounts[x][b] for b in barcodeList]
        codeList = [barcode_dicts[x][barcode] for barcode in barcodeList]
        order = argsort(barcodeCountList)[::-1]
        matches, mismatches = 0,0
        with open(barcodeCounts_files[x], "wt") as BCcounts_file:
            for i in order:
                output = [barcodeList[i], codeList[i]]
                if args.demult_file:
                    try: output.append(demult_dict[codeList[i]])
                    except: output.append("-")
                output.append(str(barcodeCountList[i]))
                BCcounts_file.write("\t".join(output) + "\n")
                if codeList[i] == missing[x]: mismatches += barcodeCountList[i]
                else: matches += barcodeCountList[i]
        
        print("\nBarcode " + x + ":")
        print("Matches =", matches)
        print("Misatches =", mismatches)
