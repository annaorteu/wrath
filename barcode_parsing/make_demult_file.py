#!/usr/bin/env python

# ---------------------------------------------------------
#                  Simon Martin, 2022
#                simon.martin@ed.ac.uk
#     Script to make demultiplexing input file for
#               parse_haptag_barcodes.py
# ---------------------------------------------------------

# example usage:
# python make_demult_file.py --C_barcode_file BC_C.txt --sample_barcode_file sample_barcodes_example.txt -o demult_file_example.txt

import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--C_barcode_file", help="File giving the barcode for each C code", default="BC_C.txt")
parser.add_argument("-s", "--sample_barcode_file", help="File giving the barcode or each sample", required=True)
parser.add_argument("-o", "--output_demult_file", help="Output demultiplex file")

args = parser.parse_args()

with open(args.C_barcode_file, "rt") as cbf:
    barcode_dict = dict([line.split() for line in cbf])

with open(args.sample_barcode_file, "rt") as bsf:
    sample_dict = dict([line.split() for line in bsf])

with open(args.output_demult_file, "wt") if args.output_demult_file else sys.stdout as df:
    for code in barcode_dict.keys():
        try: df.write(code + "\t" + sample_dict[barcode_dict[code]] + "\n")
        except: continue
