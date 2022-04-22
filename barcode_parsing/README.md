# Parsing barcodes in haplotagging sequence data

Illumina data from a haplotagging run consists of four reads:
* forward and reverse **sequencing** reads
* forward and reverse **index** reads

(If you are starting from raw illumina files that have not yet been separated into the four read types, see XXX)

The unique barcodes that tag reads from the same molecule are stored in the index reads. `parse_haptag_barcodes.py` will identify the barcodes in the index reads and add this information into the fastq file headers of the forward and reverse sequencing reads. This information can then be used by toos such as *WRATH*.

Each read will be labelled with four codes, one for each of the A, B, C and D barcodes. These are stored in text files `BC_A.txt` etc (see in this repository).

## Assumptions

This code currently assumes the following:

* Both index files are 13 bp long
* Index read 1 carries the C and A barcodes in the forward orientation, arranged as `CCCCCCNAAAAAA`
* Index read 2 carries the D and B barcodes in the forward orientation, arranged as `DDDDDDNBBBBBB`

Currently, there are no built-in switches to change this, but these can be easily implemented by editing the following lines in 

```
    I1 = get_read(index_files[0], reverse_complement=False)
    I2 = get_read(index_files[1], reverse_complement=False)
```
and
```
    barcodes = {"C":I1.seq[:6], "A":I1.seq[7:], "D": I2.seq[:6], "B":I2.seq[7:]}
```

## Parse barcode information without demultiplexing separate individuals

To generate a single pair of fastq files with the barcodes tagged, run the script as follows:

```
python parse_haptag_barcodes.py -R read1.fq.gz read2.fq.gz -I index1.fq.gz index2.fq.gz --output_label my_experiment --output_dir /path/to/ouput/
```

This will produce four output files in the specified output directory:
```
my_experiment.R1.fastq.gz
my_experiment.R2.fastq.gz
unassigned.my_experiment.R1.fastq.gz
unassigned.my_experiment.R2.fastq.gz
```

The two "unassigned" files represent reads in which one or more of the barcode sequences could not be unambiguosly identified. By default this script will allow one base mismatches in barcode sequences, provided they are unambiguous. Add option `--exact_macth_only` to use only exact matches (this will reduce the total number of assigned reads).

Add option `--count_barcodes` to produce an addition output giving the number of times each barcode was seen.

If your barcode files are not named `BC_A.txt` etc. you will have to add the option `--barcode_files` followd by the file names of the four barcodes files.

## Parse barcodes and demultiplexing separate individuals

To generate a separate pair of fastq files for each individual, we need a **demultiplexing file**. This gives the C code (`C01-C96`) for each individual. If you multiplexed fewer han 96 individuals on a lane, you might have multiple C codes for each individual.

To generate the demultiplexing file, you need a file giving the barcode(s) sequence(s) for each individual. See the example `sample_barcodes_example.txt` in this repository. Note how some individuals in this example have multiple barcodes whereas others have just one. This is just to demonstrate that such an arrangement is allowed.

Generate the demultiplexing file using the script `make_demult_file.py`:

```
python make_demult_file.py --C_barcode_file BC_C.txt --sample_barcode_file sample_barcodes_example.txt -o demult_file_example.txt
```

Now we are ready to parse the barcodes in the fastq reads:

```
python parse_haptag_barcodes.py -R read1.fq.gz read2.fq.gz -I index1.fq.gz index2.fq.gz --output_label my_experiment_label --demult_file demult_file.txt --output_dir /path/to/ouput/
```

This will produce four output files per individual. For example:

```
sample1.my_experiment.R1.fastq.gz
sample1.my_experiment.R2.fastq.gz
sample1.unassigned.my_experiment.R1.fastq.gz
sample1.unassigned.my_experiment.R2.fastq.gz
```

Here, the "unassigned" reads are those that could be assigned to an individual (the C barcode could be matched unambiguously), but one of the other barcodes was ambiguous, so the read could not be assigned to a molecule.

There will also be two files for unassigned reads that could not be unambiguously matched to a C barcode and therefore an individual:

```
unassigned.my_experiment.R1.fastq.gz
unassigned.my_experiment.R2.fastq.gz
```

## Processing raw Illumina files

The scripts above use separate sequencing and index reads. If you need to generate these yourself from the raw Illumina data, you need to first ensure that you understand how many bases to expect in each read and index read.

The files are processed using Illumina's [`bcl2fastq`](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html).

For 150 bp reads with 13 bp index reads, the command would look something like:

```
bcl2fastq --barcode-mismatches 1 --use-bases-mask 'Y150n,I13,I13,Y150n' --mask-short-adapter-reads 1 --create-fastq-for-index-reads -R /path/to/data -o /path/to/output/
```



