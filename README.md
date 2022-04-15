# wrath
WRATH: WRapped Analysis of Tagged Haplotypes

<img src="wrath_logo.png" alt="logo" width="50%"/>


---
## Running wrath

The main script `wrath.sh` runs the main version of *WRATH*.

A typical command looks like:

```bash
wrath.sh -g genome_file.fa -c chromosome_name  -w 50000  -s list_of_bam_files.txt -t 15
```

Input options are:

```
wrath: wrapped analysis of tagged haplotypes

DESCRIPTION:
 Program produces a jaccard matrix camparing the barcode content between all pairs windows whithin a chromosome.

wrath [-h] [-g GENOMEFILE] [-c CHROMOSOMENAME] [-w WINDOWSIZE] [-s FILELIST] [-t THERADS] [-p] [-v] [-x STEP]

OPTIONS:
  -h                show this help text
  -g GENOMEFILE     reference genome
  -c CHROMOSOMENAME chromosome
  -w WINDOWSIZE     window size
  -s FILELIST       list of bam files with paths of the individuals of the population/phenotype of interest
  -t THERADS        threads to use
  -p                skip plotting the heatmap
  -x STEP           start from a given step. Note that this only works if filenames match those expected by wrath. Possible step options are: makewindows, getbarcodes, matrix, outliers or plot
  -v                verbose (only for the matrix generating step)
```


## Requirements

Command line programs:
- [samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)

Python:
- NumPy
- [Seaborn](https://seaborn.pydata.org/installing.html)
- [matplotlib](https://matplotlib.org/)
- [pandas](https://pandas.pydata.org/)
- [sklearn](https://scikit-learn.org/stable/index.html)

```bash
pip install -U numpy seaborn matplotlib pandas scikit-learn
```

R:
- [ggplot2](https://ggplot2.tidyverse.org/)
- [tidyr](https://tidyr.tidyverse.org/)
- [dplyr](https://dplyr.tidyverse.org/)
- [nlraa](https://github.com/femiguez/nlraa)

In R, run

```R
# The easiest way to get ggplot2, tidyr and dplyr is to install the whole tidyverse:
install.packages("tidyverse")

#and then nlraa
install.packages("nlraa")
```

## Input files

The input necessary is a genome file in fasta format and a list of the sample bam files that need to be analysed including their paths.

Like:
```
/home/samples/bams/group1/sample1.bam
/home/samples/bams/group1/sample2.bam
/home/samples/bams/group2/sample3.bam
/home/samples/bams/group2/sample4.bam
```


## Output
