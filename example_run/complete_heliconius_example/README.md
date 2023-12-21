# Complete test example of Wrath

This folder contains all the code, data and results of a small test example using butterfly (Heliconius erato) data. 

Contents:
1. bams: _Heliconius erato lativitta_ bam files used for the example.\
2. samples.txt: list of sample names with paths.
3. code_line.txt: line of code to run
4. nohup.out: std-out output of the run
5. wrath_out: output of the example


Code to run: 

```{bash}
wrath -g Heliconius_erato_demophoon_v1_-_scaffolds.fa -c Herato0204 -w 10000 -a samples.txt -t 15 -l 
```

Note that if you run the example in this folder the `wrath_out` directory will be overwritten. Changing the name of the folder would prevent that. 
