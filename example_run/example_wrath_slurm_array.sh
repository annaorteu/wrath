#!/bin/bash

#SBATCH -J jaccard
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=15
#SBATCH -o slurm_jaccard_%j.out
#SBATCH -e slurm_jaccard_%j.err
#SBATCH -p cclake
#SBATCH -A JIGGINS-SL3-CPU
#SBATCH --array=1-195

#load the necessary modules
module load bedtools
module load python-3.6.1-gcc-5.4.0-23fr5u4 #for slurm
module load R/4.0.3 
source ~/.bashrc

#set variables
chr=$(sed -n "$SLURM_ARRAY_TASK_ID"p Hera_chr)
winSize=50000
threads=10
start=700000
end=850000
group=malleti.txt

#run wrath
wrath -g ~/genomes/Hmel/Hmel2.5.fa -c ${chr}  -w ${winSize}  -a ${group} -t ${threads} -l -s ${start} -e ${end} 