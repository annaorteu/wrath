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

chr=`cat Hera_chr |awk 'NR == '$SLURM_ARRAY_TASK_ID' {print $1}'`

module load python-3.6.1-gcc-5.4.0-23fr5u4

~/wrath/wrath -g ~/genomes/Hera/Heliconius_erato_demophoon_v1_-_scaffolds.fa -c ${chr}  -w 10000  -s all_erato.txt -t 15 

