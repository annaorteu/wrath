#!/bin/bash
# ------------------------------------------------------------------
#         Anna Farre Orteu, 2022
#             af658@cam.ac.uk
#     Script to detect SVs in haplotagging data
# ------------------------------------------------------------------

bold=$(tput bold)
normal=$(tput sgr0)

subject=wrath
usage="
${bold}wrath: wrapped analysis of tagged haplotypes

${bold}DESCRIPTION:
${normal} Program produces a jaccard matrix camparing the barcode content between all pairs windows whithin a chromosome.

wrath [-h] [-g GENOMEFILE] [-c CHROMOSOMENAME] [-w WINDOWSIZE] [-s FILELIST] [-t THERADS] [-p] [-v] [-x STEP]

${bold}OPTIONS: ${normal}
  -h                show this help text
  -g GENOMEFILE     reference genome
  -c CHROMOSOMENAME chromosome
  -w WINDOWSIZE     window size
  -s FILELIST       list of bam files with paths of the individuals of the population/phenotype of interest
  -t THERADS        threads to use
  -p                skip plotting the heatmap
  -x STEP           start from a given step. Note that this only works if filenames match those expected by wrath. Possible step options are: makewindows, getbarcodes, matrix, outliers or plot
  -v                verbose (only for the matrix generating step)

"


# --- Option processing --------------------------------------------
if [ $# == 0 ] ; then
    echo "$usage"
    exit 1;
fi
while getopts "g:c:w:s:t:pvx:h" optname
  do
    case "$optname" in
      "g") genome="$OPTARG" ;;
      "c") chromosome="$OPTARG" ;;
      "w") winSize="$OPTARG" ;;
      "s") group="$OPTARG" ;;
      "t") threads="$OPTARG" ;;
      "p") noplot=1 ;;
      "v") verbose="--verbose" ;;
      "x") step="$OPTARG" ;;
      "h")
        echo "$usage"
        exit 0;
        ;;
      "?")
        >&2 echo "Unknown option $OPTARG"
        exit 1;
        ;;
      ":")
        >&2 echo "No argument value for option $OPTARG"
        exit 1;
        ;;
      *)
        >&2 echo "Unknown error while processing options"
        exit 1;
        ;;
    esac
  done

if [ $OPTIND -eq 1 ]; then
  printf "\nNo options were passed
  $usage"
  exit 1;
fi

shift $(($OPTIND - 1))


if [ ! -z "$step" ]
then
    case $step in
        makewindows) makewindows=1 ;;
        getbarcodes) getbarcodes=1 ;;
        matrix) matrix=1 ;;
        outliers) outliers=1 ;;
        plot) plot=1 ;;
        *) echo "Wrong step specified!"
        echo "$usage"
        exit 1
    esac
fi

# -----------------------------------------------------------------
#  SCRIPT LOGIC GOES HERE
# -----------------------------------------------------------------

echo "running wrath_out with options:

genome = ${genome}
chromosome = ${chromosome}
window size = ${winSize}
sample bams file = ${group}
threads = ${threads}

"
######################################################################
# Set working paths

#Get the directory where the scipt is saved
#the code resolves symlinks
SOURCE=${BASH_SOURCE[0]}
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )


######################################################################
# Modules
#create the output directory (if it doesn't exist)
mkdir -p wrath_out


######################################################################
# Make genomic windows

if [ -z ${step+x} ] || [ ! -z ${makewins+x} ]; then

  mkdir -p wrath_out/beds
  #check if the file with genome sizes already exists
  if [ ! -f wrath_out/size.genome ]; then
    # first get chromosome sizes from the reference
    echo "Getting chromsome sizes"
    faidx ${genome} -i chromsizes > wrath_out/size.genome || { >&2 echo 'Getting chromosome sizes from genome file failed' ; exit 1; }
  fi

  #then make windows
  echo "Getting ${chromosome} size from genome file"
  awk -v pat="${chromosome}" '$1 == pat {print $0}' wrath_out/size.genome > wrath_out/size.${chromosome} || { >&2 echo "Getting ${chromosome} size from genome file failed" ; exit 1; }

  echo "Making ${winSize} windows of ${chromosome}"
  bedtools makewindows -g wrath_out/size.${chromosome} -w ${winSize} >  wrath_out/beds/windows_${winSize}_${chromosome}.bed || { >&2 echo "Making ${winSize} windows of ${chromosome} failed" ; exit 1; } #need to split the bed file in multiple beds one per chr

  rm wrath_out/size.${chromosome}
fi


######################################################################
# Get barcodes

if [ -z ${step+x} ] || [ ! -z ${getbarcodes+x} ]; then

  #get barcodes by phenotype
  for sample in $(cat ${group})
  do
    echo "Getting ${sample} barcodes from ${chromosome}"
    samtools view -q 1 -@ ${threads} ${sample} ${chromosome} | grep -o -P "${chromosome}.*BX:Z:[0-9A-Z]*\t" | awk '{print $1"\t"$2"\t"$2"\t"$NF}' > wrath_out/beds/barcodes_${chromosome}_$(basename "$group" .txt)_$(basename $sample .bam).bed ; done || { >&2 echo "Getting ${sample} barcodes from ${chromosome} failed" ; exit 1; }

  #sort barcodes
  echo "Sorting of $(basename "$group" .txt) barcodes bed files from ${chromosome}"
  cat wrath_out/beds/barcodes_${chromosome}_$(basename "$group" .txt)_*bed | bedtools sort -i - | bgzip -@ ${threads} > wrath_out/beds/barcodes_${chromosome}_sorted_$(basename "$group" .txt).bed.gz && tabix wrath_out/beds/barcodes_${chromosome}_sorted_$(basename "$group" .txt).bed.gz && rm wrath_out/beds/barcodes_${chromosome}_$(basename "$group" .txt)_*bed || { >&2 echo "Sorting of $(basename "$group" .txt) barcodes bed files from ${chromosome} failed" ; exit 1; }

fi


######################################################################
# Generate similarity matrix

if [ -z ${step+x} ] || [ ! -z ${matrix+x} ]; then

  mkdir -p wrath_out/matrices
  # compute the jaccard index and save it in a matrix
  echo "Computing of jaccard index matrix for chromsome ${chromosome} of $(basename "$group" .txt) of window size ${winSize}"
  python ${DIR}/jaccard_matrix_simplequeue.py  --threads ${threads} -w wrath_out/beds/windows_${winSize}_${chromosome}.bed -b wrath_out/beds/barcodes_${chromosome}_sorted_$(basename "$group" .txt).bed.gz  -o wrath_out/matrices/jaccard_matrix_${winSize}_${chromosome}_$(basename "$group" .txt).txt ${verbose} || { >&2 echo  "Computing of jaccard index matrix for chromsome ${chromosome} of $(basename "$group" .txt) of window size ${winSize} failed" ; exit 1; }

  #edit the output (remove a colon at the end of the line)
  echo "Editing of jacard index matrix for chromsome ${chromosome} of $(basename "$group" .txt) of window size ${winSize}"
  sed -i 's/,$//' wrath_out/matrices/jaccard_matrix_${winSize}_${chromosome}_$(basename "$group" .txt).txt || { >&2 echo "Editing of jacard index matrix for chromsome ${chromosome} of $(basename "$group" .txt) of window size ${winSize} failed" ; exit 1; }

fi


######################################################################
# Detect outliers

if [ -z ${step+x} ] || [ ! -z ${outliers+x} ]; then

  echo "Detecting outliers"
  mkdir -p wrath_out/outliers
  Rscript ${DIR}/outlier_detection.R wrath_out/matrices/jaccard_matrix_${winSize}_${chromosome}_$(basename "$group" .txt).txt wrath_out/outliers/outliers_${winSize}_${chromosome}_$(basename "$group" .txt).csv || { >&2 "Detecting outliers from matrix wrath_out/matrices/jaccard_matrix_${winSize}_${chromosome}_$(basename "$group" .txt).txt step failed"; exit 1; }

fi


######################################################################
# Detect SVs and plot results

#if the option is given to plot it, then do
if [ -z ${step+x} ] || [ ! -z ${plot+x} ]; then # -z asks if ${plot+x} is empty. Thus, [ ! -z ${plot+x} ] asks if ${plot+x} is not empty

  #plot the optput
  mkdir -p wrath_out/plots
  python ${DIR}/sv_detection_and_heatmap.py --matrix wrath_out/matrices/jaccard_matrix_${winSize}_${chromosome}_$(basename "$group" .txt).txt -o wrath_out/outliers/outliers_${winSize}_${chromosome}_$(basename "$group" .txt).csv  -p wrath_out/plots/heatmap_${winSize}_${chromosome}_$(basename "$group" .txt).png -s wrath_out/SVs/sv_${winSize}_${chromosome}_$(basename "$group" .txt).txt -w ${winSize} || { >&2 "Detecting SVs and plotting of matrix wrath_out/matrices/jaccard_matrix_${winSize}_${chromosome}_$(basename "$group" .txt).txt step failed"; exit 1; }

fi


######################################################################
# Detect SVs without plotting the results

#if the option is given to plot it, then do
if [ -z ${step+x} ] || [ ! -z ${plot+x} ]; then # -z asks if ${plot+x} is empty. Thus, [ ! -z ${plot+x} ] asks if ${plot+x} is not empty

  #plot the optput
  mkdir -p wrath_out/plots
  mkdir -p wrath_out/SVs
  python ${DIR}/sv_detection.py --matrix wrath_out/matrices/jaccard_matrix_${winSize}_${chromosome}_$(basename "$group" .txt).txt -o wrath_out/outliers/outliers_${winSize}_${chromosome}_$(basename "$group" .txt).csv -s wrath_out/SVs/sv_${winSize}_${chromosome}_$(basename "$group" .txt).txt -w ${winSize} || { >&2 "Detecting SVs in matrix wrath_out/matrices/jaccard_matrix_${winSize}_${chromosome}_$(basename "$group" .txt).txt step failed"; exit 1; }

fi