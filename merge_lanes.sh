#!/bin/bash

{
usage="$(basename "$0") [-h] [-i <working_directory>] [-g reference genome directory] [-o <output_directory>] This program merges bam files that represented the same reads in different lanes
    -h  show this help text
    -i  Path to the working directory (the main directory for the raw reads)    
    -o  Path to the output directory (the main directory for the clean reads)   
    -g is the name of the reference directory"
options=':h:i:g:o:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    i) i=$OPTARG;;
    g) g=$OPTARG;;
    o) o=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done


echo ""
echo "Working directory: " $i
echo "Output directory: " $o
echo "Reference directory: " $g
echo ""


# mandatory arguments
if [ ! "$i" ] || [ ! "$g" ] || [ ! "$o" ]; then
  echo "arguments -o, and -i  must be provided"
  echo "$usage" >&2; exit 1
fi

#module load bwa/2020_03_19
module load samtools/1.16

cd $i

samtools merge KLC098_${g}_merged.bam KLC098_USD16091388L_HKFJFDSXX_L4.merged_sorted.bam KLC098_USD16091388L_HKG5MDSXX_L1.merged_sorted.bam

samtools merge RLK019_${g}_merged.bam RLK019_USD16091390L_HJMVKDSXX_L3.merged_sorted.bam RLK019_USD16091390L_HKG5MDSXX_L1.merged_sorted.bam

#module unload bwa/2020_03_19
module unload samtools/1.16
}
