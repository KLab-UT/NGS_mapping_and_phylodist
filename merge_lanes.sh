#!/bin/bash

{
usage="$(basename "$0") [-h] [-i <working_directory>] [-g reference genome directory] [-o <output_directory>] This program merges bam files that were sequenced on lanes and change the names of bam files sequences on single lanes
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
#merge (merged.bam) files sequenced across multiple lanes
# KLC098 and RLK019
# I use "[ -f <file_name> ] && do thing" to check if file exists before doing tasks
[ -f <KLC098_${g}_merged.bam> ] && samtools merge KLC098_${g}_merged.bam KLC098_*_L4.merged_sorted.bam KLC098_*_L1.merged_sorted.bam

[ -f <RLK019_${g}_merged.bam> ] && samtools merge RLK019_${g}_merged.bam RLK019_*_L3.merged_sorted.bam RLK019_*_L1.merged_sorted.bam

#merge (unmerged.bam) files that were sequenced across multiple lanes
[ -f <KLC098_${g}_unmerged.bam> ] && samtools merge KLC098_${g}_unmerged.bam KLC098_*_L4.unmerged_sorted.bam KLC098_*_L1.unmerged_sorted.bam

[ -f <RLK019_${g}_unmerged.bam> ] && samtools merge RLK019_${g}_unmerged.bam RLK019_*_L3.unmerged_sorted.bam RLK019_*_L1.unmerged_sorted.bam

#change names of (merged_sorted.bam) and (unmerged_sorted.bam) files sequenced on a single lane
# RLK004 and RLK034
[ -f <RLK004_${g}_merged.bam> ] && cp RLK004*.merged_sorted.bam RLK004_${g}_merged.bam
[ -f <RLK034_${g}_merged.bam> ] && cp RLK034*.merged_sorted.bam RLK034_${g}_merged.bam

[ -f <RLK004_${g}_unmerged.bam> ] && cp RLK004*.unmerged_sorted.bam RLK004_${g}_unmerged.bam
[ -f <RLK034_${g}_unmerged.bam> ] && cp RLK034*.unmerged_sorted.bam RLK034_${g}_unmerged.bam

#files look like this...
#RLK004_USD16091389L_HJNHCDSXX_L3.unmerged_sorted.bam
#RLK034_USD16091387L_HJNHCDSXX_L2.unmerged_sorted.bam

#module unload bwa/2020_03_19
module unload samtools/1.16
}
