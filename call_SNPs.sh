#!bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input files] [-g mapped reads directory]
This program will map reads onto the reference directory
        -h show help text
        -i directory name where input files are located
        -g mapped reads directory"

options=':h:i:g:'
while getopts $options option; do
        case "$option" in
                h) echo "$usage"; exit;;
                i) i=$OPTARG;;
                g) g=$OPTARG;;
                :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
                \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
        esac
done

echo ""
echo "working directory: " $i
echo "mapped reads directory: " $g
module load BCFtools/1.3.1
module load SAMtools/1.3.1

# List all BAM files in current directory
bam_files=$(ls *.bam | tr '\n' ' ')

# Create a space-separated string of bam files
BAM_FILES=$(echo $bam_files | tr ' ' '\n' | paste -s -d ' ')

# Run mpileup on all bam files and pipe the output to BCFtools
samtools mpileup -g -B $BAM_FILES | bcftools view -Ou - > all_samples.bcf
