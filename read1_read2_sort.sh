#!/bin/bash

{
usage="$(basename "$0") [-h] [-i <working_directory] [-o <output_directory>] This program will trim raw reads.
    -h  show this help text
    -i  Path to the working directory (the main directory for the raw reads)
    -o  Path to the output directory (the main directory for the clean reads)"
options=':h:i:o:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    i) i=$OPTARG;;
    o) o=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

echo ""
echo "Working directory: " $i
echo "Output directory: " $o
echo ""

# Directory containing the files
# Get a list of files in the directory
files=$(ls "$i")

# Iterate over the files
for file in $files; do
    # Extract the unmerged number using cut
    unmerged=$(echo "$file" | cut -d '.' -f 2 )

    if [ "$unmerged" = "unmerged1" ]; then
        echo "$file" >> $o/read1.txt
    elif  [ "$unmerged" = "unmerged2" ]; then
        echo "$file" >> $o/read2.txt
    else
	echo "neither read1 or read2"
    fi
done

} | tee outfile
