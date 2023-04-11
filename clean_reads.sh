#!/bin/bash

{
usage="$(basename "$0") [-h] [-d <working_directory] [-o <output_directory>]
This program will trim raw reads.
    -h  show this help text
    -d  Path to the working directory (the main directory for the raw reads)
    -o	Path to the output directory (the main directory for the clean reads)
    -t  Number of CPU processors"
options=':h:t:d:o:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    t) t=$OPTARG;;
    d) d=$OPTARG;;
    o) o=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done


echo ""
echo "Working directory: " $d
echo "Threads:           " $t
echo "Output directory: " $o
echo ""


# mandatory arguments
if [ ! "$t" ] || [ ! "$d" ] || [ ! "$o" ]; then
  echo "arguments -t, -o, and -d  must be provided"
  echo "$usage" >&2; exit 1
fi

begin=`date +%s`

####################################################
# Trimming downloaded Illumina datasets with fastp #
####################################################
module load fastp/0.20.1

echo "Trimming downloaded datasets with fastp."
echo ""

function trim_reads {
cd ${d}
pwd
ls *.fq.gz | cut -d "_" -f "1,2,3,4" | sort | uniq > fastq_list
while read sample; do 
	fastp 	-w ${t} \
		-i ${sample}_1.fq.gz -I ${sample}_2.fq.gz \
		-m --merged_out ${o}/merged_reads/"$sample"_merged.fq / --out1 ${o}/unmerged_reads/"$sample"_unmerged1.fq --out2 ${o}/unmerged_reads/"$sample"_unmerged2.fq \

cd ../cleaned_reads/merged_reads
gzip "$sample".merged.fq
cd ../cleaned_reads/unmerged1.fq
gzip "$sample".unmerged1.fq
cd ../cleaned_reads/unmerged2.fq
gzip "$sample".unmerged2.fq

cd ../../raw_reads

done<fastq_list
cd ..
}

trim_reads

module unload fastp/0.20.1

########## Done ###########

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

mdate=`date +'%d/%m/%Y %H:%M:%S'`

} | tee outfile
