#!/bin/bash

time parallel -j{}
{
usage="$(basename "$0") [-h] [-i <working_directory] [-o <output_directory>]
This program will trim raw reads. Raw reads are compressed fastq files (.fq.gz)
    -h  show this help text
    -i  Path to the working directory (the main directory for the raw reads)
    -o	Path to the output directory (the main directory for the clean reads)
    -t  Number of CPU processors"
options=':h:t:i:o:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    t) t=$OPTARG;;
    i) i=$OPTARG;;
    o) o=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done


echo ""
echo "Working directory: " $i
echo "Threads:           " $t
echo "Output directory: " $o
echo ""


# mandatory arguments
if [ ! "$t" ] || [ ! "$i" ] || [ ! "$o" ]; then
  echo "arguments -t, -o, and -i  must be provided"
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
cd ${i}
pwd
#ls *.fq.gz | cut -d "_" -f "1,2,3,4" | sort | uniq | time parallel -j${t} --eta --bar "fastp -w 2 -i {}_1.fq.gz -I {}_2.fq.gz -m --merged_out ${o}/merged_reads/{}.merged.fq --out1 ${o}/unmerged_reads/{}.unmerged1.fq --out2 ${o}/unmerged_reads/{}.unmerged2.fq"
ls *.fq.gz | cut -d "_" -f "1,2,3,4" | sort | uniq > fastq_list.txt
while read sample; do
	fastp 	-w ${t} \
		-i ${sample}_1.fq.gz -I ${sample}_2.fq.gz \
		-m --merged_out ${o}/merged_reads/"$sample".merged.fq / --out1 ${o}/unmerged_reads/"$sample".unmerged1.fq --out2 ${o}/unmerged_reads/"$sample".unmerged2.fq \

cd ../trimmed_reads/merged_reads
gzip "$sample".merged.fq
cd ../unmerged_reads
gzip "$sample".unmerged1.fq
gzip "$sample".unmerged2.fq

cd ../../raw_reads

done<fastq_list.txt
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
