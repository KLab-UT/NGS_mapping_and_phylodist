#!/bin/bash

{
usage="$(basename "$0") [-h] [-i <working_directory] [-o <output_directory>]This program will trim raw reads.
    -h  show this help text
    -i  Path to the working directory (the main directory for the raw reads)    -o  Path to the output directory (the main directory for the clean reads)"
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


# mandatory arguments
if [ ! "$i" ] || [ ! "$o" ]; then
  echo "arguments -o, and -i  must be provided"
  echo "$usage" >&2; exit 1
fi


#Get number of unmerged reads from fastp for each species
#a. Count number of unmerged for read1, then count number of unmerged for read2

cd $i

#A. septemvittatus 
    #KLC098_USD16091388L_HKFJFDSXX_L4_1.fq.gz
    #KLC098_USD16091388L_HKFJFDSXX_L4_2.fq.gz
    #KLC098_USD16091388L_HKG5MDSXX_L1_1.fq.gz
    #KLC098_USD16091388L_HKG5MDSXX_L1_2.fq.gz

#whole species read count
echo "KLC098 A.septemvittatus total read count (read1_count + read2_count)" >> unmerged_read_count.txt
$($(zcat KLC098*unmerged[12].fq.gz | wc -l) + $(zcat KLC098*unmerged[12].fq.gz | wc -l))/4|bc >> $o/unmerged_read_count.txt

#read1 count
echo "KLC098 read_1 count" >> unmerged_read_count.txt
$(zcat KLC098*unmerged1.fq.gz read1.txt | wc -l)4|bc >> $o unmerged_read_count.txt

#read2 count
echo "KLC098 read_1 count" >> unmerged_read_count.txt
$(zcat KLC098*unmerged2.fq.gz read1.txt | wc -l)4|bc >> $o unmerged_read_count.txt

# A. gularis
    #RLK004_USD16091389L_HJNHCDSXX_L3_1.fq.gz
    #RLK004_USD16091389L_HJNHCDSXX_L3_2.fq.gz


#whole species read count
echo "RLK004 A.gularis total read count (read1_count + read2_count)" >> unmerged_read_count.txt
$($(zcat RLK004*unmerged[12].fq.gz | wc -l) + $(zcat RLK004*unmerged[12].fq.gz | wc -l))/4|bc >> $o/unmerged_read_count.txt

#read1 count
echo "RLK004 read_1 count" >> unmerged_read_count.txt
$(zcat RLK004*unmerged1.fq.gz read1.txt | wc -l)4|bc >> $o/unmerged_read_count.txt

#read2 count
echo "RLK004 read_1 count" >> unmerged_read_count.txt
$(zcat RLK004*unmerged2.fq.gz read1.txt | wc -l)4|bc >> $o/unmerged_read_count.txt

# A. inornatus
    #RLK019_USD16091390L_HJMVKDSXX_L3_1.fq.gz
    #RLK019_USD16091390L_HJMVKDSXX_L3_2.fq.gz
    #RLK019_USD16091390L_HKG5MDSXX_L1_1.fq.gz
    #RLK019_USD16091390L_HKG5MDSXX_L1_2.fq.gz

#whole species read count
echo "RLK019 A.inornatus total read count (read1_count + read2_count)" >> unmerged_read_count.txt
$($(zcat RLK019*unmerged[12].fq.gz | wc -l) + $(zcat RLK019*unmerged[12].fq.gz | wc -l))/4|bc >> $o/unmerged_read_count.txt

#read1 count
echo "RLK019 read_1 count" >> unmerged_read_count.txt
$(zcat RLK019*unmerged1.fq.gz read1.txt | wc -l)4|bc >> $o/unmerged_read_count.txt

#read2 count
echo "RLK019 read_1 count" >> unmerged_read_count.txt
$(zcat RLK019*unmerged2.fq.gz read1.txt | wc -l)4|bc >> $o/unmerged_read_count.txt

# A. marmoratus
    #RLK034_USD16091387L_HJNHCDSXX_L2_1.fq.gz
    #RLK034_USD16091387L_HJNHCDSXX_L2_2.fq.gz

#whole species read count
echo "RLK034 A.marmoratus total read count (read1_count + read2_count)" >> unmerged_read_count.txt
$($(zcat RLK034*unmerged[12].fq.gz | wc -l) + $(zcat RLK034*unmerged[12].fq.gz | wc -l))/4|bc >> $o/unmerged_read_count.txt

#read1 count
echo "RLK034 read_1 count" >> unmerged_read_count.txt
$(zcat RLK034*unmerged1.fq.gz read1.txt | wc -l)4|bc >> $o/unmerged_read_count.txt

#read2 count
echo "RLK034 read_1 count" >> unmerged_read_count.txt
$(zcat RLK034*unmerged2.fq.gz read1.txt | wc -l)4|bc >> $o/unmerged_read_count.txt

} | tee outfile 
