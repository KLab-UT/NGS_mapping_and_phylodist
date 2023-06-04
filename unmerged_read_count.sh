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
echo "KLC098 A. septemvittatus total read count (read1_count + read2_count)" >> unmerged_read_count_working.txt
zcat KLC098*fq.gz | wc -l >> unmerged_read_count_working.txt

#read1 count
echo "KLC098 read_1 count" >> unmerged_read_count_working.txt
while IFS= read -r line || [[ -n "$line" ]]; do
	if [[ $line == KLC098* ]]; then
		zcat "$line" | wc -l >> unmerged_read_count_working.txt
	fi
done < read1.txt

#read2 count
echo "KLC098 read_2 count" >> unmerged_read_count_working.txt
while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $line == KLC098* ]]; then
                zcat "$line" | wc -l >> unmerged_read_count_working.txt
        fi
done < read2.txt

echo "" >> unmerged_read_count_working.txt



# A. gularis
    #RLK004_USD16091389L_HJNHCDSXX_L3_1.fq.gz
    #RLK004_USD16091389L_HJNHCDSXX_L3_2.fq.gz

#whole species read count
echo "RLK004 A.gularis total read count (read1_count + read2_count)" >> unmerged_read_count_working.txt
zcat RLK004*fq.gz | wc -l >> unmerged_read_count_working.txt

#read1 count
echo "RLK004 read_1 count" >> unmerged_read_count_working.txt
while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $line == RLK004* ]]; then
                zcat "$line" | wc -l >> unmerged_read_count_working.txt
        fi
done < read1.txt

#read2 count
echo "RLK004 read_2 count" >> unmerged_read_count_working.txt
while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $line == RLK004* ]]; then
                zcat "$line" | wc -l >> unmerged_read_count_working.txt
        fi
done < read2.txt


echo "" >> unmerged_read_count_working.txt



# A. inornatus
    #RLK019_USD16091390L_HJMVKDSXX_L3_1.fq.gz
    #RLK019_USD16091390L_HJMVKDSXX_L3_2.fq.gz
    #RLK019_USD16091390L_HKG5MDSXX_L1_1.fq.gz
    #RLK019_USD16091390L_HKG5MDSXX_L1_2.fq.gz

#whole species read count
echo "RLK019 A.inornatus total read count (read1_count + read2_count)" >> unmerged_read_count_working.txt
zcat RLK019*fq.gz | wc -l >> unmerged_read_count_working.txt

#read1 count
echo "RLK019 read_1 count" >> unmerged_read_count_working.txt
while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $line == RLK019* ]]; then
                zcat "$line" | wc -l >> unmerged_read_count_working.txt
        fi
done < read1.txt

#read2 count
echo "RLK019 read_2 count" >> unmerged_read_count_working.txt
while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $line == RLK019* ]]; then
                zcat "$line" | wc -l >> unmerged_read_count_working.txt
        fi
done < read2.txt

echo "" >> unmerged_read_count_working.txt



# A. marmoratus
    #RLK034_USD16091387L_HJNHCDSXX_L2_1.fq.gz
    #RLK034_USD16091387L_HJNHCDSXX_L2_2.fq.gz

#whole species read count
echo "RLK034 A.marmoratus total read count (read1_count + read2_count)" >> unmerged_read_count_working.txt
zcat RLK034*fq.gz | wc -l >> unmerged_read_count_working.txt

#read1 count
echo "RLK034 read_1 count" >> unmerged_read_count_working.txt
while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $line == RLK034* ]]; then
                zcat "$line" | wc -l >> unmerged_read_count_working.txt
        fi
done < read1.txt

#read2 count
echo "RLK034 read_2 count" >> unmerged_read_count_working.txt
while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $line == RLK034* ]]; then
                zcat "$line" | wc -l >> unmerged_read_count_working.txt
        fi
done < read2.txt

} | tee outfile
