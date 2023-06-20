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

module load samtools/1.16

cd ${i}

ls *_${g}_* | cut -d "_" -f "1,2" | uniq > id.txt

for id in $(cat "id.txt"); do
	echo "${id}" >> ${i}/mapping_percentage.txt

	#merged
	a=$(samtools view -c ${id}_merged.bam) 
	b=$(samtools view -c -F 260 ${id}_merged.bam) 
	echo "Merged Mapped Percentage: " >> ${i}/mapping_percentage.txt
	echo "scale=4; ${b} / ${a}" | bc >> ${i}/mapping_percentage.txt
	echo " " >> ${i}/mapping_percentage.txt

	#unmerged
	c=$(samtools view -c ${id}_unmerged.bam) 
	d=$(samtools view -c -F 260 ${id}_unmerged.bam) 
	echo "Unmerged Mapped Percentage: " >> ${i}/mapping_percentage.txt
	echo "scale=4; ${d} / ${c}" | bc >> ${i}/mapping_percentage.txt
	echo " " >> ${i}/mapping_percentage.txt

	#total
	e=$(samtools view -c ${id}) 
	f=$(samtools view -c -F 260 ${id})
	echo "Merged and Unmerged Mapped Percentage: " >> ${i}/mapping_percentage.txt
	echo "scale=4; ${f} / ${e}" | bc >> ${i}/mapping_percentage.txt
	echo " " >> ${i}/mapping_percentage.txt
	echo " " >> ${i}/mapping_percentage.txt

done

module unload samtools/1.16
}


