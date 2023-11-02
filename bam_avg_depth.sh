#!bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input files] [-o output file where mapped reads go] [-t number of threds]
This program will map reads onto the reference directory
        -h show help text
        -i directory name where input files are located
        -o output file where mapped reads go"

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
echo "working directory:  $i"
echo "output file:  $o"

if [ ! "$i" ] || [ ! "$o" ]; then
	echo "arguments -o, and -i  must be provided"
	echo "$usage" >&2; exit 1
fi

module load samtools/1.16

# sample_ID, Ref_name, #_of_reads, avg_depth, map_percentage
# make depth_percentage.txt in sbatch file q.avg_depth.sh
# echo "#sample_ID, Ref_name,           number_of_reads, avg_depth, map_percentage" > "$2"/depth_percentage.txt

depth() {
    g=$1
    output=$2
    echo -e "sample: ${g}\n"
	# get average depth
	avg_depth=$(samtools depth -a "${g}" | awk '{sum+=$3} END {print sum/NR}')
    echo "avg_depth: $avg_depth"
	#sepeate sample_ID and Ref_name using IFS
	IFS="_." read sample_ID ref_name1 ref_name2 merge_status <<< ${g}
	ref_name="${ref_name1}_${ref_name2}"
	# get percentage of reads mapped to each reference
	#denominator=$(samtools view -c ${g}_merged.bam)
	#numerator=$(samtools view -c -F 260 ${g}_merged.bam)
	#percentage=$((numerator/denominator))
	#percentage=$(echo "scale=2; $numerator / $denominator * 100" | bc)

	#https://sarahpenir.github.io/bioinformatics/awk/calculating-mapping-stats-from-a-bam-file-using-samtools-and-awk/
	aligned_percentage=$(samtools flagstat "${g}" | awk -F "[(|%]" 'NR == 3 {print $2}')
    # map_recentage is number of alignments divided by total number of reads, but mapped_reads is number of reads.
    # 0x4    UNMAP   segment unmapped
    # 0x100 SECONDARY   secondary alignment
    # 0x800   SUPPLEMENTARY   supplementary alignment
    mapped_reads=$(samtools view -F 0x4 ${g} | cut -f 1 | sort | uniq | wc -l)
    unmapped_reads=$(samtools view -f 0x4 ${g} | cut -f 1 | sort | uniq | wc -l)
	#total_reads=$(samtools flagstat "${g}" | awk -F " " 'NR == 1 {print $1}')
    total_reads=$(($mapped_reads+$unmapped_reads)) 
    mapped_percentage=$(echo "scale=6; $mapped_reads / $total_reads" | bc)
    echo "Aligned_percentage: $aligned_percentage" 
    echo "mapped_reads: $mapped_reads"
    echo "unmapped_reads: $unmapped_reads"
    echo "total_reads $total_reads"
    #mapping_fraction=$(($mapped_reads / $total_reads))
	#used commas as delimiters, could use spaces instead if prefered
	echo "$sample_ID,$ref_name,$merge_status,$total_reads,$avg_depth,$aligned_percentage,$mapped_reads,$mapped_percentage" >> ${output}/depth_percentage.txt
}
export -f depth

echo "Reading depth."
# in *.bam '*' is turned into the variable ${1}
# RLK004_Aspidoscelis_marmoratus_merged.bam
# [A-Z]\+ looks for 1 or more capital letters
######################################################### if not using parallel do "xargs -I {}" instead
#ls | grep '^[A-Z]\+[0-9]\+_[A-Za-z]\+_[a-z]\+_merged.bam' | cut -d "_" -f "1,2,3" | parallel depth "{}" "$o"

# files look like this 
# KLC098_Alligator_mississippiensis_merged.bam
# KLC098_Alligator_mississippiensis_unmerged.bam
cd $i
genome=$(ls | grep -E '^[A-Z]+[0-9]+_[A-Za-z]+_[a-z]+_[a-z]+.bam' | cut -d "_" -f "1-4")
#(-e enables interpretation of backslash escapes)
echo -e "Genome: $genome\nOutput: $o"

# make header
echo "#sample_ID, Ref_name, merge_status, total_reads, avg_depth, aligned_percentage, #_of_mapped_reads, mapped_percentage" > mapped_percentage.txt
# get info from .bam files
echo "$genome" | parallel depth "{}" "$o"

echo "Done"
module unload samtools/1.16
}

