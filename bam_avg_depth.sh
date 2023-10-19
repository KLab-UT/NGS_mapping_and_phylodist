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
	# get average depth
	avg_depth=$(samtools depth -a ${1}_merged.bam | awk '{sum+=$3} END {print sum/NR}')
	#sepeate sample_ID and Ref_name using IFS
	IFS=_ read sample_ID ref_name1 ref_name2 <<< ${1}
	ref_name="${ref_name1}_${ref_name2}"
	# get percentage of reads mapped to each reference
	#denominator=$(samtools view -c ${1}_merged.bam)
	#numerator=$(samtools view -c -F 260 ${1}_merged.bam)
	#percentage=$((numerator/denominator))
	#percentage=$(echo "scale=2; $numerator / $denominator * 100" | bc)

	#https://sarahpenir.github.io/bioinformatics/awk/calculating-mapping-stats-from-a-bam-file-using-samtools-and-awk/
	percentage=$(samtools flagstat ${1}_merged.bam | awk -F "[(|%]" 'NR == 3 {print $2}')
#	samtools view -F 0x4 foo.sorted.bam | cut -f 1 | sort | uniq | wc -l
	total_reads=$(samtools flagstat ${1}_merged.bam | awk -F " " 'NR == 1 {print $1}')
	#used commas as delimiters, could use spaces instead if prefered
	echo "$sample_ID,$ref_name,$total_reads,$avg_depth,$percentage" >> "$2"/depth_percentage.txt
}
export -f depth

echo "Reading depth."
cd $i
# in *.bam '*' is turned into the variable ${1}
# RLK004_Aspidoscelis_marmoratus_merged.bam
# [A-Z]\+ looks for 1 or more capital letters
######################################################### if not using parallel do "xargs -I {}" instead
#ls | grep '^[A-Z]\+[0-9]\+_[A-Za-z]\+_[a-z]\+_merged.bam' | cut -d "_" -f "1,2,3" | parallel depth "{}" "$o"

genome=$(ls | grep '^[A-Z]\+[0-9]\+_[A-Za-z]\+_[a-z]\+_merged.bam' | cut -d "_" -f "1,2,3")
#(-e enables interpretation of backslash escapes)
echo -e "Genome: $genome\nOutput: $o"

echo "$genome" | parallel depth "{}" "$o"

echo "Done"
module unload samtools/1.16
}

