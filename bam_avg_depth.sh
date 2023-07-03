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
echo "working directory: " $i
echo "output file: " $o	


module load samtools/1.16

depth() {
	#clear file contents
	> "$2"/depth_percentage.txt

	# get average depth
	echo ${1} >> "$2"/depth_percentage.txt
	echo "average depth" >> "$2"/depth_percentage.txt
	samtools depth -a ${1}.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' >> "$2"/depth_percentage.txt
	# get percentage of reads mapped to each reference
	echo "percentage" >> "$2"/depth_percentage.txt
	echo "denominator" >> "$2"/depth_percentage.txt
	denominator=samtools view -c >> "$2"/depth_percentage.txt
	echo "numerator" >> "$2"/depth_percentage.txt
	numerator=samtools view -c -F 260 >> "$2"/depth_percentage.txt
	percentage=$numerator/$denominator
	echo $percentage >> "$2"/depth_percentage.txt
}
export -f depth

echo "Read depth."
cd $i
# in *.bam '*' is turned into the variable ${1}
ls *[0-9]_[a-z][a-z].bam | cut -d "." -f "1" | parallel depth {} $o

echo "Done"
module unload samtools/1.16
}

