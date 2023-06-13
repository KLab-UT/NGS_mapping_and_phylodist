#!bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input files] [-g mapped reads directory][-o output file where mapped reads go] [-t number of threds]
This program will map reads onto the reference directory
        -h show help text
        -i directory name where input files are located
        -o output file where mapped reads go
	-g mapped reads directory"

options=':h:i:g:o:'
while getopts $options option; do
        case "$option" in
                h) echo "$usage"; exit;;
                i) i=$OPTARG;;
		g) g=$OPTARG;;
		o) o=$OPTARG;;
                :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
                \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
        esac
done

echo ""
echo "working directory: " $i
echo "output file: " $o
echo "mapped reads directory: " $g

module load samtools/1.16

depth() {
        samtools depth -a ${1}.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > "$2"/${1}.depth.txt

}
export -f depth

echo "Read depth."
cd $i
ls *_${g}_* | cut -d "." -f "1" | parallel depth {} $o

echo "Done"
module unload samtools/1.16
}
