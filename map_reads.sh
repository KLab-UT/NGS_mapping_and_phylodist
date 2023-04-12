#!bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input files] [-g reference genome file] [-o output file where mapped reads go] [-t 
This program will map reads onto the reference directory
	-h show help text
	-i directory name where input files are located
	-g name of reference file
	-o output file where mapped reads go
	-t number of threads"
options=':h:i:r:p:d:'
while getopts $options option; do
	case "$option" in
		h) echo "$usage"; exit;;
		i) i=$OPTARG;;
		g) g=$OPTARG;;
		o) o=$OPTARG;;
		t) t=$optarg;;
		:) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
		\?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
	esac
done

echo ""
echo "working directory: " $i
echo "reference genome: " $g
echo "output file: " $o
echo "number of threads: " $t

# mandatory arguements
if [ ! "$i" ] || [ ! "$g" ] || [ ! "$o" ] || [ ! "$t" ]; then
	echo "arguments -i, -R, -o, and -t  must be provided"
	echo "$usage" >&2; exit 1
fi

# map reads --------------------------------------------------------------------------------

######################
# Indexing Reference #
######################

module load bwa/2020_03_19

echo "Indexing Reference"
echo ""
bwa index ${g}

module unload bwa/2020_03_19

###########################################################################################
# Aligning datasets againts reference with minimap #
###########################################################################################

module load bwa/2020_03_19
module load samtools/1.16

echo "Aligning reads with reference with bwa mem."
echo ""

#make file of file names the run through then
cd ${i}
ls *fq.gz | cut -d '.' -f '1' > cleaned_reads.txt

while read sample; do 
	bwa mem ${sample}.fq.gz > ${o}/${sample}_mapped.sam
	samtools sort ${o}/${sample}_mapped.sam > ${o}/${sample}_sorted.bam
done<cleaned_reads.txt

module unload bwa/2020_03_19
module unload samtools/1.16

}
