#!bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input files] [-g reference genome file] [-o output file where mapped reads go] [-t number of threds]
This program will map reads onto the reference directory
        -h show help text
        -i directory name where input files are located
        -g path and name of reference file
        -o output file where mapped reads go
	-a  A number between (0-1) indicating Viral Frequency for Freebayes variant calling (i.e.: 0.5 = 50% viral frequency).
	-t  Number of CPU processors
	"
options=':h:i:g:o:a:t:'
while getopts $options option; do
        case "$option" in
                h) echo "$usage"; exit;;
                i) i=$OPTARG;;
                g) g=$OPTARG;;
                o) o=$OPTARG;;
		a) a=$OPTARG;;
		t) t=$OPTARG;;
                :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
                \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
        esac
done

echo ""
echo "working directory: " $i
echo "reference genome: " $g
echo "output file: " $o

#################################################
### Performing Variant Calling with freebayes ###
#################################################

module load freebayes/1.3.4

echo "Performing Variant Calling with freebayes:"
echo ""

#x= ls -1 *.bam
freebayes -f ${g} ${i} > ${o}/RLK034_USD16091387L_HJNHCDSXX_L2.merged_sorted.bam.vcf
#x=${i}
#for x in *.bam; do freebayes-parallel <(fasta_generate_regions.py ${g}.fai 2000) ${t} -f ${g} -F ${a} -b ${x} > ${o}/RLK034_USD16091387L_HJNHCDSXX_L2.merged_sorted.bam.freebayes.vcf
#done

module unload freebayes/1.3.4


##################
### Count SNPS ###
##################

module load vcftools/0.1.15-6

vcfallelicprimitives ${o}/RLK034_USD16091387L_HJNHCDSXX_L2.merged_sorted.bam.vcf > output_primitives.vcf

vcftools --vcf output_primitives.vcf --keep-only-indels --count-alleles

############################
### Remove SNPs depth <5 ###
############################

vcffilter -f "DP >= 5" input.vcf > output_filtered.vcf

module unload vcftools/0.1.15-6
}

