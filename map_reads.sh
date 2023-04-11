#!bin/bash

{
usage="$(basename "$0") [-h] [-i directory for input files] [-r reference genome file] [-p path to reference genome] [-d working directory]
This program will map reads onto the reference directory
	-h show help text
	-i directory name where input files are located
	-r name of reference file
	-p path to reference file
	-d working directory"
options=':h:i:r:p:d:'
while getopts $options option; do
	case "$option" in
		h) echo "$usage"; exit;;
		i) i=$OPTARG;;
		r) r=$OPTARG;;
		p) p=$OPTARG;;
		d) d=$OPTARG;;
		:) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
		\?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
	esac
done

echo ""
echo "input directory: " $i
echo "reference genome: " $r
echo "path to reference genome: " $p
echo "working directory: " $d

# mandatory arguements
if [ ! "$I" ] || [ ! "$N" ] || [ ! "$g" ] || [ ! "$t" ] || [ ! "$d" ]; then
	echo "arguments -I, -N, -g, -t, and -d  must be provided"
	echo "$usage" >&2; exit 1
fi

# map reads --------------------------------------------------------------------------------

######################
# Indexing Reference #
######################

module load bwa/2020_03_19

echo "Indexing Reference"
echo ""
cd ${d}
bwa index ${d}/${p}/${r}

module unload bwa/2020_03_19

###########################################################################################
# Aligning datasets againts reference with minimap #
###########################################################################################

module load bwa/2020_03_19
module load samtools/1.16

echo "Aligning reads with reference with bwa mem."
echo ""

cd ${d}
while read cleaned_file; do 
	do bwa mem ${cleaned_file} > ${d}/mapped_reads/${cleaned_file}_mapped.sam
	samtools sort ${d}/mapped_reads/${cleaned_file}_mapped.sam > ${d}/mapped_reads/${cleaned_file}_sorted.bam
done<${i}

module unload bwa/2020_03_19
module unload samtools/1.16

########## Done ###########

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

mdate=`date +'%d/%m/%Y %H:%M:%S'`
# NOTE: If you are running a mac and having trouble with the code below,
# ----- try using 'vm_stat' instead of 'vmstat'
mcpu=$[100-$(vmstat 1 2|tail -1|awk '{print $15}')]%
mmem=`free | grep Mem | awk '{print $3/$2 * 100.0}'`
echo "$mdate | $mcpu | $mmem" >> ./stats-cpu
###############################################################
#
| tee logfile
}
