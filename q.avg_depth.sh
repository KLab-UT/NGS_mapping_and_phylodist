#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N


# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data/mapped_reads

#mkdir depth
# */ is list of directories in current directory
for reference in */; do
	# -o is $wd because I want one text file for all reference diretories.
	bash bam_avg_depth.sh -i $wd/$reference -o $wd
done 

