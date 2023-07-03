#!/bin/sh
#SBATCH --account=utu-biol4310
#SBATCH --partition=lonepeak
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N


# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data/mapped_reads

#mkdir depth
for reference in */; do
	bash bam_avg_depth.sh -i $wd/$reference -o $wd/$reference
done 

