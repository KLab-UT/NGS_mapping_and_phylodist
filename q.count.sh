#!/bin/sh
#SBATCH --account=utu-biol4310
#SBATCH --partition=lonepeak
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N


# create working environment
#wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data/cleaned_reads/unmerged_reads

#bash read1_read2_sort.sh -i $wd -o $wd
bash unmerged_read_count.sh -i $wd -o $wd
