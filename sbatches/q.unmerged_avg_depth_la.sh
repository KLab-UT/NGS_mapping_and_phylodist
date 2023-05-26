#!/bin/sh
#SBATCH --account=utu-biol4310
#SBATCH --partition=lonepeak
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N


wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data/mapped_reads/la 
# Make working directory and run bash script
bash ~/Biol_4310/whiptail_nmt_variation/parallel_bam.sh -i $wd -o $wd/depth/
