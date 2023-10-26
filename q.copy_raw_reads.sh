#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N

module load parallel

# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data

# download reference genomes: Randy
cd $wd/raw_reads

ls *fq.gz | parallel -j 10 cp {} /scratch/general/nfs1/utu_4310/whiptail_nmt_variation_data/raw_reads
