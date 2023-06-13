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

bash bam_avg_depth.sh -i $wd/hc -o $wd/hc/depth

bash bam_avg_depth.sh -i $wd/la -o $wd/la/depth

bash bam_avg_depth.sh -i $wd/pb -o $wd/pb/depth

bash bam_avg_depth.sh -i $wd/pm -o $wd/pm/depth

bash bam_avg_depth.sh -i $wd/pr -o $wd/pr/depth

bash bam_avg_depth.sh -i $wd/su -o $wd/su/depth
