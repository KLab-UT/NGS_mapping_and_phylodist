#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=kingspeak
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N


# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data/mapped_reads
home=~/Biol_4310/NGS_mapping_and_phylodist
cd $wd
#mkdir depth
echo "#sample_ID, Ref_name, merge_status, total_reads, avg_depth, aligned_percentage, #_of_mapped_reads, mapped_percentage" > depth_percentage.txt
# */ is list of directories in current directory
for reference in */; do
	# -o is $wd because I want one text file for all reference diretories.
	bash $home/bam_avg_depth.sh -o $wd -i $wd/$reference
done 

