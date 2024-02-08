#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=kingspeak
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N

cd /scratch/general/nfs1/utu_4310/whiptail_shared_data/trimmed_reads/
echo "file_name, total_reads" > count_reads.csv
cd merged_reads



for file in *.fq.gz; do
	total_lines=$(zcat "$file" | wc -l)


	read_counts=$((total_lines /4))


	echo "$file , $read_counts" >> ../count_reads.csv

cd ../unmerged_reads


for file in *.fq.gz; do
        total_lines=$(zcat "$file" | wc -l)


        read_counts=$((total_lines /4))

        echo "$file , $read_counts" >> ../count_reads.csv


done
