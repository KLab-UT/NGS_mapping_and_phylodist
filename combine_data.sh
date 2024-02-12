#!/bin/sh
##SBATCH --account=utu
##SBATCH --partition=lonepeak
##SBATCH --nodes=1
##SBATCH --ntasks=1
##SBATCH --time=1:00:00
##SBATCH -o slurm-%j.out-%N
##SBATCH -e slurm-%j.err-%N


destination_file="corrected_total.csv"
temp_file=".temp_corrected_total.txt"
source_file_1="mapped_percentage.txt"
source_file_2="count_reads.txt"

> "$temp_file"
#echo "#sample_ID, Ref_name, total_reads, #_of_mapped_reads, mapped_percentage" > $destination_file

declare -A sums

while IFS=',' read -r col1 col2 col3 col4 col5 col6 col7; do
	key="$col1,$col2"

	if [[ -n "${sums[$key]}" ]]; then
		((sums["$key"] += col6))
	else
		sums["$key"]=$col6

	fi


done < "$source_file_1"

for key in "${!sums[@]}"; do
    echo "$key,${sums[$key]}" >> "$temp_file"
done




