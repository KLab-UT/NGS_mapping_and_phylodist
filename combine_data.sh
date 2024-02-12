#!/bin/sh
##SBATCH --account=utu
##SBATCH --partition=lonepeak
##SBATCH --nodes=1
##SBATCH --ntasks=1
##SBATCH --time=1:00:00
##SBATCH -o slurm-%j.out-%N
##SBATCH -e slurm-%j.err-%N


destination_file="corrected_total.csv"
temp_file_1=".temp_corrected_total_1.txt"
temp_file_2=".temp_corrected_total_2.txt"
temp_file_3=".temp_corrected_total_3.txt"
sorted_temp_file=".temp_sorted_corrected.txt"
source_file_1="mapped_percentage.txt"
source_file_2="count_reads.txt"

> "$temp_file_1"
> "$temp_file_2"
> "$sorted_temp_file"
> "$destination_file"
echo "#sample_ID, Ref_name, total_reads,#_of_mapped_reads, mapped_percentage" > "$destination_file"



while IFS=',' read -r col1 col2 col3 col4 col5 col6 col7; do
	echo "$col1,$col2,$col6" >> "$temp_file_1"
done < "$source_file_1"

declare -A sums

line_number=0

while IFS=',' read -r col1 col2 col3; do
    ((line_number++))

    # Ignore the first and last lines
    if [ "$line_number" -gt 1 ] && [ "$line_number" -lt "$(wc -l < "$temp_file_1")" ]; then
        key="$col1,$col2"

        if [[ -n "${sums[$key]}" ]]; then
            ((sums["$key"] += col3))
        else
            sums["$key"]=$col3
        fi
    fi


done < "$temp_file_1"

for key in "${!sums[@]}"; do
    echo "$key,,${sums[$key]}," >> "$temp_file_2"
done

sort -t',' -k2,2 "$temp_file_2" > "$sorted_temp_file"

awk -F',' 'NR==FNR{a[$1]=$2; next} $1 in a {$3=a[$1]}1' OFS=',' "$source_file_2" "$sorted_temp_file" > "$temp_file_3"


awk -F',' '{ $5 = $4 / $3; print }' OFS=',' "$temp_file_3" >> "$destination_file"





