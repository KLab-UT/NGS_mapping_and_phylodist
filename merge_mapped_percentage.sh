#!/bin/bash

destination_file="Merged_mapped_percentage.txt"
temp_file=".temp_Merged_mapped_percentage.txt"
source_file_1="mapped_percentage.txt"

> $temp_file

while IFS="," read -r col_1a col_2a col_3a col_4a col_5a col_6a col_7a; do
    while IFS="," read -r col_1b col_2b col_3b col_4b col_5b col_6b col_7b; do
        if [[ "$col_1a" == "$col_1b" ]]; then
            if [[ "$col_2a" == "$col_2b" ]]; then
                if [[ "$col_3a" != "$col_3b" ]]; then
                    total_reads=$((${col_4a} + ${col_4b}))
                    mapped_reads=$((${col_6a} + ${col_6b}))
                    mapped_percentage=$(echo "scale=9; $mapped_reads / $total_reads" | bc)
                    echo "$col_2b,$col_1b,$total_reads,$mapped_reads,$mapped_percentage" >> $temp_file
                fi
            fi
        fi
    done < ${source_file_1}
done < ${source_file_1}

echo "#Ref_name, sample_ID, Total_reads, Mapped_reads, Mapped_percentage" > $destination_file
sort -u $temp_file >> $destination_file

