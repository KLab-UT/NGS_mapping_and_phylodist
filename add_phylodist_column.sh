#!/bin/bash

destination_file="Map_percentage_and_phylodist.csv"
source_file_1="sample_dist.txt"
source_file_2="total_reads_avg_depth.csv"
#source_file_2=$(sed '1d' depth_percentage.txt)

echo "Sample_ID, Ref_name, Total_reads, Mapped_reads, Mapped_proportion, Avg_depth, Mapped_percentage, Patristic_distance" > $destination_file

while IFS="," read -r col_1a col_2a col_3a; do
    while IFS="," read -r col_1b col_2b col_3b col_4b col_5b col_6b col_7b; do
        if [[ "$col_1a" == "$col_1b" ]]; then
            if [[ "$col_2a" == "$col_2b" ]]; then
                echo "$col_1b,$col_2b,$col_3b,$col_4b,$col_5b,$col_6b,$col_7b,$col_3a" >> $destination_file
            fi
        fi

    done < ${source_file_2}
done < ${source_file_1}

