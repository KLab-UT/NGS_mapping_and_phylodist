#!/bin/bash

destination_file="depth_and_dist.txt"
source_file_1="sample_dist.txt"
source_file_2="depth_percentage.txt"
#source_file_2=$(sed '1d' depth_percentage.txt)

echo "#sample_ID, Ref_name, number_of_reads, avg_depth, map_percentage, patristic_distance" > $destination_file

while IFS="," read -r col_1a col_2a col_3a; do
    while IFS="," read -r col_1b col_2b col_3b col_4b col_5b; do
        if [[ "$col_1a" == "$col_1b" ]]; then
            if [[ "$col_2a" == "$col_2b" ]]; then
                echo "$col_1b,$col_2b,$col_3b,$col_4b,$col_5b,$col_3a" >> $destination_file
            fi
        fi

    done < depth_percentage.txt
done < sample_dist.txt

