#! /bin/bash

sed 's/\r$//g' depth_and_dist.txt > .temp_remove_M_line_endings.txt 
cp .remove_M_line_endings_intermediate.txt depth_and_dist.txt

# \r: This represents the carriage return character (commonly used for line endings in Windows text files).
# $: This matches the end of a line.
