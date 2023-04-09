#!/bin/bash

{
usage="$(basename "$0") [-h] [-I <SRA_list_Illumina>] [-N <SRA_list_Nanopore>] [-g <reference_genome>] [-d <working_directory]
This program downloads the reference sequences that are being used to test read mapping to divergent genomes
    -h  show this help text
    -d  path to reference directory (where references will be saved)
    -l  path to list of ftp reference paths"
options=':h:d:l::'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    d) d=$OPTARG;;
    l) l=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

echo ""
echo "Working directory:    " $d
echo "path to ref list:     " $l
echo ""

# mandatory arguments
if [ ! "$d" ] || [ ! "$l" ]; then
  echo "arguments -d and -l  must be provided"
  echo "$usage" >&2; exit 1
fi

# download reference genomes
function download_ref {
	cd ${d}
	while read ref; do
		id=$( echo "$ref" | cut -d '/' -f 7 )
		genome_file="${id}_genome.fna.gz"
		annotation_file="${id}_genome.gff.gz"
		rsync --copy-links --times --verbose rsync://"$ref"/genome_file .
		rsync --copy-links --times --verbose rsync://"$ref"/annotation_file .
	done<${l}


