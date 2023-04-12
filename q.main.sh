#!/bin/sh
#SBATCH --account=utu_4310
#SBATCH --partition=lonepeak
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N


# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data
mkdir -p $wd/cleaned_reads
mkdir -p $wd/cleaned_reads/merged_reads
mkdir -p $wd/cleaned_reads/unmerged_reads
mkdir -p $wd/references
mkdir -p $wd/mapped_reads/pm
mkdir -p $wd/mapped_reads/pr
mkdir -p $wd/mapped_reads/hc
mkdir -p $wd/mapped_reads/su
mkdir -p $wd/mapped_reads/la
mkdir -p $wd/mapped_reads/pb

# download reference genomes
bash get_refs.sh -d $wd/references -l ref_list.txt

# check read quality of raw reads
bash check_qc.sh -i $wd/raw_reads -t 10

# clean reads
bash clean_reads.sh -i $wd/raw_reads -o $wd/cleaned_reads -t 10

# check read quality of cleaned reads
bash check_qc.sh -i $wd/cleaned_reads -t 10

# map cleaned reads to reference
## Candice
bash map_reads.sh -i $wd/cleaned_reads -r $wd/references/Podarcis_muralis.fna -o $wd/mapped_reads/pm/ -t 10
bash map_reads.sh -i $wd/cleaned_reads -r $wd/references/Podarcis_raffonei.fna -o $wd/mapped_reads/pr/ -t 10

## Dante
bash map_reads.sh -i $wd/cleaned_reads -g $wd/references/Hemicordylus_capensis.fna -o $wd/mapped_reads/hc/ -t 10
bash map_reads.sh -i $wd/cleaned_reads -g $wd/references/Sceloporus_undulatus.fna -o $wd/mapped_reads/su/ -t 10

## Seun
bash map_reads.sh -i $wd/cleaned_reads -r $wd/references/Python_bivittatus.fna -o $wd/mapped_reads/pb/ -t 10
bash map_reads.sh -i $wd/cleaned_reads -r $wd/references/Lacerta_agilis.fna -o $wd/mapped_reads/la/ -t 10

# examine mapping results
bash mapping_analysis.sh -i $wd/mapped_reads -l $wd/mapped_reads/ref_species_list.txt -t 10
