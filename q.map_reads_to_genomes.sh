#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --nodes=1
#SBATCH --ntasks=27
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N

# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data
mkdir -p $wd/cleaned_reads
mkdir -p $wd/cleaned_reads/merged_reads
mkdir -p $wd/cleaned_reads/unmerged_reads
mkdir -p $wd/references
mkdir -p $wd/mapped_reads/Aspidoscelis_marmoratus
mkdir -p $wd/mapped_reads/Alligator_mississippiensis
mkdir -p $wd/mapped_reads/Ambystoma_mexicanum
mkdir -p $wd/mapped_reads/Bison_bison
mkdir -p $wd/mapped_reads/Calyptommatus_sinebrachiatus
mkdir -p $wd/mapped_reads/Danio_rerio
mkdir -p $wd/mapped_reads/Eptatretus_burgeri
mkdir -p $wd/mapped_reads/Euleptes_europaea
mkdir -p $wd/mapped_reads/Fundulus_heteroclitus
mkdir -p $wd/mapped_reads/Falco_peregrinus
mkdir -p $wd/mapped_reads/Gekko_japonicus
mkdir -p $wd/mapped_reads/Hemicordylus_capensis
mkdir -p $wd/mapped_reads/Homo_sapiens
mkdir -p $wd/mapped_reads/Lacerta_agilis
mkdir -p $wd/mapped_reads/Mus_musculus
mkdir -p $wd/mapped_reads/Podarcis_muralis
mkdir -p $wd/mapped_reads/Podarcis_raffonei
mkdir -p $wd/mapped_reads/Python_bivittatus
mkdir -p $wd/mapped_reads/Protopterus_annectens
mkdir -p $wd/mapped_reads/Pristis_pectinata
mkdir -p $wd/mapped_reads/Petromyzon_marinus
mkdir -p $wd/mapped_reads/Rana_temporaria
mkdir -p $wd/mapped_reads/Stegostoma_fasciatum
mkdir -p $wd/mapped_reads/Sceloporus_undulatus
mkdir -p $wd/mapped_reads/Salvator_merianae
mkdir -p $wd/mapped_reads/Sphenodon_punctatus
mkdir -p $wd/mapped_reads/Tretioscincus_oriximinensis

# map cleaned reads to reference:
MapReads() {
	bash map_reads.sh -i $wd/cleaned_reads/merged_reads -g $wd/references/${2} -o $wd/mapped_reads/${1} -t 27
	bash map_reads.sh -i $wd/cleaned_reads/unmerged_reads -g $wd/references/${2} -o $wd/mapped_reads/${1} -t 27
}

# use references in ref_genome in the MapReads function and run each line in parallel
# {1} is the directory name and {2} is the file name of the genome reference mathcing that directory
grep -v '^#' ref_genomes.txt | parallel --colsep ' ' MapReads {1} {2}

