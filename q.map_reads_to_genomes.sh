#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --time=72:00:00
#SBATCH --nodes=3
#SBATCH --ntasks=36
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N

# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data
mkdir -p $wd/trimmed_reads
mkdir -p $wd/trimmed_reads/merged_reads
mkdir -p $wd/trimmed_reads/unmerged_reads
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
#mkdir -p $wd/mapped_reads/Protopterus_annectens
mkdir -p $wd/mapped_reads/Pristis_pectinata
mkdir -p $wd/mapped_reads/Petromyzon_marinus
mkdir -p $wd/mapped_reads/Rana_temporaria
mkdir -p $wd/mapped_reads/Stegostoma_fasciatum
mkdir -p $wd/mapped_reads/Sceloporus_undulatus
mkdir -p $wd/mapped_reads/Salvator_merianae
mkdir -p $wd/mapped_reads/Sphenodon_punctatus
mkdir -p $wd/mapped_reads/Tretioscincus_oriximinensis
mkdir -p $wd/mapped_reads/Latimeria_chalumnae

# map trimmed reads to reference:
echo "Beggining mapping"
MapReads() {
	wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data
	echo "############################"
	echo -e "Species: ${1}\nGenome: ${2}"
	# you're passing in n threads where n is number of reads (6) multiplied by number of threads used by functions in map_reads.sh (2).
	# This is done for every genome you want to map (27 genomes listed in ref_genomes.txt).
	#bash map_reads.sh -i $wd/trimmed_reads/merged_reads -g $wd/references/${2} -o $wd/mapped_reads/${1} -t 12
	echo "merched read mapped"
	bash map_reads_unmerged.sh -i $wd/trimmed_reads/unmerged_reads -g $wd/references/${2} -o $wd/mapped_reads/${1} -t 12
	echo "unmerged read mapped"
}
export -f MapReads

# use references in ref_genome in the MapReads function and run each line in parallel
# {1} is the directory name and {2} is the file name of the genome reference mathcing that directory
# grep -v filter out lines starting with #

#grep -v '^#' ref_genomes.txt | cut -d " " -f 1,2 | parallel MapReads {1} {2}
grep -v '^#' ref_genomes.txt | parallel --colsep ' ' MapReads {1} {2}
echo "mapping done"

# Keep in mind that if your MapReads function itself uses multiple threads (e.g., with the -t option set to 4 as in your example), the total number of threads used will be a combination of those used by MapReads and those used by parallel.

