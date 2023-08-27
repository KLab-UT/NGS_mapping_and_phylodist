#!/bin/sh
#SBATCH --account=utu
#SBATCH --partition=lonepeak
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N


# create working environment
wd=/scratch/general/nfs1/utu_4310/whiptail_shared_data
mkdir -p $wd/trimmed_reads
mkdir -p $wd/trimmed_reads/merged_reads
mkdir -p $wd/trimmed_reads/unmerged_reads
mkdir -p $wd/references
mkdir -p $wd/mapped_reads/Aspidoscelis_marmoratus
mkdir -p $wd/mapped_reads/Podarcis_muralis
mkdir -p $wd/mapped_reads/Podarcis_raffonei
mkdir -p $wd/mapped_reads/Hemicordylus_capensis
mkdir -p $wd/mapped_reads/Sceloporus_undulatus
mkdir -p $wd/mapped_reads/Lacerta_agilis
mkdir -p $wd/mapped_reads/Python_bivittatus

# download reference genomes: Randy
bash get_refs.sh -d $wd/references -l ref_list.txt

# check read quality of raw reads: Seun
bash check_qc.sh -i $wd/raw_reads -t 10

# clean reads: Candice
bash clean_reads.sh -i $wd/raw_reads -o $wd/trimmed_reads -t 10

# check read quality of trimmed reads: Seun
bash check_qc.sh -i $wd/trimmed_reads -t 10

# map trimmed reads to reference: Dante
bash map_reads.sh -i $wd/trimmed_reads/merged_reads -g $wd/references/GCA_014337955.1_AspMar1.0_genomic.fna.gz -o $wd/mapped_reads/Aspidoscelis_marmoratus -t 24
bash map_reads.sh -i $wd/trimmed_reads/merged_reads -g $wd/references/GCF_004329235.1_PodMur_1.0_genomic.gff.gz -o $wd/mapped_reads/Podarcis_muralis/ -t 24
bash map_reads.sh -i $wd/trimmed_reads/merged_reads -g $wd/references/GCF_027172205.1_rPodRaf1.pri_genomic.gff.gz -o $wd/mapped_reads/Podarcis_raffonei/ -t 24
bash map_reads.sh -i $wd/trimmed_reads/merged_reads -g $wd/references/GCF_027244095.1_rHemCap1.1.pri_genomic.fna.gz -o $wd/mapped_reads/Hemicordylus_capensis/ -t 24
bash map_reads.sh -i $wd/trimmed_reads/merged_reads -g $wd/references/GCF_019175285.1_SceUnd_v1.1_genomic.fna.gz -o $wd/mapped_reads/Sceloporus_undulatus/ -t 24
bash map_reads.sh -i $wd/trimmed_reads/merged_reads -g $wd/references/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.fna.gz -o $wd/mapped_reads/Python_bivittatus/ -t 24
bash map_reads.sh -i $wd/trimmed_reads/merged_reads -g $wd/references/GCF_009819535.1_rLacAgi1.pri_genomic.fna.gz -o $wd/mapped_reads/Lacerta_agilis/ -t 24

# examine mapping results
bash mapping_analysis.sh -i $wd/mapped_reads -l $wd/mapped_reads/ref_species_list.txt -t 10
