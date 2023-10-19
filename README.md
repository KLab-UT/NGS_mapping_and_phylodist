# whiptail_mapping_and_phylogenetic_distance
Scripts used to examine the difference in depth percentage between 28 vertibrate genomes and sexual species of whiptail lizards (genus Aspidoscelis) using whole-genome sequencing data.

---

# Contents

- [Background](#background)
- [Starting Files](#starting-files)
- [General Pipeline](#general-pipeline)

---

# <a name="background"></a>
# Background


Parthenogenetic whiptail lizards have been shown to have lower endurance capacity 
([Cullum, 1997](https://www.journals.uchicago.edu/doi/abs/10.1086/286055?casa_token=q8DOEvxRkccAAAAA:nmq4l99bzJ7XY8vxokdkj0eRg6816F4_zQ9VSFx7sstxB_qBfty9GAPVe1uUGPgpuMU7CZL4ySIZ); [Klabacka et al., 2022](https://www.journals.uchicago.edu/doi/full/10.1086/719014?casa_token=_E1ccM7e3WkAAAAA%3A1JC_ft2sxeGGwmoiBGjjWNGuLMJX-gXpmfMZsWbjXGbXV4iFVCKvK1R8vbg92gTPLfhYSnbAPYNt))
 and mitochondrial respiration ([Klabacka et al., 2022](https://www.journals.uchicago.edu/doi/full/10.1086/719014?casa_token=_E1ccM7e3WkAAAAA%3A1JC_ft2sxeGGwmoiBGjjWNGuLMJX-gXpmfMZsWbjXGbXV4iFVCKvK1R8vbg92gTPLfhYSnbAPYNt)) compared to closely related sexual species. 

 We hypothesize that lower endurance and mitochondrial respiration in these lizards is a result of reduced genetic compatibility in Nuc-Mt\* genes due to either 
 (A) inter-genomic interactions or (B) intra-genomic interactions. Inter- and intra- prefixes are in reference to the parental genomes of the parthenogenetic whiptails, which are of hybrid origin (i.e., the crossing of two divergent species resulted in the new evolutionary lineage of the parthenogenetic species). Inter-genomic refers to interactions between parental genomes (between gene products of maternal ancestry and gene products of paternal ancestry) and intra-genomic refers to interactions within a parental genome (e.g., between nuclear gene products of maternal ancestry and mitochondrial gene products of maternal ancestry). 

Examining the Nuc-Mt genetic variability between the sexual parental species provides a base for testing the inter-genomic hypothesis: If we find variation in Nuc-Mt genes between hybridizing species, this variation is likely present in the 
 "frozen"\*\* 
 genomes 
of the parthenogenetic species. If no variation is present in Nuc-Mt genes between hybridizing species, then a source other than reduced compatibility between the divergent genomes is responsible for the reduced performance in parthenogens (e.g., intra-genomic interactions). 
As a first step in this project, we will simply be examining which reference genome is the best to map to. We will make quantitative comparisons to 6 annottated reference genomes: The common wall lizard (*Podarcis muralis*), the aeolian wall lizard (*Podarcis raffonei*), the sand lizard (*Lacerta agilis*), the false girdled lizard (*Hemicordylus capensis*), the eastern fence lizard (*Sceloporus undulatus*), and the Burmese python (*Python bivittatus*).

---

# <a name="starting-files"></a>
# Starting Files
We performed paired-end Illumina sequencing for four parental sexual species that are closely related and are all involved in hybridization that has resulted in parthenogenetic lineages: *Aspidoscelis inornatus*, *Aspidoscelis gularis*, *A. marmoratus*, and *A. septemvittatus* (a.k.a. *A. gularis septemvittatus*). Details for the ```.fastq``` files for these species are below:

|                Filename                  | SampleID |      Species           | Size            |
|:----------------------------------------:|:--------:|:----------------------:|:---------------:|
| KLC098_USD16091388L_HKFJFDSXX_L4_1.fq.gz |  KLC098  | A. septemvittatus      | 7.2 GB          |
| KLC098_USD16091388L_HKFJFDSXX_L4_2.fq.gz |  KLC098  | A. septemvittatus      | 7.4 GB          |
| KLC098_USD16091388L_HKG5MDSXX_L1_1.fq.gz |  KLC098  | A. septemvittatus      | 3.1 GB          |
| KLC098_USD16091388L_HKG5MDSXX_L1_2.fq.gz |  KLC098  | A. septemvittatus      | 3.2 GB          |
| RLK004_USD16091389L_HJNHCDSXX_L3_1.fq.gz |  RLK004  | A. gularis             | 11 GB           |
| RLK004_USD16091389L_HJNHCDSXX_L3_2.fq.gz |  RLK004  | A. gularis             | 12 GB           |
| RLK019_USD16091390L_HJMVKDSXX_L3_1.fq.gz |  RLK019  | A. inornatus           | 8.2 GB          |
| RLK019_USD16091390L_HJMVKDSXX_L3_2.fq.gz |  RLK019  | A. inornatus           | 8.6 GB          |
| RLK019_USD16091390L_HKG5MDSXX_L1_1.fq.gz |  RLK019  | A. inornatus           | 1.7 GB          |
| RLK019_USD16091390L_HKG5MDSXX_L1_2.fq.gz |  RLK019  | A. inornatus           | 1.8 GB          |
| RLK034_USD16091387L_HJNHCDSXX_L2_1.fq.gz |  RLK034  | A. marmoratus          | 11 GB           |
| RLK034_USD16091387L_HJNHCDSXX_L2_2.fq.gz |  RLK034  | A. marmoratus          | 12 GB           |

You will have to obtain these files for your own working environment. The raw *Aspidoscelis* reads are available on the CHPC lonepeak cluster here:

```
/scratch/general/nfs1/utu_4310/whiptail_shared_data/raw_data
```

You should copy all of the files into a scratch directory of your own that uses your uNID. For example, if my uNID was "u0123456" I would do the following:

```
cd /scratch/general/nfs1/
mkdir -p u0123456/whiptail_nmt_variation_wd/raw_data
```

You can then copy the files to your directory. You will need to do this either as a batch submission or using an interactive job. You can copy files in parallel (i.e., using multiple threads) using the shell tool [GNU parallel](https://www.gnu.org/software/parallel/). In order to do this, you need to make sure you specify the number of cores (executing units within a processor) you plan to use. Your cores will be able to perform tasks simultaneously, this is known as "in parallel". You will generally specify the number of "tasks" rather than the number of "cores", which is essentially the same concept (specifying the number of tasks is telling the computer "this is how many things I want to be done in parallel"). If my username was u0123456, I could do the following (either in my batch script or interactive job):

```
# Setup interactive job with 10 hr time limit and 12 tasks (cores)
salloc --time=10:00:00 -N 1 -n 12
# Move into the directory that currently has the genomic data that needs to be copied
cd /scratch/general/nfs1/utu_4310/whiptail_shared_data/raw_data
# Copy the genomic files in parallel to my working directory
ls *.fq.gz | time parallel -j12 --eta --bar cp {} /scratch/general/nfs1/u0612750/whiptail_nmt_variation_wd/raw_data/
```
Because there is no annotated reference genome for Aspidoscelis, we will map the raw data to the annotated genomes of species with varying levels of relatedness. You will copy the annotated reference files from Genbank to your working directory using the (NCBI FTP)[https://ftp.ncbi.nlm.nih.gov/genomes/] with the ```rsync``` command. If my username was u0123456, I could do the following (either in my batch script or interactive job):

```
# Move into the working directory and create a directory for the Reference
cd /scratch/general/nfs1/u0612750/whiptail_nmt_variation_wd
mkdir -p references
cd references
```

Then you can copy the reference genome you plan to use. For this project, we will be examining read mapping success with three reference genomes:

1. Podarcis muralis
```
### --- fasta
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Podarcis_muralis/latest_assembly_versions/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_genomic.fna.gz .
### --- annotation
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Podarcis_muralis/latest_assembly_versions/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_genomic.gff.gz .
```
2. Podarcis raffonei
```
### --- fasta
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Podarcis_raffonei/latest_assembly_versions/GCF_027172205.1_rPodRaf1.pri/GCF_027172205.1_rPodRaf1.pri_genomic.fna.gz .
### --- annotation
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Podarcis_raffonei/latest_assembly_versions/GCF_027172205.1_rPodRaf1.pri/GCF_027172205.1_rPodRaf1.pri_genomic.gff.gz .
```
3. Lacerta agilis
```
### --- fasta
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Lacerta_agilis/latest_assembly_versions/GCF_009819535.1_rLacAgi1.pri/GCF_009819535.1_rLacAgi1.pri_genomic.fna.gz
### --- annotation
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Lacerta_agilis/latest_assembly_versions/GCF_009819535.1_rLacAgi1.pri/GCF_009819535.1_rLacAgi1.pri_genomic.gff.gz .
```
4. Hemicordylus capensis
```
### --- fasta
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Hemicordylus_capensis/latest_assembly_versions/GCF_027244095.1_rHemCap1.1.pri/GCF_027244095.1_rHemCap1.1.pri_genomic.fna.gz .
### --- annotation
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Hemicordylus_capensis/latest_assembly_versions/GCF_027244095.1_rHemCap1.1.pri/GCF_027244095.1_rHemCap1.1.pri_genomic.gff.gz .
```
5. Sceloporus undulatus
```
### --- fasta
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Sceloporus_undulatus/latest_assembly_versions/GCF_019175285.1_SceUnd_v1.1/GCF_019175285.1_SceUnd_v1.1_genomic.fna.gz .
### --- annotation
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Sceloporus_undulatus/latest_assembly_versions/GCF_019175285.1_SceUnd_v1.1/GCF_019175285.1_SceUnd_v1.1_genomic.gff.gz .
```
6. Python bivittatus
```
### --- fasta
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Python_bivittatus/latest_assembly_versions/GCF_000186305.1_Python_molurus_bivittatus-5.0.2/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.fna.gz .
### --- annotation
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Python_bivittatus/latest_assembly_versions/GCF_000186305.1_Python_molurus_bivittatus-5.0.2/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.gff.gz .
```

## BIOL 4310 Project
As a class, we will create a pipeline that contains quality check, read cleaning, mapping, and analysis steps. Responsibility for the coding for these steps was randomly assigned and are as follows:

- Quality checking steps: Seun
- Read cleaning step: Candice
- Read mapping step: Dante
- Analysis: Everyone

You will each contribute to the pipeline contained within the bash file ```fastq_2_bam.sh```. Your sections are marked by comments- you should only add code to your respective section (this way we can avoid merge conflicts). We will work on sections that require everyone's input in class. IF you would like to move forward, you should make a copy of the ```fastq_2_bam.sh``` file that you can use to create your own pipeline. I encourage this! It will allow you to work on all of the pieces together. However, for the official ```fastq_2_bam.sh``` file you should only add code to your assigned section.




