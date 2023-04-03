# whiptail_nmt_variation
Scripts used to examine the amount of genetic variation in nuclear-encoded mitochondrial-targetting (NucMt) genes between sexual species of whiptail lizards (genus Aspidoscelis) using whole-genome sequencing data

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
<!--- and mitochondrial respiration ([Klabacka et al., 2022](https://www.journals.uchicago.edu/doi/full/10.1086/719014?casa_token=_E1ccM7e3WkAAAAA%3A1JC_ft2sxeGGwmoiBGjjWNGuLMJX-gXpmfMZsWbjXGbXV4iFVCKvK1R8vbg92gTPLfhYSnbAPYNt)). We hypothesize that lower endurance and mitochondrial respiration in these lizards is a result of reduced genetic compatibility in Nuc-Mt\* genes due to either (A) inter-genomic interactions or (B) intra-genomic interactions. Inter- and intra- prefixes are in reference to the parental genomes of the parthenogenetic whiptails, which are of hybrid origin (i.e., the crossing of two divergent species resulted in the new evolutionary lineage of the parthenogenetic species). Inter-genomic refers to interactions between parental genomes (between gene products of maternal ancestry and gene products of paternal ancestry) and intra-genomic refers to interactions within a parental genome (e.g., between nuclear gene products of maternal ancestry and mitochondrial gene products of maternal ancestry). Examining the Nuc-Mt genetic variability between the sexual parental species provides a base for testing the inter-genomic hypothesis: If we find variation in Nuc-Mt genes between hybridizing species, this variation is likely present in the "frozen"\*\* genomes of the parthenogenetic species. If no variation is present in Nuc-Mt genes between hybridizing species, then a source other than reduced compatibility between the divergent genomes is responsible for the reduced performance in parthenogens (e.g., intra-genomic interactions). --->
---

# <a name="starting-files"></a>
# Starting Files
We performed paired-end Illumina sequencing for four parental sexual species that are closely related and are all involved in hybridization that has resulted in parthenogenetic lineages: *Aspidoscelis inornatus*, *Aspidoscelis gularis*, *A. marmoratus*, and *A. septemvittatus* (a.k.a. *A. gularis septemvittatus*). Details for the ```.fastq``` files for these species are below:

|                Filename                  | SampleID |      Species      | Size    |
|:----------------------------------------:|:--------:|:-----------------:|:-------:|
| KLC098_USD16091388L_HKFJFDSXX_L4_1.fq.gz |  KLC098  | A. septemvittatus | 7.2 GB  |
| KLC098_USD16091388L_HKFJFDSXX_L4_2.fq.gz |  KLC098  | A. septemvittatus | 7.4 GB  |
| KLC098_USD16091388L_HKG5MDSXX_L1_1.fq.gz |  KLC098  | A. septemvittatus | 3.1 GB  |
| KLC098_USD16091388L_HKG5MDSXX_L1_2.fq.gz |  KLC098  | A. septemvittatus | 3.2 GB  |
| RLK004_USD16091389L_HJNHCDSXX_L3_1.fq.gz |  RLK004  | A. gularis        | 11 GB   |
| RLK004_USD16091389L_HJNHCDSXX_L3_2.fq.gz |  RLK004  | A. gularis        | 12 GB   |
| RLK019_USD16091390L_HJMVKDSXX_L3_1.fq.gz |  RLK019  | A. inornatus      | 8.2 GB  |
| RLK019_USD16091390L_HJMVKDSXX_L3_2.fq.gz |  RLK019  | A. inornatus      | 8.6 GB  |
| RLK019_USD16091390L_HKG5MDSXX_L1_1.fq.gz |  RLK019  | A. inornatus      | 1.7 GB  |
| RLK019_USD16091390L_HKG5MDSXX_L1_2.fq.gz |  RLK019  | A. inornatus      | 1.8 GB  |
| RLK034_USD16091387L_HJNHCDSXX_L2_1.fq.gz |  RLK034  | A. marmoratus     | 11 GB   |
| RLK034_USD16091387L_HJNHCDSXX_L2_2.fq.gz |  RLK034  | A. marmoratus     | 12 GB   |


