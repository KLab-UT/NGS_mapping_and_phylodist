# We are using the consensus squamate tree from Tonini et al. (2016) - doi: 10.1016/j.biocon.2016.03.039
# We saved the tree (from the supplementary materials) in this repository) ('IBC06772-mmc6.tre')

# Set your working directory
# Use this if not in this directory: setwd(<path-to-dir-with-tre-file>)

# Import ape library
## Install if you haven't using install.packages('ape'=Analysis of Phylogenetics and Evolution)
library(ape)

# Read tree (which is a Newick file)
squa <- read.tree('IBC06772-mmc6.tre')
# Verify the 9755 species
squa

# Create vector with species of interest
#species_of_interest <- c("Aspidoscelis_gularis","Aspidoscelis_scalaris","Aspidoscelis_inornata","Aspidoscelis_marmorata","Podarcis_muralis","Podarcis_raffoneae","Lacerta_agilis","Hemicordylus_capensis","Sceloporus_undulatus","Python_bivittatus","Calyptommatus_sinebrachiatus","Euleptes_europaea","Gekko_japonicus","Salvator_merianae","Sphenodon_punctatus","Tretioscincus_oriximinensis")


# Prune squa
squa_pruned <- drop.tip(squa,squa$tip.label[-match(species_of_interest, squa$tip.label)])
squa_pruned
# Set plot size
par(mar = c(5, 5, 4, 0) + 0.1)  # Adjust the margin values as needed
# Make branches longer so text doesn't cover species name
#plot(squa_pruned, edge.width = 2, xlim = c(100, 1100), ylim = c(-10, 110))
plot(ladderize(squa_pruned, right = TRUE), edge.width = 2)  # ladderize() makes the branches longer
# The right = TRUE parameter is used to specify that the branches should be lengthened on the right side of the tree
#laddersize() isn't working as expected. The main issue is the text boxes change size as I increase the plot size, 
#so when I make the plot bigger the problem (the text box covering the species name) gets worse.

# Check tree
plot(squa_pruned, edge.width = 1)
squa_pruned$edge.length <- round(squa_pruned$edge.length, digits = 2)
edgelabels(squa_pruned$edge.length, bg="black", col="white", font=2, cex = 0.6, box.lwd = 0, xpad = 5, ypad = 5)

# See tree class (should be 'phylo')
class(squa_pruned)
# Get distance from ref tips to whip tips using adephylo (install using install.packages('adephylo') if not installed)
library(adephylo)
## Two approaches: get pairwise distances between all tips
# Extract data for each pair of interest
#library('usedist')

####################################################################
ref_genomes <- c("Aspidoscelis_gularis","Aspidoscelis_scalaris","Aspidoscelis_inornata","Aspidoscelis_marmorata","Podarcis_muralis","Podarcis_raffoneae","Lacerta_agilis","Hemicordylus_capensis","Sceloporus_undulatus","Python_bivittatus","Calyptommatus_sinebrachiatus","Euleptes_europaea","Gekko_japonicus","Salvator_merianae","Sphenodon_punctatus","Tretioscincus_oriximinensis")
#bigger list of 27 genomes
#("Aspidoscelis_marmoratus", "Podarcis_muralis", "Podarcis_raffonei", "Lacerta_agilis", "Hemicordylus_capensis", "Sceloporus_undulatus", "Python_bivittatus", "Gekko_japonicus", "Euleptes_europaea", "Alligator_mississippiensis", "Falco_peregrinus", "Bison_bison", "Homo_sapiens", "Mus_musculus", "Rana_temporaria", "Danio_rerio", "Fundulus_heteroclitus", "Protopterus_annectens", "Stegostoma_fasciatum", "Pristis_pectinata", "Petromyzon_marinus", "Salvator_merianae", "Calyptommatus_sinebrachiatus", "Tretioscincus_oriximinensis", "phenodon_punctatus", "Ambystoma_mexicanum", "Eptatretus_burgeri")
samples <- c("Aspidoscelis_gularis","Aspidoscelis_scalaris","Aspidoscelis_inornata","Aspidoscelis_marmorata")
pairwise_dist <- distTips(squa_pruned, tips='all', 'patristic')
pd_matrix <- as.matrix(pairwise_dist)
pd_matrix

a <- pd_matrix["Sphenodon_punctatus","Aspidoscelis_gularis"]
a
pairwise_dist_df <- as.data.frame(pd_matrix)
write.csv(pairwise_dist_df, "test.csv", row.names = TRUE)


# make genome list
sample_dist <- c()

file_path <- "sample_dist.txt"
file_connection <- file(file_path, "w")
for (g in ref_genomes){
  asp_sep <- "KLC098"
  asp_gul <- "RLK004"
  asp_inor <- "RLK019"
  asp_marm <- "RLK034"
  dist <- pd_matrix["Aspidoscelis_scalaris", g]
  dist_line <- paste(asp_sep, g, dist, sep = " ")
  #sample_dist <- append(sample_dist, dist_line)
  writeLines(dist_line, con = file_connection)
  
  dist <- pd_matrix["Aspidoscelis_gularis", g]
  dist_line <- paste(asp_gul, g, dist, sep = " ")
  #sample_dist <- append(sample_dist, dist_line)
  writeLines(dist_line, con = file_connection)
  
  dist <- pd_matrix["Aspidoscelis_inornata", g]
  dist_line <- paste(asp_inor, g, dist, sep = " ")
  sample_dist <- append(sample_dist, dist_line)
  writeLines(dist_line, con = file_connection)
  
  dist <- pd_matrix["Aspidoscelis_marmorata", g]
  dist_line <- paste(asp_marm, g, dist, sep = " ")
  sample_dist <- append(sample_dist, dist_line)
  writeLines(dist_line, con = file_connection)
}
close(file_connection)
sample_dist
# file_path <- "sample_dist.txt"
# file_connection <- file(file_path, "w")
# for (dline in sample_dist) {
#   writeLines(dline, con = file_connection)
# }
#close(file_connection)

################################################################
# am_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Aspidoscelis_marmorata")))[1,2]
# am_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris","Aspidoscelis_marmorata")))[1,2]
# am_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata","Aspidoscelis_marmorata")))[1,2]
# am_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata","Aspidoscelis_marmorata")))[1,2]
# 
# # Hemicordylus capensis ref
# hc_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Hemicordylus_capensis")))[1,2]
# hc_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris","Hemicordylus_capensis")))[1,2]
# hc_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata","Hemicordylus_capensis")))[1,2]
# hc_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata","Hemicordylus_capensis")))[1,2]
# 
# # Lacerta agilis ref
# la_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Lacerta_agilis")))[1,2]
# la_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris","Lacerta_agilis")))[1,2]
# la_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata","Lacerta_agilis")))[1,2]
# la_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata","Lacerta_agilis")))[1,2]
# 
# # Python bivittatus ref
# pb_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Python_bivittatus")))[1,2]
# pb_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris","Python_bivittatus")))[1,2]
# pb_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata","Python_bivittatus")))[1,2]
# pb_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata","Python_bivittatus")))[1,2]
# 
# # Podarcis muralis ref
# pm_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Podarcis_muralis")))[1,2]
# pm_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris","Podarcis_muralis")))[1,2]
# pm_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata","Podarcis_muralis")))[1,2]
# pm_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata","Podarcis_muralis")))[1,2]
# 
# # Podarcis raffoneae ref
# pr_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Podarcis_raffoneae")))[1,2]
# pr_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris","Podarcis_raffoneae")))[1,2]
# pr_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata","Podarcis_raffoneae")))[1,2]
# pr_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata","Podarcis_raffoneae")))[1,2]
# 
# # Sceloporus undulatus ref
# su_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Sceloporus_undulatus")))[1,2]
# su_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris","Sceloporus_undulatus")))[1,2]
# su_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata","Sceloporus_undulatus")))[1,2]
# su_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata","Sceloporus_undulatus")))[1,2]
# 
# ## Randy just added this one, the other new ones need to be addedhenodon punctatushenodon punctatus
# # Sphenodon punctatus ref
# sp_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Sphenodon punctatus")))[1,2]
# sp_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris","Sphenodon punctatus")))[1,2]
# sp_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata","Sphenodon punctatus")))[1,2]
# sp_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata","Sphenodon punctatus")))[1,2]

