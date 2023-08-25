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
species_of_interest <- c("Aspidoscelis_gularis","Aspidoscelis_scalaris","Aspidoscelis_inornata","Aspidoscelis_marmorata","Podarcis_muralis","Podarcis_raffoneae","Lacerta_agilis","Hemicordylus_capensis","Sceloporus_undulatus","Python_bivittatus","Calyptommatus_sinebrachiatus","Euleptes_europaea","Gekko_japonicus","Salvator_merianae","Sphenodon_punctatus","Tretioscincus_oriximinensis")



# Prune squa
squa_pruned <- drop.tip(squa,squa$tip.label[-match(species_of_interest, squa$tip.label)])

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
pairwise_dist <- distTips(squa_pruned, tips='all', 'patristic')
# Extract data for each pair of interest
library('usedist')

# make genome list

ref_genomes <- c("Aspidoscelis_gularis","Aspidoscelis_scalaris","Aspidoscelis_inornata","Aspidoscelis_marmorata","Podarcis_muralis","Podarcis_raffoneae","Lacerta_agilis","Hemicordylus_capensis","Sceloporus_undulatus","Python_bivittatus","Calyptommatus_sinebrachiatus","Euleptes_europaea","Gekko_japonicus","Salvator_merianae","Sphenodon_punctatus","Tretioscincus_oriximinensis")
#bigger list of 27 genomes
#("Aspidoscelis_marmoratus", "Podarcis_muralis", "Podarcis_raffonei", "Lacerta_agilis", "Hemicordylus_capensis", "Sceloporus_undulatus", "Python_bivittatus", "Gekko_japonicus", "Euleptes_europaea", "Alligator_mississippiensis", "Falco_peregrinus", "Bison_bison", "Homo_sapiens", "Mus_musculus", "Rana_temporaria", "Danio_rerio", "Fundulus_heteroclitus", "Protopterus_annectens", "Stegostoma_fasciatum", "Pristis_pectinata", "Petromyzon_marinus", "Salvator_merianae", "Calyptommatus_sinebrachiatus", "Tretioscincus_oriximinensis", "phenodon_punctatus", "Ambystoma_mexicanum", "Eptatretus_burgeri")
#for (g in ref_genomes){
#print(dist_subset(pairwise_dist, c("Aspidoscelis_gularis", g)))


############################ 
asp_gul <- as.matrix(pairwise_dist, c("Aspidoscelis_gularis", ref_genomes))
asp_gul_df <- data.frame(Reference = ref_genomes, Distance = asp_gul)
print(asp_gul_df)
write.csv(asp_gul_df, "phylodist.csv", row.names = TRUE)
############################

asp_gul_list <- numeric(length(ref_genomes))
asp_sep_list <- c()
asp_inor_list <- c()
asp_marm_list <- c()
for (i in seq_along(ref_genomes)){
  asp_gul <- (dist_subset(pairwise_dist, c("Aspidoscelis_gularis", ref_genomes)))
  asp_gul_list[i] <- asp_gul
  asp_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris", ref_genomes)))[1,2]
  asp_sep_list[i] <- asp_sep
  asp_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata", ref_genomes)))[1,2]
  asp_inor_list[i] <- asp_inor
  asp_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata", ref_genomes)))[1,2]
  asp_marm_list[i] <- asp_marm
}
print(paste(asp_gul_list))
# for (ref_genome in ref_genomes){
#   asp_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis", ref_genome)))[1,2]
#   asp_gul_list <- c(asp_gul_list, asp_gul)
#   asp_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris", ref_genome)))[1,2]
#   asp_sep_list <- c(asp_sep_list, asp_sep)
#   asp_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata", ref_genome)))[1,2]
#   asp_inor_list <- c(asp_inor_list, asp_inor)
#   asp_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata", ref_genome)))[1,2]
#   asp_marm_list <- c(asp_marm_list, asp_marm)
# }
result_df <- data.frame(
  "Genome" = genome_list,
  "Aspidoscelis_gularis" = asp_gul_list
)
print(result_df)
  
genome_result <- data.frame(
  "Genome" = ref_genomes,
  "AG" = asp_gul_list,
  "AS" = asp_sep_list,
  "AI" = asp_inor_list,
  "AM" = asp_marm_list
)
print(genome_result)
  
#result_df <- rbind(result_df, genome_result)


# Call the function and store the results
#result_table <- Get_Dist(ref_genomes)
#Get_Dist(ref_genomes)



# #iterate through genome list to get phylogenetic distance
# Get_Dist <- function(genome_list){
# 	asp_gul_list <- list()
# 	asp_sep_list <- list()
# 	asp_inor_list <- list()
# 	asp_marm_list <- list()
# 	for (genome in genome_list){
# 	  
# 		asp_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis",genome)))[1,2]
# 		asp_gul_list <- c(asp_gul_list, asp_gul)
# 		asp_sep <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_scalaris",genome)))[1,2]
# 		asp_sep_list <- c(asp_sep_list, asp_sep)
# 		asp_inor <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_inornata",genome)))[1,2]
# 		asp_inor_list <- c(asp_inor_list, asp_inor)
# 		asp_marm <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_marmorata",genome)))[1,2]
# 		asp_marm_list <- append(asp_marm_list, asp_marm)
# 	}
#   row_names <- c("Aspidoscelis_gularis", "Aspidoscelis_scalaris", "Aspidoscelis_inornata", "Aspidoscelis_marmorata")
#   #make matrix then table of phylogenetic distances
#   phylodist_matrix <- matrix("", nrow = length(row_names), ncol = length(ref_genomes), dimnames = list(row_names, ref_genomes))
#   phylodist_table <- as.data.frame(phylodist_matrix, stringsAsFactors = FALSE)
# 	cat("ref_genomes: ", paste(ref_genomes, collapse = ", "), "\n")
# 	print("Phylogenetic Distances:")
# 	cat("AG:", paste(asp_gul_list, collapse = ", "), "\n")
# 	cat("AS:", paste(asp_sep_list, collapse = ", "), "\n")
# 	cat("AI:", paste(asp_inor_list, collapse = ", "), "\n")
# 	cat("AM:", paste(asp_marm_list, collapse = ", "), "\n")
# }
# Get_dist(ref_genomes)

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

