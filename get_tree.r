# We are using the consensus squamate tree from Tonini et al. (2016) - doi: 10.1016/j.biocon.2016.03.039
# We saved the tree (from the supplementary materials) in this repository) ('IBC06772-mmc6.tre')

# Set your working directory
# Use this if not in this directory: setwd(<path-to-dir-with-tre-file>)

# Import ape library
## Install if you haven't using install.packages('ape')
library(ape)

# Read tree (which is a Newick file)
squa <- read.tree('IBC06772-mmc6.tre')
# Verify the 9755 species
squa

# Create vector with species of interest
species_of_interest <- c("Aspidoscelis_gularis","Aspidoscelis_scalaris","Aspidoscelis_inornata","Aspidoscelis_marmorata","Podarcis_muralis","Podarcis_raffoneae","Lacerta_agilis","Hemicordylus_capensis","Sceloporus_undulatus","Python_bivittatus")

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
am_gul <- as.matrix(dist_subset(pairwise_dist, c("Aspidoscelis_gularis","Aspidoscelis_marmorata")))[1,2]

