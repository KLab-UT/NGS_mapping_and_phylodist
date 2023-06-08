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

# Check tree
plot(squa_pruned)
edgelabels(squa_pruned$edge.length, bg="black", col="white", font=2)

# See tree class (should be 'phylo')
class(squa_pruned)
# Get distance from ref tips to whip tips using adephylo (install using install.packages('adephylo') if not installed)
library(adephylo)
## Two approaches: get pairwise distances between all tips
pairwise_dist <- distTips(squa_pruned, tips='all', 'patristic')
## Then you will need to extract the distances from each reference to each whiptail species
## OR you can prune the tree to get the comparison you want (and do this for every comparison), example below:
hc_marm <- drop.tip(squa_pruned,squa_pruned$tip.label[-match(c("Aspidoscelis_marmorata","Hemicordylus_capensis"), squa_pruned$tip.label)])
hc_marm_dist <- distTips(hc_marm, tips='all', 'patristic')

