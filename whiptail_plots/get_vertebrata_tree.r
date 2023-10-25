# Set your working directory
# Use this if not in this directory: setwd(<path-to-dir-with-tre-file>)

# Import ape library
## Install if you haven't using install.packages('ape'=Analysis of Phylogenetics and Evolution)
library(ape)
library(phytools)

# Read tree (which is a Newick file)
vert_tree <- read.tree('vertebrata.tre')

plotTree(vert_tree, direction="upwards",edge.width = 2, label.offset = 0.5)

