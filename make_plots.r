library('ggplot2')
library('gridExtra')

# Set your working directory
# Use this if not in this directory: setwd(<path-to-dir-with-tre-file>)

# Import data as dataframe
dat <- read.csv('depth_and_dist.csv')
# look at first few rows of data
head(dat)
#dat$log_dist = log(dat$dist)

# Phylogenetic distance ~ Mapping success
mapping_dist <- ggplot(dat, aes(patristic_distance, mapping_percentage)) + geom_point(aes(colour = cat), size = 3) + scale_color_manual(values = c("total_mapping_success" = "#2C7C9F", "unmerged_mapping_success" = "#2C7C9F", "merged_mapping_success" = "#7CBEE4"))
mapping_dist <- ggplot(dat, aes(patristic_distance, map_percentage)) + geom_point(size = 3) 

grid.arrange(merged_mapping_dist, unmerged_mapping_dist, total_mapping_dist, nrow=3)
