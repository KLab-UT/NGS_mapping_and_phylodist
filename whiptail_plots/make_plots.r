library('ggplot2')
library('gridExtra')

# Set your working directory
# Use this if not in this directory: setwd(<path-to-dir-with-tre-file>)

# Import data as dataframe
dat <- read.csv('../Map_percentage_and_phylodist.csv')
# look at first few rows of data
head(dat)
#dat$log_dist = log(dat$dist)

# Phylogenetic distance ~ Mapping success
mapping_dist <- ggplot(dat, aes(Patristic_distance, Mapped_percentage)) + geom_point(size = 3) + geom_smooth() + theme_classic()
log_mapping_dist <- ggplot(dat, aes(log_patristic_distance, log(Mapped_percentage))) + geom_point(size = 3) + geom_smooth() + theme_classic()
depth_dist <- ggplot(dat, aes(patristic_distance, avg_depth)) + geom_point(size = 3) + geom_smooth(span = 0.5)
reads_dist <- ggplot(dat, aes(patristic_distance, number_of_reads)) + geom_point(size = 3) 

grid.arrange(mapping_dist, depth_dist, nrow=2)
