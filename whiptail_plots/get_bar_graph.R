# ctrl + shift + c to comment

library(ggplot2)
# geom_bar(
#   mapping = NULL,
#   data = NULL,
#   stat = "count",
#   position = "stack",
#   ...,
#   just = 0.5,
#   width = NULL,
#   na.rm = FALSE,
#   orientation = NA,
#   show.legend = NA,
#   inherit.aes = TRUE
# )

file_path <- "Map_percentage_and_phylodist_sorted.csv"
data <- read.csv(file_path) 
map_percentages <- data[, 5]
avg_depth <- data[, 4]
species_names <- data[, 1]
#geom_bar(map_percentage)
par(mar = c(14, 4, 4, 0))
par(cex.axis = 1.0)
barplot(avg_depth,
        names.arg = species_names,
        col = "green",
        main = "Comparison of Species to Whiptails",
        xlab = "",
        ylab = "Average Depth",
        las = 2)  # Rotates x-axis labels vertically for readability

file_path <- "top_3_map_percentage.csv"
data <- read.csv(file_path) 
map_reads <- data[, 4]
species_names <- data[, 1]
#geom_bar(map_reads)
par(mar = c(4, 4, 2, 0))
par(cex.axis = 0.7)
barplot(map_reads,
        names.arg = species_names,
        col = "green",
        main = "Comparison of #_of_Mapped_reads",
        xlab = "",
        ylab = "# of Mapped Reads")
        #las = 2)  # Rotates x-axis labels vertically for readability

file_path <- "top_3_map_percentage.csv"
data <- read.csv(file_path) 
avg_depth <- data[, 6]
species_names <- data[, 1]
#geom_bar(avg_depth)
par(mar = c(4, 4, 2, 0))
par(cex.axis = 0.7)
barplot(avg_depth,
        names.arg = species_names,
        col = "green",
        main = "Comparison of Avg_depth",
        xlab = "",
        ylab = "Average Depth")
        #las = 2)  # Rotates x-axis labels vertically for readability
