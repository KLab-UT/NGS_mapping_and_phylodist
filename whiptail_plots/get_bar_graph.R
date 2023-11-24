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
mapped_percentages <- data[, 5]
species_names <- data[, 1]
#geom_bar(depth_percentage)
par(mar = c(14, 4, 4, 0))
barplot(mapped_percentages,
        names.arg = species_names,
        col = "green",
        main = "Comparison of Species to Whiptails",
        xlab = "",
        ylab = "Mapped Percentage",
        las = 2)  # Rotates x-axis labels vertically for readability
