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

file_path <- "../Map_percentage_and_phylodist.csv"
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

#data_sorted <- data[order(data, $map_percentages), ]


file_path <- "../full_corrected_data.csv"
data <- read.csv(file_path) 
#data_sorted <- data[order(data$Mapped_percentage), ]
map_percentages <- data[, 7]
species_names <- data[, 2]
#geom_bar(map_reads)
pdf("map_percentage_bar.pdf")
plot <- par(mar = c(2, 4, 4, 0))
par(cex.axis = 0.3)
barplot(map_percentages,
        names.arg = species_names,
        col = "blue",
        main = "Comparison of Mapped Percentage",
        xlab = "",
        ylab = "Mapped Percentage",
		las = 2)  # Rotates x-axis labels vertically for readability
	
print(plot)
dev.off()


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
