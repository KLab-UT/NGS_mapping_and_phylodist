# This script was used to add the average depth data to the corrected_total.csv (which Laisha created using combine_data.sh)

correct_tr <- read.csv('corrected_total.csv')
correct_ad <- read.csv('Merged_mapped_percentage.txt')

correct_ad$unique_id <- paste0(correct_ad$Sample_ID, correct_ad$X.Ref_name)
correct_tr$unique_id <- paste0(correct_tr$X.sample_ID, correct_tr$Ref_name)

total_reads_avg_depth <- merge(correct_tr, correct_ad, "unique_id", all=T)

total_reads_avg_depth = subset(total_reads_avg_depth, select = c(X.sample_ID, Ref_name, total_reads, X._of_mapped_reads, mapped_percentage, Avg_depth))

names(total_reads_avg_depth)[names(total_reads_avg_depth) == "X.sample_ID"] <- "sample_ID"
names(total_reads_avg_depth)[names(total_reads_avg_depth) == "mapped_percentage"] <- "mapped_proportion"
names(total_reads_avg_depth)[names(total_reads_avg_depth) == "X._of_mapped_reads"] <- "n_mapped_reads"
names(total_reads_avg_depth)[names(total_reads_avg_depth) == "mapped_percentage"] <- "mapped_proportion"
names(total_reads_avg_depth)[names(total_reads_avg_depth) == "Avg_depth"] <- "avg_depth"

total_reads_avg_depth$mapped_percentage <- total_reads_avg_depth$mapped_proportion * 100

write.csv(total_reads_avg_depth, file="total_reads_avg_depth.csv",quote=F,row.names = F)



