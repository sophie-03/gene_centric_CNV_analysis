##this script creates a file that maps sample number from log and BAF files to UKB ID

setwd("~/../../data/UKB/sophie")

library(dplyr)

#read in fam file
fam <- read.table("fam_files/ukb22431_c22_b0_v2_s488146.fam", header=FALSE)

#create new data frame containing ID and sample number
map <- data.frame(fam$V1)

#add corresponding log/baf sample number to df
map$sample <- seq.int(nrow(map))

#rename columns
map <- rename(map, UKB_id = fam.V1)

#export
write.table(map, "map_file.txt", sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = TRUE)
