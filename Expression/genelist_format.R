## this script formats the list of human genes into a bed format

setwd("~/GTEX")

library(data.table)
library(tidyr)
library(dplyr)

genelist <- fread("gencode.v41.annotation.gtf", header = FALSE)
genelist <- as_tibble(genelist)

#filter for genes only
genelist <- genelist[genelist$V3 ==  'gene',]
#separate ensembl geneID into a different column
genelist <- separate(genelist, V9, into = c(NA,"gene_id", NA),"\"")
#reorder columns to bed format
genelist <- genelist[,c(1,4,5,3,2,6,7,8,9)]
#keep only columns that are needed
genelist <- genelist %>% select(V1, V4, V5, gene_id, V7)

#export genelist
write.table(genelist, "hg38_genelist.bed", quote=FALSE, col.names=FALSE, row.names=FALSE, sep = "\t")

#add 1kb buffer onto the gene start and end positions
buffer <- genelist
buffer$V4 <- buffer$V4 - 1000
buffer$V5 <- buffer$V5 + 1000

write.table(buffer, "hg38_genelist_buffer.bed", quote=FALSE, col.names=FALSE, row.names=FALSE, sep ="\t")
