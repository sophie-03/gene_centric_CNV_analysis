setwd("~/UKB/genes_in_cnvs")

library(tidyr)
library(dplyr)

#read in genes
genes <- read.table("biomart_grch37_genes.txt", header = F, skip = 1)

#order genes by chr
genes <- genes[order(genes$V1),]
#remove scaffolds
genes <- filter(genes, V1 %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                         11, 12, 13, 14, 15, 16, 17, 18,
                         19, 20, 21, 22, "X", "Y"))
#add chr to chr column
genes$V1 <- paste("chr", genes$V1, sep="")

#export as bed
write.table(genes, "biomart_grch37_genes.bed", quote=FALSE, col.names=FALSE, row.names=FALSE, sep ="\t")
