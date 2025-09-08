setwd("~/../../data2/smatthews/UKB/genes_in_CNVs")

#this script expands the gene locations by 1kb each side to try and include promoters

library(tidyr)

#import data
genes <- read.table("hg19_ensembl_genes.bed")
genes <- as_tibble(genes)

#add 1000bp to the end position and subtract 1000bp from start position
buffer <- genes
buffer$V2 <- buffer$V2 - 1000
buffer$V3 <- buffer$V3 + 1000

#keep only genes in chrs
buffer <- buffer[buffer$V1 %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                      "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                      "chr17","chr18","chr19","chr20","chr21","chr22"),]

write.table(buffer, "hg19_ensembl_genes_buffer.bed", quote=FALSE, col.names=FALSE, row.names=FALSE, sep ="\t")
