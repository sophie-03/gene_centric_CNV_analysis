setwd("~/../../data2/smatthews/UKB/genes_in_CNVs/TSGs")

#this script takes a list of tumour supressor genes (obtained from cosmic gene census) and turns them into bedfiles

library(tidyr)
library(dplyr)

#read in list of cancer genes from cosmic gene database
cancer_genes <- read.table("Cosmic_Gene_Census_23-3-2022.tsv", sep = '\t', header = TRUE, fill = TRUE)
cancer_genes <- as_tibble(cancer_genes)

#filter list of cancer genes for TSGs only
tsgs <- cancer_genes[cancer_genes$Role_in_Cancer %like% "TSG",]

#export list of tsgs
write.table(TSGs, "TSGs.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#read in TSG list
tsgs <- read.table("TSGs.txt", header=TRUE, sep = "\t")
tsgs <- as_tibble(tsgs)
tsgs

#seperate location into chr, start and end
tsgs <- separate(tsgs, Genome_Location, into =c("chr","geneStart","geneEnd"), sep="[:,-]")
tsgs

#move columns to be in bed format
tsgs <- tsgs %>% select(chr, geneStart, geneEnd, Gene_Symbol, Entrez_Gene_Id)
tsgs

#order list by chr
tsgs <- arrange(tsgs,chr)

#export bed file
write.table(tsgs, "TSGs.bed", quote=FALSE, col.names=FALSE, row.names=FALSE, sep ="\t")
