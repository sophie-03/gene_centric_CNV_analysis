#this script takes oncogenes from the cosmic gene database file and turns nto a bedfile list of genes

##on cloud

setwd("~/UKB/genes_in_cnvs/Oncogenes")

library(tidyr)
library(dplyr)
library(data.table)

#read in list of cancer genes from cosmic gene database
cancer_genes <- read.table("Cosmic_Gene_Census_23-3-2022.tsv", sep = '\t', header = TRUE, fill = TRUE)
cancer_genes <- as_tibble(cancer_genes)

#filter list of cancer genes for oncogenes only
oncogenes <- cancer_genes[cancer_genes$Role_in_Cancer %like% "oncogene",]

#export list of oncogenes
write.table(oncogenes, "oncogenes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#read in oncogene list
oncogenes <- read.table("oncogenes.txt", header=TRUE, sep = "\t")
oncogenes <- as_tibble(oncogenes)

#seperate location into chr, start and end
oncogenes <- separate(oncogenes, Genome_Location, into =c("chr","geneStart","geneEnd"), sep="[:,-]")

#move columns to be in bed format
oncogenes <- oncogenes %>% select(chr, geneStart, geneEnd, Gene_Symbol, Entrez_Gene_Id)

#order list by chr
oncogenes <- arrange(oncogenes,chr)

#export bed file
write.table(oncogenes, "oncogenes.bed", quote=FALSE, col.names=FALSE, row.names=FALSE, sep ="\t")
