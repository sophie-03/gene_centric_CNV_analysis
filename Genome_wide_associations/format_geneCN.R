# script to format gene CN data into the correct format for input into ParseCNV2

setwd("/data4/smatthews/pheWAS/gene_cnvs/any_overlap")

library(data.table)
library(dplyr)

#load geneCn data
calls <- fread("protein_coding_genes_cnvs.txt")

#name columns to avoid confusion
colnames(calls) <- c("chr","geneStart", "geneEnd", "EnsemblID", "strand", "gene",
                     "chro", "cnvStart", "cnvEnd", "cn", "sample", "UKB_id")

#filter for genes only (ENSG not ENST - which are transcripts)
calls <- calls[!grepl("ENST", calls$EnsemblID),]

# merge first 3 columns to make column format gene:4-7
calls$col1 <- paste(calls$chr, calls$geneStart, sep = ":")
calls$col1 <- paste(calls$col1, calls$geneEnd, sep = "-")

## add numsnp column
calls$col2 <- "numsnp=1"

#add length column
calls$col3 <- paste("length=", calls$geneEnd - calls$geneStart, sep = "")
#add state/cn column
calls$col4 <-case_when(
  calls$cn == "0" ~ "state1,cn=0",
  calls$cn == "1" ~ "state2,cn=1",
  calls$cn == "2" ~ "state5,cn=3", #this only happens in chrX (changes in accordance with https://parsecnv.sourceforge.net/)
  calls$cn == "3" ~ "state5,cn=3",
  calls$cn == "4" ~ "state6,cn=4"
)
#id columns
calls$col5 <- calls$UKB_id
#start/end snp column
calls$col6 <- paste("startsnp=", calls$chr, sep = "")
calls$col6 <- paste(calls$col6, calls$geneStart, sep = "_")
calls$col7 <- paste("endsnp=", sep = "")
calls$col7 <- paste(calls$col7, calls$chr, sep = "")
calls$col7 <- paste(calls$col7, calls$geneEnd, sep = "_")

#create new dataframe with just the columns we want
formatted <- select(calls, col1, col2, col3, col4, col5, col6, col7)

#order based on gene and id
formatted <- formatted[order(formatted$col1, formatted$col5),]

#export
write.table(formatted, "protein_coding_geneCN.rawcnv", sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = FALSE)

#split calls into files for each chr
for (i in 1:22){
  chr <- formatted[grepl(paste("chr", i ,":",sep = ""), formatted$col1),]

  write.table(chr, paste("protein_coding_chr", i, "_geneCN.rawcnv", sep=""), sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = FALSE)}


chrX <- formatted[grepl("chrX:", formatted$col1),]
write.table(formatted, "protein_coding_chrX_geneCN.rawcnv", sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = FALSE)
