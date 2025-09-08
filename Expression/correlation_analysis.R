library(tidyr)
library(plyr)
library(data.table)
library(dplyr)
library(tidyverse)

### FORMAT EXPRESSION DATA
#read in expression file
expression <- fread("~/Documents/PhD/Chap1/GTEX/Data/gene_tpm_whole_blood.gct", sep = "\t")
#remove numbers after . in gene names
#expression <- separate(expression, Name, into = c("Name"), sep = ("\\..*"))
#find the number of duplicated genes and delete duplications
#expression <- expression[!duplicated(expression$Name), ]
#make gene names rownames
expression <- expression %>% column_to_rownames(var="Name")

#delete irrelevant columns
expression <- expression %>% dplyr::select(-c(id,Description))

#keep only ID number for individual IDs
n <- names(expression)
a <- data.frame(n)
a <- separate(a, n, into = c(NA, "individual"), sep="-")
colnames(expression) <- a$individual
rm(a)

#transpose so genes are columns
h <- t(expression)
expression <- as.data.frame(h)
rm(h)

#order rows and columns
expression <- expression[ order(row.names(expression)), ]
expression <- expression[, order(names(expression))]


#label genes
gene_labels <- read.table("~/Documents/PhD/Chap1/GTEX/Data/gene_list_biomart.txt", sep = "\t", header=TRUE)


## FORMAT COPY NUMBER DATA
#read in gene/copy number data
bed_output <- read.table("~/Documents/PhD/Chap1/GTEX/Data/all_samples_buffer_intersect.out")
bed_output <- as_tibble(bed_output)
#give column names
bed_output <- bed_output %>% rename(c( "chr" = V1,
                                       "gene_start" = V2,
                                       "gene_end" = V3,
                                       "gene" = V4,
                                       "direction" = V5,
                                       "individual" = V6,
                                       "cnv_start" = V8,
                                       "cnv_end" = V9,
                                       "cnv_type" = V10,
                                       "read_depth" = V11))
#just keep individual ID numbers 
bed_output <- separate(bed_output, individual, into = c( "individual"), sep = "\\..*")  #remove file extension from ID 

#make cn df so genes = col and rows = ID
column_names <- unique(bed_output$individual) #get list of individuals 
row_names <- unique(bed_output$gene) #get list of genes
cn <- data.frame(matrix(ncol=length(column_names),nrow=length(row_names))) # make empty df with correct length
colnames(cn) <- column_names # name columns after individuals
rownames(cn) <- row_names # name rows as genes
cn[cbind(match(bed_output$gene, rownames(cn)),
         match(bed_output$individual, colnames(cn)))] <- bed_output$read_depth # populate new data frame with read depth (copy number) from no_XY

#fill in individuals with no copy number as 1 (the reference)
cn <- replace(cn, is.na(cn),"1")

#format individual IDs so we only have the ID number
n <- names(cn)
a <- data.frame(n)
a <- separate(a, n, into = c(NA, "individual"), sep="-")
colnames(cn) <- a$individual
#remove duplicated columns
dups <- duplicated(colnames(cn))
cn <- cn[!dups]

#format gene IDs to remove the version number
#remove numbers after . on genes
#cn <- tibble::rownames_to_column(cn, "gene") #first have to convert row names to a column to edit
#cn <- separate(cn, gene, into = c("gene"), sep = ("\\..*"))
#find the number of duplicated genes and delete duplications
#cn <- cn[!duplicated(cn$gene),]
#convert column back to row names
#rownames(cn) <- cn$gene
#cn$gene <- NULL

#transpose so genes are columns
g <- t(cn)
cn <- as.data.frame(g)
rm(g)

#order rows and columns
cn <- cn[ order(row.names(cn)), ]
cn <- cn[, order(names(cn))]


## FILTER
#filter expression for genes that have median TPM > 0.5
meds <- apply(expression, 2, median)
expressed_genes <- expression[, meds >= 0.5]

#filter expression and cn dataframes to only contain individuals and genes that are found in both data
#filter individuals
indv_list <- intersect(rownames(cn),rownames(expressed_genes))
cn <- cn[(row.names(cn) %in% indv_list),]
ex <- expressed_genes[(row.names(expressed_genes) %in% indv_list),]
#remove version number of genes
#colnames(cn) <- sub("\\..*", "", colnames(cn))
#colnames(ex) <- sub("\\..*", "", colnames(ex))
#filter genes
overlap_genes <- intersect(names(cn),names(ex))
cn <- cn[,(colnames(cn) %in% overlap_genes)]
ex <- ex[,(colnames(ex) %in% overlap_genes)]

#order both dataframes so they are in the same order
ex <- ex[ order(row.names(ex)), ]
cn <- cn[ order(row.names(cn)), ]
ex <- ex[, order(colnames(ex)) ]
cn <- cn[, order(colnames(cn)) ]


## CORRELATIONS - DELETIONS
#empty dataframe
cor_results_del <- data.frame(matrix(ncol=3))
cor_results_del <- cor_results_del[0,]
#remove duplications from the dataframe so they don't influence the correlations, we're only looking at deletions here
cn_del <- cn
cn_del[cn_del > 1] <- NA
##correlation test - if RD is less than 1 in > 5% of individuals
for(i in 1:ncol(cn_del)) {
  x = as.numeric(cn_del[,i])
  
  # Check if more than 5% of the values are < 1
  if(length(x[x < 1]) > 0.05*length(x)) {
    
    # Use tryCatch to handle potential errors in correlation test
    tryCatch({
      cor = cor.test(as.numeric(cn_del[,i]), ex[,i], method="spearman", use = "complete.obs")
      df  <- data.frame(colnames(cn_del)[i], cor$estimate, cor$p.value)
      cor_results_del <- rbind(cor_results_del, df)
    }, error = function(e) {
      # If an error occurs, print a message and continue to the next iteration
      message(paste("Error in correlation test for column", colnames(cn_del)[i], ":", e$message))
    })
  }
}


# Remove rows where cor.estimate is NA
cor_results_del <- cor_results_del[!is.na(cor_results_del$cor.estimate), ]

colnames(cor_results_del) <- c("Gene_stable_ID_version","cor.estimate","p_value")
cor_results_del$Gene_stable_ID <- sub("\\..*$", "", cor_results_del$Gene_stable_ID_version)

cor_results_del <- left_join(cor_results_del, gene_labels, by = "Gene_stable_ID")

# correction for multiple tests
cor_results_del$p.adjusted <- p.adjust(cor_results_del$p_value, method = "BH")
sig_results_del <- cor_results_del[cor_results_del$p_value < 0.05, ]

#histogram
hist(cor_results_del$cor.estimate, main = "Distribution of correlation coefficients")
hist(sig_results_del$cor.estimate, main = "Distribution of correlation coefficients (p.adjust > 0.05")

# how many of the genes are positively correlated
pos_cor_del <- filter(sig_results_del, cor.estimate > 0)
#negatively correlated
neg_cor_del <- filter(sig_results_del, cor.estimate < 0)

## IMMUNE GENES
## are immune genes enriched in the significant results?

#read in immune genes from innateDB 28/07/25
immune_genes <- read.table("~/Documents/PhD/Chap1/GTEX/Data/InnateDB_genes.txt", head = TRUE, sep = "\t", fill = TRUE)

#get immune genes in correlation genes
intersect(immune_genes$ensembl, cor_results_del$Gene_stable_ID)
#subset for immune genes
immune_subset <- cor_results_del %>%
  filter(Gene_stable_ID %in% immune_genes$ensembl)

no.correlation_genes <- nrow(cor_results_del)
no.significant_genes <- nrow(sig_results_del)
no.correlation_immune_genes <- nrow(immune_subset)
no.correlation_sig_immune_genes <- nrow(immune_subset[immune_subset$p.adjusted < 0.05, ])

#enrichment test
# Calculate non-immune counts
no.correlation_nonimmune_genes <- no.correlation_genes - no.correlation_immune_genes
no.significant_nonimmune_genes <- no.significant_genes - no.correlation_significant_genes

# Calculate non-significant counts for immune and non-immune
immune_not_significant <- no.correlation_immune_genes - no.correlation_significant_genes
nonimmune_not_significant <- no.correlation_nonimmune_genes - no.significant_nonimmune_genes

# Build contingency table
contingency_table <- matrix(
  c(
    no.correlation_significant_genes, no.significant_nonimmune_genes,
    immune_not_significant, nonimmune_not_significant
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(contingency_table) <- c("Significant", "Not Significant")
colnames(contingency_table) <- c("Immune", "Non-immune")

print(contingency_table)

fisher.test(contingency_table)


## CORRELATIONS - DUPLICATIONS
empty dataframe
cor_results_dup <- data.frame(matrix(ncol=3))
cor_results_dup <- cor_results_dup[0,]
#remove duplications from the dataframe so they don't influence the correlations, we're only looking at deletions here
cn_dup <- cn
cn_dup[cn_dup < 1] <- NA
for(i in 1:ncol(cn_dup)) {
  x = as.numeric(cn_dup[,i])
  
  # Check if more than 5% of the values are > 1
  if(length(x[x > 1]) > 0.05*length(x)) {
    
    # Use tryCatch to handle potential errors in correlation test
    tryCatch({
      cor = cor.test(as.numeric(cn_dup[,i]), ex[,i], method="spearman", use = "complete.obs")
      df  <- data.frame(colnames(cn_dup)[i], cor$estimate, cor$p.value)
      cor_results_dup <- rbind(cor_results_dup, df)
    }, error = function(e) {
      # If an error occurs, print a message and continue to the next iteration
      message(paste("Error in correlation test for column", colnames(cn_dup)[i], ":", e$message))
    })
  }
}


# Remove rows where cor.estimate is NA
cor_results_dup <- cor_results_dup[!is.na(cor_results_dup$cor.estimate), ]

colnames(cor_results_dup) <- c("Gene_stable_ID_version","cor.estimate","p_value")
cor_results_dup$Gene_stable_ID <- sub("\\..*$", "", cor_results_dup$Gene_stable_ID_version)

cor_results_dup <- left_join(cor_results_dup, gene_labels, by = "Gene_stable_ID")

# correction for multiple tests
cor_results_dup$p.adjusted <- p.adjust(cor_results_dup$p_value, method = "BH")
sig_results_dup <- cor_results_dup[cor_results_dup$p_value < 0.05, ]

#histogram
hist(cor_results_dup$cor.estimate, main = "Distribution of correlation coefficients")
hist(sig_results_dup$cor.estimate, main = "Distribution of correlation coefficients (p.adjust > 0.05")

# how many of the genes are positively correlated
pos_cor_dup <- filter(sig_results_dup, cor.estimate > 0)
#negatively correlated
neg_cor_dup <- filter(sig_results_dup, cor.estimate < 0)

## IMMUNE GENES
## are immune genes enriched in the significant results?

#read in immune genes from innateDB 28/07/25
immune_genes <- read.table("~/Documents/PhD/Chap1/GTEX/Data/InnateDB_genes.txt", head = TRUE, sep = "\t", fill = TRUE)

#get immune genes in correlation genes
intersect(immune_genes$ensembl, cor_results_dup$Gene_stable_ID)
#subset for immune genes
immune_subset <- cor_results_dup %>%
  filter(Gene_stable_ID %in% immune_genes$ensembl)

no.correlation_genes <- nrow(cor_results_dup)
no.significant_genes <- nrow(sig_results_dup)
no.correlation_immune_genes <- nrow(immune_subset)
no.correlation_sig_immune_genes <- nrow(immune_subset[immune_subset$p.adjusted < 0.05, ])

#enrichment test
# Calculate non-immune counts
no.correlation_nonimmune_genes <- no.correlation_genes - no.correlation_immune_genes
no.significant_nonimmune_genes <- no.significant_genes - no.correlation_significant_genes

# Calculate non-significant counts for immune and non-immune
immune_not_significant <- no.correlation_immune_genes - no.correlation_significant_genes
nonimmune_not_significant <- no.correlation_nonimmune_genes - no.significant_nonimmune_genes

# Build contingency table
contingency_table <- matrix(
  c(
    no.correlation_significant_genes, no.significant_nonimmune_genes,
    immune_not_significant, nonimmune_not_significant
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(contingency_table) <- c("Significant", "Not Significant")
colnames(contingency_table) <- c("Immune", "Non-immune")

print(contingency_table)

fisher.test(contingency_table)

