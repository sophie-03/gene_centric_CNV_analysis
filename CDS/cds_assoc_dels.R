#library
library(dplyr)
library(tidyr)
library(speedglm)

# Set the working directory
setwd("/data4/smatthews/pheWAS/CDS/")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
group_id <- args[2]
phenotype <- args[3]

# Print the boundary_group_id
print(paste("Processing Group:", group_id))

# Read in data
covariates <- read.table("/data4/smatthews/pheWAS/covariates.txt", header = TRUE)
cases <- read.table(paste("/data4/smatthews/pheWAS/pheno_", phenotype, "/", phenotype, "_cases.txt", sep=""))
cn <- read.table(input_file, header = FALSE)
colnames(cn) <- c("UKB_id","gene","chr","has_cds_deletion")

covariates$cases <- ifelse(covariates$IID %in% cases$V1, 1, 0)

# Merge 'covariates' and 'cn' dataframes based on the UKB_id column
merged_df <- merge(covariates, cn, by.x = "IID", by.y = "UKB_id", all.x = TRUE)

# Replace NA values in the 'has_cds_deletion' column with 0
merged_df$has_cds_deletion <- ifelse(is.na(merged_df$has_cds_deletion), 0, merged_df$has_cds_deletion)

# Perform the logistic regression
logistic_model <- speedglm(cases ~ has_cds_deletion + age + sex + PC1 + PC2 + PC3 + PC4 + batch + smoke, data = merged_df, family = binomial())

# Get the summary of the logistic regression model
model_summary <- summary(logistic_model)

# extract results and put in df
cds_results <- data.frame(
  gene = unique(cn$gene),  # Add gene
  model_summary$coefficients["has_cds_deletion", ],
  stringsAsFactors = FALSE
)
cds_results <- cds_results[, c("gene", "Estimate", "Std.Error", "z.value", "Pr(>|z|)")]

# Generate output filename for each group
output_file <- paste("group_", group_id, "/temp_logistic.txt", sep = "")

# Write output to file (append to the same file for each group)
write.table(cds_results, output_file, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = FALSE, append = TRUE)
