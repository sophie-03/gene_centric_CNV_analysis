library(speedglm)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
group_id <- args[2]
phenotype_file <- args[3]
covariate_file <- args[4]
type <- args[5]


# Read in data
covariates <- read.table(covariate_file, header = TRUE)
cases <- read.table(phenotype_file)
cn <- read.table(input_file, header = FALSE)
colnames(cn) <- c("chr","geneStart", "geneEnd", "EnsemblID", "strand", "gene",
                     "chro", "cnvStart", "cnvEnd", "cn", "IID")

# print the gene
print(paste("Processing gene:", unique(cn$gene))

# add phenotype column to covariates df
covariates$cases <- ifelse(covariates$IID %in% cases$V1, 1, 0)

# Merge covariates and cn dataframes based on the 'sample_id' column
merged_df <- merge(covariates, cn[, c("IID", "cn")], by = "IID", all.x = TRUE)

# Replace NA values in the 'cn' column with 2 (normal copy number)
merged_df$cn <- ifelse(is.na(merged_df$cn), 2, merged_df$cn)

## filter for either deletions or duplications
if (type == "deletion") {
    merged_df <- merged_df[merged_df$cn <= 2, ]
} else if (type == "duplication") {
    merged_df <- merged_df[merged_df$cn >= 2, ]
}


# get covariates
covariate_names <- setdiff(names(merged_df), c("cases", "cn", "IID"))
# create formula
formula <- paste("cases ~ cn +", paste(covariate_names, collapse = " + "))
# Perform the logistic regression
logistic_model <- speedglm(formula, data = merged_df, family = binomial())

# Get the summary of the logistic regression model
model_summary <- summary(logistic_model)

cds_results <- data.frame(
  gene = unique(cn$gene),  # Add gene
  model_summary$coefficients["cn", ],
  stringsAsFactors = FALSE
)

# Generate output filename for each group
output_file <- paste("group_", group_id, "/temp_logistic.txt", sep = "")

# Write output to file (append to the same file for each group)
write.table(cn_results, output_file, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = FALSE, append = TRUE)



