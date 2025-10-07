# Set the working directory
setwd(paste("/define/home/smatthews/pheWAS/binary/pheno_", commandArgs(trailingOnly = TRUE)[3], sep=""))

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
boundary_group_id <- args[2]
phenotype <- args[3]

# Print the boundary_group_id
print(paste("Processing Boundary Group:", boundary_group_id))

# Read in data
covariates <- read.table("/define/home/smatthews/pheWAS/covariates.txt", header = TRUE)
cases <- read.table(paste(phenotype, "_cases.txt", sep=""))
cn <- read.table(input_file, header = FALSE)

covariates$cases <- ifelse(covariates$IID %in% cases$V1, 1, 0)

# Extract the value after '=' in the 4th column of 'cn'
cn$cn_value <- as.numeric(sub(".*cn=(\\d+)", "\\1", cn$V4))

# Merge 'covariates' and 'cn' dataframes based on the 'V5' column
merged_df <- merge(covariates, cn[, c("V5", "cn_value")], by.x = "IID", by.y = "V5", all.x = TRUE)

# Replace NA values in the 'cn' column with 2
merged_df$cn_value <- ifelse(is.na(merged_df$cn_value), 2, merged_df$cn_value)

# remove individuals with a cn less than 2
merged_df <- merged_df[merged_df$cn_value >= 2, ]

# Perform the logistic regression
logistic_model <- glm(cases ~ cn_value + age + sex + PC1 + PC2 + PC3 + PC4 + batch + smoke, data = merged_df, family = binomial)

# Get the summary of the logistic regression model
model_summary <- summary(logistic_model)

# Extract the results for the variable 'cn_value' and put it into a dataframe
cn_value_results <- as.data.frame(t(model_summary$coefficients["cn_value", ]))

# Add gene location to the dataframe
cn_value_results <- cbind(cn_value_results, V1 = unique(cn$V1))

# Generate output filename for each group
output_file <- paste("group_", boundary_group_id, "/temp_logistic.txt", sep = "")

# Write output to file (append to the same file for each group)
write.table(cn_value_results, output_file, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = FALSE, append = TRUE)
