library(CNVRanger)
library(GenomicRanges)
library(data.table)
library(dplyr)

# ---- Get command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)

# Phenotype name
pheno_name <- args[1]
part_num <- as.integer(args[2])

#read in cnv calls
cnv_file <- paste0("/data4/smatthews/pheWAS/cnv_GWAS/cnv_regions_splits/cnv_regions_part", part_num, ".txt")
cnvs_with_regions <- fread(cnv_file)

# get data for GWAS
pheno_file <- paste0("/data4/smatthews/pheWAS/pheno_", pheno_name, "/", pheno_name, "_cases.txt")
pheno <- read.table(pheno_file)
covar <- read.table("/data4/smatthews/pheWAS/covariates.txt", header = TRUE)

# add pheno data to covar df
covar$pheno <- ifelse(covar$IID %in% pheno$V1, 1, 0)

# Get all unique region_ids
all_regions <- unique(cnvs_with_regions$region_id)

# Initialize results dataframe
results_df <- data.frame(
  region_id = character(),
  beta = numeric(),
  se = numeric(),
  z_value = numeric(),
  p_value = numeric(),
  n_carriers = numeric(),
  n_total = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each region
for(region_no in all_regions) {
  tryCatch({
  
  # Subset data for current region
  region <- subset(cnvs_with_regions, region_id == region_no)
  
  # Merge with covariates
  data <- left_join(covar, region, by = c("IID" = "UKB_id"))
  
  # Set individuals with no CNV to state = 2 (normal)
  data$state <- ifelse(is.na(data$state), 2, data$state)
  
  # Remove individuals with CN > 2 (duplications) - focusing on deletions
  data <- data[data$state <= 2, ]
  
  # Perform logistic regression
  logistic_model <- glm(pheno ~ state + age + sex + PC1 + PC2 + PC3 + PC4 + batch + smoke, 
                        data = data, family = binomial)
  
  # Extract statistics for the 'state' coefficient
  model_summary <- summary(logistic_model)
  stats <- model_summary$coefficients["state",]
    
    # Add results to dataframe
    results_df <- rbind(results_df, data.frame(
      region_id = region_no,
      beta = stats["Estimate"],
      se = stats["Std. Error"],
      z_value = stats["z value"],
      p_value = stats["Pr(>|z|)"],
      n_carriers = sum(data$state < 2),  # Number of deletion carriers
      n_total = nrow(data),
      stringsAsFactors = FALSE
    ))

    print(region_no)
  }, error = function(e) {
    # This block runs only if there's an error
    warning(paste("Error in region", region_no, ":", e$message))
    # The loop will continue to the next region
  })
}


## WRITE OUTPUT
# ---- Safe append to master output ----
output_file <- paste0("/data4/smatthews/pheWAS/cnv_GWAS/", pheno_name, "_del_logistic_results.txt")
lock_file <- paste0(output_file, ".lock")

# Write results to a temporary file first
tmp_file <- tempfile(pattern = paste0("part", part_num, "_"), fileext = ".txt")
fwrite(results_df, tmp_file, sep = "\t", col.names = (part_num == 1))  # header only for first

# Append safely using a lock
lock <- lock(file = lock_file, exclusive = TRUE, timeout = 6000)
cat("Appending results from part", part_num, "to master file...\n")

if (!file.exists(output_file)) {
  # if master file doesn’t exist, include header
  fwrite(results_df, output_file, sep = "\t", col.names = TRUE, append = FALSE)
} else {
  # append without header
  fwrite(results_df, output_file, sep = "\t", col.names = FALSE, append = TRUE)
}

unlock(lock)
