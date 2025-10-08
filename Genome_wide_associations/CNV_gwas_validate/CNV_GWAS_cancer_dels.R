library(CNVRanger)
library(GenomicRanges)
library(data.table)
library(dplyr)

#read in cnv calls
cnvs <- read.table("/data4/smatthews/pheWAS/cnv_GWAS/cnv_calls.bed")
colnames(cnvs) <- c("chr","start","end","state","sample_id","UKB_id")

#group the calls by sample ID and convert them to a GRangesList
grl <- GenomicRanges::makeGRangesListFromDataFrame(cnvs, 
                                            split.field="sample_id", keep.extra.columns=TRUE)
#sort
grl <- GenomicRanges::sort(grl)

#summarize cnv calls to cnv regions
cnvrs <- populationRanges(grl, density=0.1)

#convert regions to df
cnvrs_df <- as.data.frame(cnvrs)
cnvrs_df$region_id <- paste0("region_",1:nrow(cnvrs_df))

# Convert both to data.tables
cnvs_dt <- as.data.table(cnvs)
cnvrs_dt <- as.data.table(cnvrs_df)


# combine cnv calls with cnv regions
cnvs_with_regions <- cnvs_dt[cnvrs_dt, 
                             on = .(chr == seqnames, 
                                    start >= start, 
                                    end <= end),
                             nomatch = NULL,  # Remove non-matches
                             .(chr = seqnames, start = i.start, end = i.end, 
                               state, sample_id, UKB_id, region_id, freq, type)]

# get data for GWAS
pheno <- read.table("/data4/smatthews/pheWAS/pheno_cancer/cancer_cases.txt")
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
  
}

write.table(results_df, "/data4/smatthews/pheWAS/cnv_GWAS/cancer_del_logistic_results.txt", col.names = TRUE, row.names = FALSE, quotes = FALSE)
