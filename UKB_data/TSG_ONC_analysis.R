## Setup


### association between TSG/Onc copy number and cancer from UKB data - Fisher's and logistic regression

#libraries
library(data.table)
library(dplyr)
#library(ggplot2)
library(statmod)

#setwd
setwd("/data4/smatthews/pheWAS/UKB_tsg_onc")


#list of all individual IDs
all_IDs <- read.table("map_file.txt", header = TRUE)

#read in list of tumour suppressors
tsgs <- read.table("TSGs.bed")
colnames(tsgs) <- c("chr", "start", "end", "gene", "entrez_id")

#read in list of oncogenes
oncogenes <- read.table("oncogenes.bed")
colnames(oncogenes) <- c("chr", "start", "end", "gene", "entrez_id")


#read in cancer data
cancer <- read.table("40009.tab", header=TRUE, fill = TRUE)
#keep only individuals who have had reported incidences of cancer
cancer <- na.omit(cancer)
#list of ids
cancer <- data.frame(cancer$f.eid)
#rename column
colnames(cancer) <- "UKB_id"
#subset to individuals who we have cnvs called for 
cancer <- filter(cancer, UKB_id %in% all_IDs$UKB_id)
#add another column for when merge with other dataframes
cancer$cancer <- "yes"


## DRIVER MUTATIONS
drivers <- read.csv("driver_genes_Bailey2018.csv", sep = ",")
#subset for just list of genes and if tsg or oncogene
drivers <- drivers %>% dplyr::select(gene, Tumor.suppressor.or.oncogene.prediction..by.20.20..)

#subset driver genes for tsgs
driver_tsg <- filter(drivers, Tumor.suppressor.or.oncogene.prediction..by.20.20.. == "tsg"  | Tumor.suppressor.or.oncogene.prediction..by.20.20.. == "possible tsg")
#remove repeats from this list
driver_tsg <- driver_tsg[!duplicated(driver_tsg$gene),]

#subset driver genes for oncogenes
driver_oncogene <- filter(drivers, Tumor.suppressor.or.oncogene.prediction..by.20.20.. == "oncogene"  | Tumor.suppressor.or.oncogene.prediction..by.20.20.. == "possible oncogene")
#remove repeats from this list
driver_oncogene <- driver_tsg[!duplicated(driver_oncogene$gene),]


### Whole genes

#read in data
genes <- fread("intersect_genes_cnvs.txt")
colnames(genes) <- c("chr", "geneStart", "geneEnd", "ensembl_id", "strand",
                     "gene", "cnv_chr", "cnvStart", "cnvEnd", "cn", "sample", "UKB_id")
# remove all IDs from redacted individuals (negative IDs)
genes <- genes %>%
  filter(UKB_id >= 0)


# get list of all individuals we have cnvs for and add cancer data
all_individuals <- data.frame(UKB_id = unique(genes$UKB_id))
master_df <- left_join(all_individuals, cancer, by = "UKB_id")
# set all NA cancer data to no
master_df$cancer[is.na(master_df$cancer)] <- "no"

## TSGs
#subset genes in cnvs list to tsgs only
tsgs_in_cnvs <- merge(tsgs, genes, by = "gene")
## DELETIONS
# get list of individuals with deletion of TSG
indvs_with_tsg_del <- tsgs_in_cnvs$UKB_id[tsgs_in_cnvs$cn < 2]
# now add to master_df
master_df$tsg_del <- ifelse(master_df$UKB_id %in% indvs_with_tsg_del, "yes", "no")
## DUPLICATIONS
# get list of individuals with deletion of TSG
indvs_with_tsg_dup <- tsgs_in_cnvs$UKB_id[tsgs_in_cnvs$cn > 2]
# now add to master_df
master_df$tsg_dup <- ifelse(master_df$UKB_id %in% indvs_with_tsg_dup, "yes", "no")

## DRIVER TSGs
#subset tsgs in cnvs for driver mutations only
drivers_in_cnvs <- inner_join(driver_tsg, genes, by = "gene", multiple = 'all')
## DELETIONS
# get list of individuals with deletion of TSG
indvs_with_dtsg_del <- drivers_in_cnvs$UKB_id[drivers_in_cnvs$cn < 2]
# now add to master_df
master_df$driver_tsg_del <- ifelse(master_df$UKB_id %in% indvs_with_dtsg_del, "yes", "no")
## DUPLICATIONS
# get list of individuals with deletion of TSG
indvs_with_dtsg_dup <- drivers_in_cnvs$UKB_id[drivers_in_cnvs$cn > 2]
# now add to master_df
master_df$driver_tsg_dup <- ifelse(master_df$UKB_id %in% indvs_with_dtsg_dup, "yes", "no")

## ONCOGENES
#subset genes in cnvs list to oncs only
oncs_in_cnvs <- merge(oncogenes, genes, by = "gene")
## DELETIONS
# get list of individuals with deletion of onc
indvs_with_onc_del <- oncs_in_cnvs$UKB_id[oncs_in_cnvs$cn < 2]
# now add to master_df
master_df$onc_del <- ifelse(master_df$UKB_id %in% indvs_with_onc_del, "yes", "no")
## DUPLICATIONS
# get list of individuals with deletion of onc
indvs_with_onc_dup <- oncs_in_cnvs$UKB_id[oncs_in_cnvs$cn > 2]
# now add to master_df
master_df$onc_dup <- ifelse(master_df$UKB_id %in% indvs_with_onc_dup, "yes", "no")

## DRIVER ONCs
#subset oncs in cnvs for driver mutations only
driver_oncs_in_cnvs <- inner_join(driver_oncogene, genes, by = "gene", multiple = 'all')
## DELETIONS
# get list of individuals with deletion of onc
indvs_with_donc_del <- driver_oncs_in_cnvs$UKB_id[driver_oncs_in_cnvs$cn < 2]
# now add to master_df
master_df$driver_onc_del <- ifelse(master_df$UKB_id %in% indvs_with_donc_del, "yes", "no")
## DUPLICATIONS
# get list of individuals with deletion of TSG
indvs_with_donc_dup <- driver_oncs_in_cnvs$UKB_id[driver_oncs_in_cnvs$cn > 2]
# now add to master_df
master_df$driver_onc_dup <- ifelse(master_df$UKB_id %in% indvs_with_donc_dup, "yes", "no")

# make df to store results
results <- data.frame(
  test = character(0), 
  cancer_cnv = numeric(0),    # cancer yes & deletion yes
  cancer_nocnv = numeric(0),  # cancer yes & deletion no  
  nocancer_cnv = numeric(0),  # cancer no & deletion yes
  nocancer_nocnv = numeric(0), # cancer no & deletion no
  p_value = numeric(0),
  odds_ratio = numeric(0),
  power = numeric(0),
  stringsAsFactors = FALSE
)


###################
### TSG - Deletions

# Make contingency table
contingency_table <- table(deletion = factor(master_df$tsg_del, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Extract counts from the table
cancer_cnv <- contingency_table["yes", "yes"]    # deletions YES & cancer YES
cancer_nocnv <- contingency_table["no", "yes"]   # deletions NO & cancer YES  
nocancer_cnv <- contingency_table["yes", "no"]   # deletions YES & cancer NO
nocancer_nocnv <- contingency_table["no", "no"]  # deletions NO & cancer NO

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_table, alternative = "greater")

## POWER ANALYSIS
# counts
n1 <- sum(contingency_table["yes", ])
n2 <- sum(contingency_table["no", ])
# proportions of cancer
p1 <- contingency_table["yes", "yes"] / n1  # tsg_del = yes
p2 <- contingency_table["no", "yes"] / n2   # tsg_del = no

power.fisher.test(
  p1 = p1,
  p2 = p2,
  n1 = n1,
  n2 = n2,
  alpha = 0.05,
  nsim = 10000,
  alternative = "greater"  # matches your Fisher's test alternative
)

#effect size
p1-p2

# Create new row with results
new_row <- data.frame(
  test = "TSG_deletion",
  cancer_cnv = cancer_cnv,
  cancer_nocnv = cancer_nocnv,
  nocancer_cnv = nocancer_cnv,
  nocancer_nocnv = nocancer_nocnv,
  p_value = fisher_result$p.value,
  odds_ratio = unname(fisher_result$estimate),
  power = power,
  stringsAsFactors = FALSE
)
# Add to dataframe
results <- rbind(results, new_row)


######################
### TSGs - Duplications

# Make contingency table
contingency_table <- table(duplication = factor(master_df$tsg_dup, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Extract counts from the table
cancer_cnv <- contingency_table["yes", "yes"]    # deletions YES & cancer YES
cancer_nocnv <- contingency_table["no", "yes"]   # deletions NO & cancer YES  
nocancer_cnv <- contingency_table["yes", "no"]   # deletions YES & cancer NO
nocancer_nocnv <- contingency_table["no", "no"]  # deletions NO & cancer NO

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_table, alternative = "less")

## POWER ANALYSIS
# Counts per group
n1 <- sum(contingency_table["yes", ])  # deletion = yes
n2 <- sum(contingency_table["no", ])   # deletion = no

# Proportion of cancer cases
p1 <- contingency_table["yes", "yes"] / n1
p2 <- contingency_table["no", "yes"] / n2

power <- power.fisher.test(
  p1 = p1,
  p2 = p2,
  n1 = n1,
  n2 = n2,
  alpha = 0.05,
  nsim = 10000,
  alternative = "less"  # matches your Fisher's test alternative
)

#effect size
p1-p2

# Create new row with results
new_row <- data.frame(
  test = "TSG_duplication",
  cancer_cnv = cancer_cnv,
  cancer_nocnv = cancer_nocnv,
  nocancer_cnv = nocancer_cnv,
  nocancer_nocnv = nocancer_nocnv,
  p_value = fisher_result$p.value,
  odds_ratio = unname(fisher_result$estimate),
  power = power,
  stringsAsFactors = FALSE
)
# Add to dataframe
results <- rbind(results, new_row)



##########################
### Driver TSG - Deletions

# Make contingency table
contingency_table <- table(deletion = factor(master_df$driver_tsg_del, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Extract counts from the table
cancer_cnv <- contingency_table["yes", "yes"]    # deletions YES & cancer YES
cancer_nocnv <- contingency_table["no", "yes"]   # deletions NO & cancer YES  
nocancer_cnv <- contingency_table["yes", "no"]   # deletions YES & cancer NO
nocancer_nocnv <- contingency_table["no", "no"]  # deletions NO & cancer NO

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_table, alternative = "greater")

## POWER ANALYSIS
# counts
n1 <- sum(contingency_table["yes", ])
n2 <- sum(contingency_table["no", ])
# proportions of cancer
p1 <- contingency_table["yes", "yes"] / n1  # tsg_del = yes
p2 <- contingency_table["no", "yes"] / n2   # tsg_del = no

power <- power.fisher.test(
  p1 = p1,
  p2 = p2,
  n1 = n1,
  n2 = n2,
  alpha = 0.05,
  nsim = 10000,
  alternative = "greater"  # matches your Fisher's test alternative
)

#effect size
p1-p2

# Create new row with results
new_row <- data.frame(
  test = "driverTSG_deletion",
  cancer_cnv = cancer_cnv,
  cancer_nocnv = cancer_nocnv,
  nocancer_cnv = nocancer_cnv,
  nocancer_nocnv = nocancer_nocnv,
  p_value = fisher_result$p.value,
  odds_ratio = unname(fisher_result$estimate),
  power = power,
  stringsAsFactors = FALSE
)
# Add to dataframe
results <- rbind(results, new_row)



###############################
### Driver TSG - Duplications

# Make contingency table
contingency_table <- table(duplication = factor(master_df$driver_tsg_dup, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Extract counts from the table
cancer_cnv <- contingency_table["yes", "yes"]    # deletions YES & cancer YES
cancer_nocnv <- contingency_table["no", "yes"]   # deletions NO & cancer YES  
nocancer_cnv <- contingency_table["yes", "no"]   # deletions YES & cancer NO
nocancer_nocnv <- contingency_table["no", "no"]  # deletions NO & cancer NO

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_table, alternative = "less")

## POWER ANALYSIS
# Counts per group
n1 <- sum(contingency_table["yes", ])  # deletion = yes
n2 <- sum(contingency_table["no", ])   # deletion = no

# Proportion of cancer cases
p1 <- contingency_table["yes", "yes"] / n1
p2 <- contingency_table["no", "yes"] / n2

power <- power.fisher.test(
  p1 = p1,
  p2 = p2,
  n1 = n1,
  n2 = n2,
  alpha = 0.05,
  nsim = 10000,
  alternative = "less"  # matches your Fisher's test alternative
)

#effect size
p1-p2

# Create new row with results
new_row <- data.frame(
  test = "driverTSG_duplication",
  cancer_cnv = cancer_cnv,
  cancer_nocnv = cancer_nocnv,
  nocancer_cnv = nocancer_cnv,
  nocancer_nocnv = nocancer_nocnv,
  p_value = fisher_result$p.value,
  odds_ratio = unname(fisher_result$estimate),
  power = power,
  stringsAsFactors = FALSE
)
# Add to dataframe
results <- rbind(results, new_row)



########################
### Oncogenes - Deletions

# Make contingency table
contingency_table <- table(deletion = factor(master_df$onc_del, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Extract counts from the table
cancer_cnv <- contingency_table["yes", "yes"]    # deletions YES & cancer YES
cancer_nocnv <- contingency_table["no", "yes"]   # deletions NO & cancer YES  
nocancer_cnv <- contingency_table["yes", "no"]   # deletions YES & cancer NO
nocancer_nocnv <- contingency_table["no", "no"]  # deletions NO & cancer NO

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_table, alternative = "less")

## POWER ANALYSIS
# counts
n1 <- sum(contingency_table["yes", ])
n2 <- sum(contingency_table["no", ])
# proportions of cancer
p1 <- contingency_table["yes", "yes"] / n1  # tsg_del = yes
p2 <- contingency_table["no", "yes"] / n2   # tsg_del = no

power <- power.fisher.test(
  p1 = p1,
  p2 = p2,
  n1 = n1,
  n2 = n2,
  alpha = 0.05,
  nsim = 10000,
  alternative = "less"  # matches your Fisher's test alternative
)

#effect size
p1-p2

# Create new row with results
new_row <- data.frame(
  test = "ONC_deletion",
  cancer_cnv = cancer_cnv,
  cancer_nocnv = cancer_nocnv,
  nocancer_cnv = nocancer_cnv,
  nocancer_nocnv = nocancer_nocnv,
  p_value = fisher_result$p.value,
  odds_ratio = unname(fisher_result$estimate),
  power = power,
  stringsAsFactors = FALSE
)
# Add to dataframe
results <- rbind(results, new_row)

###########################
### Oncogenes - Duplications

# Make contingency table
contingency_table <- table(deletion = factor(master_df$onc_dup, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Extract counts from the table
cancer_cnv <- contingency_table["yes", "yes"]    # deletions YES & cancer YES
cancer_nocnv <- contingency_table["no", "yes"]   # deletions NO & cancer YES  
nocancer_cnv <- contingency_table["yes", "no"]   # deletions YES & cancer NO
nocancer_nocnv <- contingency_table["no", "no"]  # deletions NO & cancer NO

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_table, alternative = "less")

## POWER ANALYSIS
# counts
n1 <- sum(contingency_table["yes", ])
n2 <- sum(contingency_table["no", ])
# proportions of cancer
p1 <- contingency_table["yes", "yes"] / n1  # tsg_del = yes
p2 <- contingency_table["no", "yes"] / n2   # tsg_del = no

power <- power.fisher.test(
  p1 = p1,
  p2 = p2,
  n1 = n1,
  n2 = n2,
  alpha = 0.05,
  nsim = 10000,
  alternative = "less"  # matches your Fisher's test alternative
)

#effect size
p1-p2

# Create new row with results
new_row <- data.frame(
  test = "ONC_duplication",
  cancer_cnv = cancer_cnv,
  cancer_nocnv = cancer_nocnv,
  nocancer_cnv = nocancer_cnv,
  nocancer_nocnv = nocancer_nocnv,
  p_value = fisher_result$p.value,
  odds_ratio = unname(fisher_result$estimate),
  power = power,
  stringsAsFactors = FALSE
)
# Add to dataframe
results <- rbind(results, new_row)


###############################
### Driver Oncogenes - Deletions

# Make contingency table
contingency_table <- table(deletion = factor(master_df$driver_onc_del, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Extract counts from the table
cancer_cnv <- contingency_table["yes", "yes"]    # deletions YES & cancer YES
cancer_nocnv <- contingency_table["no", "yes"]   # deletions NO & cancer YES  
nocancer_cnv <- contingency_table["yes", "no"]   # deletions YES & cancer NO
nocancer_nocnv <- contingency_table["no", "no"]  # deletions NO & cancer NO

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_table, alternative = "less")

## POWER ANALYSIS
# counts
n1 <- sum(contingency_table["yes", ])
n2 <- sum(contingency_table["no", ])
# proportions of cancer
p1 <- contingency_table["yes", "yes"] / n1  # tsg_del = yes
p2 <- contingency_table["no", "yes"] / n2   # tsg_del = no

power <- power.fisher.test(
  p1 = p1,
  p2 = p2,
  n1 = n1,
  n2 = n2,
  alpha = 0.05,
  nsim = 10000,
  alternative = "less"  # matches your Fisher's test alternative
)

#effect size
p1-p2

# Create new row with results
new_row <- data.frame(
  test = "driverONC_deletion",
  cancer_cnv = cancer_cnv,
  cancer_nocnv = cancer_nocnv,
  nocancer_cnv = nocancer_cnv,
  nocancer_nocnv = nocancer_nocnv,
  p_value = fisher_result$p.value,
  odds_ratio = unname(fisher_result$estimate),
  power = power,
  stringsAsFactors = FALSE
)
# Add to dataframe
results <- rbind(results, new_row)


##################################
### Driver Oncogenes - Duplications

# Make contingency table
contingency_table <- table(deletion = factor(master_df$driver_onc_dup, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Extract counts from the table
cancer_cnv <- contingency_table["yes", "yes"]    # deletions YES & cancer YES
cancer_nocnv <- contingency_table["no", "yes"]   # deletions NO & cancer YES  
nocancer_cnv <- contingency_table["yes", "no"]   # deletions YES & cancer NO
nocancer_nocnv <- contingency_table["no", "no"]  # deletions NO & cancer NO

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_table, alternative = "less")

## POWER ANALYSIS
# counts
n1 <- sum(contingency_table["yes", ])
n2 <- sum(contingency_table["no", ])
# proportions of cancer
p1 <- contingency_table["yes", "yes"] / n1  # tsg_del = yes
p2 <- contingency_table["no", "yes"] / n2   # tsg_del = no

power <- power.fisher.test(
  p1 = p1,
  p2 = p2,
  n1 = n1,
  n2 = n2,
  alpha = 0.05,
  nsim = 10000,
  alternative = "less"  # matches your Fisher's test alternative
)

#effect size
p1-p2

# Create new row with results
new_row <- data.frame(
  test = "driverONC_duplication",
  cancer_cnv = cancer_cnv,
  cancer_nocnv = cancer_nocnv,
  nocancer_cnv = nocancer_cnv,
  nocancer_nocnv = nocancer_nocnv,
  p_value = fisher_result$p.value,
  odds_ratio = unname(fisher_result$estimate),
  power = power,
  stringsAsFactors = FALSE
)
# Add to dataframe
results <- rbind(results, new_row)


#################
##  WRITE RESULTS
write.table(results, "tsg_onc_fisher_results.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)



################################################
########### LOGISTIC REGRESSION
###############################################

## TUMOUR SUPPRESSORS
tsg_summary <- tsgs_in_cnvs %>%
  group_by(UKB_id) %>%
  summarise(
    # Count different types of events
    n_hom_del = sum(cn == 0, na.rm = TRUE),
    n_het_del = sum(cn == 1, na.rm = TRUE),
    n_dup = sum(cn >= 3, na.rm = TRUE),
    n_high_dup = sum(cn >= 4, na.rm = TRUE),
    # Binary indicators
    has_hom_del = any(cn == 0),
    has_het_del = any(cn == 1),
    has_del = any(cn<= 1),
    has_dup = any(cn >= 3),
    
    # Total burden scores
    total_del_burden = sum(cn < 2, na.rm = TRUE),
    total_dup_burden = sum(cn > 2, na.rm = TRUE),
    
    # Genes affected (for sensitivity analyses)
    genes_deleted = paste(unique(gene[cn < 2]), collapse = ";"),
    genes_duplicated = paste(unique(gene[cn > 2]), collapse = ";")
  )

tsg_summary <- left_join(tsg_summary,cancer, by = "UKB_id")
tsg_summary$cancer[is.na(tsg_summary$cancer)] <- "no"
tsg_summary <- left_join(tsg_summary, covariates, by = c("UKB_id"="IID"))
tsg_summary$cancer <- factor(tsg_summary$cancer, levels = c("no", "yes"))


## DELETIONS
model <- glm(cancer ~ has_del+ age + sex + PC1 + PC2 + PC3 + PC4 + 
               batch + smoke + bmi, 
             data = tsg_summary, 
             family = binomial)
summary(model)
# Get the coefficient summary
coef_summary <- summary(model)$coefficients
# Add the row for tsg_del_count
results_df[1, ] <- list(
  test = "tsg_whole_deletions",
  estimate = coef_summary["has_delTRUE", "Estimate"],
  std_error = coef_summary["has_delTRUE", "Std. Error"],
  z_value = coef_summary["has_delTRUE", "z value"],
  p_value = coef_summary["has_delTRUE", "Pr(>|z|)"]
)

## DUPLICATIONS
model <- glm(cancer ~ has_dup + age + sex + PC1 + PC2 + PC3 + PC4 + 
               batch + smoke + bmi, 
             data = tsg_summary, 
             family = binomial)
summary(model)
# Get the coefficient summary
coef_summary <- summary(model)$coefficients
# Add the row for tsg_del_count
results_df[2, ] <- list(
  test = "tsg_whole_duplications",
  estimate = coef_summary["has_dupTRUE", "Estimate"],
  std_error = coef_summary["has_dupTRUE", "Std. Error"],
  z_value = coef_summary["has_dupTRUE", "z value"],
  p_value = coef_summary["has_dupTRUE", "Pr(>|z|)"]
)


driver_tsg_summary <- drivers_in_cnvs %>%
  group_by(UKB_id) %>%
  summarise(
    # Count different types of events
    n_hom_del = sum(cn == 0, na.rm = TRUE),
    n_het_del = sum(cn == 1, na.rm = TRUE),
    n_dup = sum(cn >= 3, na.rm = TRUE),
    n_high_dup = sum(cn >= 4, na.rm = TRUE),
    # Binary indicators
    has_hom_del = any(cn == 0),
    has_het_del = any(cn == 1),
    has_del = any(cn<= 1),
    has_dup = any(cn >= 3),
    
    # Total burden scores
    total_del_burden = sum(cn < 2, na.rm = TRUE),
    total_dup_burden = sum(cn > 2, na.rm = TRUE),
    
    # Genes affected (for sensitivity analyses)
    genes_deleted = paste(unique(gene[cn < 2]), collapse = ";"),
    genes_duplicated = paste(unique(gene[cn > 2]), collapse = ";")
  )


driver_tsg_summary <- left_join(driver_tsg_summary,cancer, by = "UKB_id")
driver_tsg_summary$cancer[is.na(driver_tsg_summary$cancer)] <- "no"
driver_tsg_summary <- left_join(driver_tsg_summary, covariates, by = c("UKB_id"="IID"))
driver_tsg_summary$cancer <- factor(driver_tsg_summary$cancer, levels = c("no", "yes"))


## DELETIONS
model <- glm(cancer ~ has_del + age + sex + PC1 + PC2 + PC3 + PC4 + 
               batch + smoke + bmi, 
             data = driver_tsg_summary, 
             family = binomial)
summary(model)
# Get the coefficient summary
coef_summary <- summary(model)$coefficients
# Add the row for tsg_del_count
results_df[3, ] <- list(
  test = "drivertsg_whole_deletions",
  estimate = coef_summary["has_delTRUE", "Estimate"],
  std_error = coef_summary["has_delTRUE", "Std. Error"],
  z_value = coef_summary["has_delTRUE", "z value"],
  p_value = coef_summary["has_delTRUE", "Pr(>|z|)"]
)

## DUPLICATIONS
model <- glm(cancer ~ has_dup + age + sex + PC1 + PC2 + PC3 + PC4 + 
               batch + smoke + bmi, 
             data = driver_tsg_summary, 
             family = binomial)
summary(model)
# Get the coefficient summary
coef_summary <- summary(model)$coefficients
# Add the row for tsg_del_count
results_df[4, ] <- list(
  test = "drivertsg_whole_duplications",
  estimate = coef_summary["has_dupTRUE", "Estimate"],
  std_error = coef_summary["has_dupTRUE", "Std. Error"],
  z_value = coef_summary["has_dupTRUE", "z value"],
  p_value = coef_summary["has_dupTRUE", "Pr(>|z|)"]
)



#### ONCOGENES
oncs_summary <- oncs_in_cnvs %>%
  group_by(UKB_id) %>%
  summarise(
    # Count different types of events
    n_hom_del = sum(cn == 0, na.rm = TRUE),
    n_het_del = sum(cn == 1, na.rm = TRUE),
    n_dup = sum(cn >= 3, na.rm = TRUE),
    n_high_dup = sum(cn >= 4, na.rm = TRUE),
    # Binary indicators
    has_hom_del = any(cn == 0),
    has_het_del = any(cn == 1),
    has_del = any(cn<= 1),
    has_dup = any(cn >= 3),
    
    # Total burden scores
    total_del_burden = sum(cn < 2, na.rm = TRUE),
    total_dup_burden = sum(cn > 2, na.rm = TRUE),
    
    # Genes affected (for sensitivity analyses)
    genes_deleted = paste(unique(gene[cn < 2]), collapse = ";"),
    genes_duplicated = paste(unique(gene[cn > 2]), collapse = ";")
  )

oncs_summary <- left_join(oncs_summary,cancer, by = "UKB_id")
oncs_summary$cancer[is.na(oncs_summary$cancer)] <- "no"
oncs_summary <- left_join(oncs_summary, covariates, by = c("UKB_id"="IID"))
oncs_summary$cancer <- factor(oncs_summary$cancer, levels = c("no", "yes"))

## DELETIONS
model <- glm(cancer ~ has_del + age + sex + PC1 + PC2 + PC3 + PC4 + 
               batch + smoke + bmi, 
             data = oncs_summary, 
             family = binomial)
summary(model)
# Get the coefficient summary
coef_summary <- summary(model)$coefficients
# Add the row for tsg_del_count
results_df[5, ] <- list(
  test = "onc_whole_deletions",
  estimate = coef_summary["has_delTRUE", "Estimate"],
  std_error = coef_summary["has_delTRUE", "Std. Error"],
  z_value = coef_summary["has_delTRUE", "z value"],
  p_value = coef_summary["has_delTRUE", "Pr(>|z|)"]
)

## DUPLICATIONS
model <- glm(cancer ~ has_dup + age + sex + PC1 + PC2 + PC3 + PC4 + 
               batch + smoke + bmi, 
             data = tsg_summary, 
             family = binomial)
summary(model)
# Get the coefficient summary
coef_summary <- summary(model)$coefficients
# Add the row for tsg_del_count
results_df[6, ] <- list(
  test = "onc_whole_duplications",
  estimate = coef_summary["has_dupTRUE", "Estimate"],
  std_error = coef_summary["has_dupTRUE", "Std. Error"],
  z_value = coef_summary["has_dupTRUE", "z value"],
  p_value = coef_summary["has_dupTRUE", "Pr(>|z|)"]
)


driver_oncs_summary <- driver_oncs_in_cnvs %>%
  group_by(UKB_id) %>%
  summarise(
    # Count different types of events
    n_hom_del = sum(cn == 0, na.rm = TRUE),
    n_het_del = sum(cn == 1, na.rm = TRUE),
    n_dup = sum(cn >= 3, na.rm = TRUE),
    n_high_dup = sum(cn >= 4, na.rm = TRUE),
    # Binary indicators
    has_hom_del = any(cn == 0),
    has_het_del = any(cn == 1),
    has_del = any(cn<= 1),
    has_dup = any(cn >= 3),
    
    # Total burden scores
    total_del_burden = sum(cn < 2, na.rm = TRUE),
    total_dup_burden = sum(cn > 2, na.rm = TRUE),
    
    # Genes affected (for sensitivity analyses)
    genes_deleted = paste(unique(gene[cn < 2]), collapse = ";"),
    genes_duplicated = paste(unique(gene[cn > 2]), collapse = ";")
  )


driver_oncs_summary <- left_join(driver_oncs_summary,cancer, by = "UKB_id")
driver_oncs_summary$cancer[is.na(driver_oncs_summary$cancer)] <- "no"
driver_oncs_summary <- left_join(driver_oncs_summary, covariates, by = c("UKB_id"="IID"))
driver_oncs_summary$cancer <- factor(driver_oncs_summary$cancer, levels = c("no", "yes"))

## DELETIONS
model <- glm(cancer ~ has_del + age + sex + PC1 + PC2 + PC3 + PC4 + 
               batch + smoke + bmi, 
             data = driver_oncs_summary, 
             family = binomial)
summary(model)
# Get the coefficient summary
coef_summary <- summary(model)$coefficients
# Add the row for tsg_del_count
results_df[7, ] <- list(
  test = "driveronc_whole_deletions",
  estimate = coef_summary["has_delTRUE", "Estimate"],
  std_error = coef_summary["has_delTRUE", "Std. Error"],
  z_value = coef_summary["has_delTRUE", "z value"],
  p_value = coef_summary["has_delTRUE", "Pr(>|z|)"]
)

## DUPLICATIONS
model <- glm(cancer ~ has_dup + age + sex + PC1 + PC2 + PC3 + PC4 + 
               batch + smoke + bmi, 
             data = driver_oncs_summary, 
             family = binomial)
summary(model)
# Get the coefficient summary
coef_summary <- summary(model)$coefficients
# Add the row for tsg_del_count
results_df[8, ] <- list(
  test = "driveronc_whole_duplications",
  estimate = coef_summary["has_dupTRUE", "Estimate"],
  std_error = coef_summary["has_dupTRUE", "Std. Error"],
  z_value = coef_summary["has_dupTRUE", "z value"],
  p_value = coef_summary["has_dupTRUE", "Pr(>|z|)"]
)



#### RESULTS
write.table(results_df, "~/Documents/PhD/TSGs_oncogenes/Data/log_reg_results.txt", 
            col.names=TRUE, row.names = FALSE, quote = FALSE)
