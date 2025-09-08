# make covariate file for cnv association test

setwd("/data4/smatthews/pheWAS")

library(data.table)
library(dplyr)

#read in list of all ukb IDs
all_IDs <- fread("map_file.txt", header = TRUE, select = c(1))

# covariates: age, sex, PCs 1-4, batch

#sex
sex <- fread("ukb22431_c22_b0_v2_s488146.fam", select = c(1, 5))
colnames(sex) <- c("UKB_id", "sex")

#age
age <- read.table("21022.txt", header = TRUE, fill = TRUE)
colnames(age) <- c("UKB_id", "age")

# PC
PC <- fread("22009.txt", header = TRUE, fill = TRUE, sep = "\t")
#remove empty first column
PC <- PC[,c(-1)]
#keep only the id col + the first 4 PCs
PC <- PC[,1:5]
#rename columns
colnames(PC) <- c("UKB_id", "PC1", "PC2", "PC3", "PC4")

#batch
batch <- fread("ukb22431_c22_b0_v2_s488146.fam", select = c(1,6))
colnames(batch) <- c("UKB_id", "batch")

#smoking
smoking <- fread("20116.txt", header = FALSE, select = c(2,3))
colnames(smoking) <- c("UKB_id", "smoke")
# change status of smoking (0=never, 1=used to smoke, 2=currently smokes)
# Assuming your dataframe is named 'smoking'
smoking <- smoking %>%
  mutate(smoke = case_when(
    smoke == 1 ~ "previous",
    smoke == 0 ~ "never",
    smoke == 2 ~ "current",
    smoke == -3 ~ "no_data",
    is.na(smoke) ~ "no_data",
    TRUE ~ as.character(smoke)
  ))

#merge covariates into one df
covariates <- data.frame(all_IDs$UKB_id)
colnames(covariates) <- "UKB_id"
covariates <- left_join(covariates, age, by = c("UKB_id"))
covariates <- left_join(covariates, sex, by = c("UKB_id"))
covariates <- left_join(covariates, PC, by = c("UKB_id"))
covariates <- left_join(covariates, batch, by = c("UKB_id"))
covariates <- left_join(covariates, smoking, by = c("UKB_id"))

#remove rows for IDs that have been redacted (- numbers)
covariates <- subset(covariates, UKB_id >= 0)

#change UKB_id to IID for plink
names(covariates)[names(covariates) == "UKB_id"] <- "IID"

# export data
write.table(covariates, "covariates.txt", sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = TRUE)
