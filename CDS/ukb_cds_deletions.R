
library(data.table)
library(dplyr)
library(tidyr)

setwd("/data4/smatthews/pheWAS/CDS/")

### SETUP
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

# read in gene cds overlapping with cnvs
genes <- fread("intersect_cds_cnvs.txt")
colnames(genes) <- c("chr", "cdsStart", "cdsEnd", "gene", "transcript_id",
                     "exon_num", "V5", "strand", "cnv_chr", "cnvStart",
                     "cnvEnd", "cn", "sample", "UKB_id")

#subset to deletions
genes <- genes %>% filter(cn <= 2)

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
tsgs_in_cnvs <- merge(tsgs, genes, by = "gene", all.x = FALSE)
## DELETIONS
# get list of individuals with deletion of TSG
indvs_with_tsg_del <- tsgs_in_cnvs$UKB_id[tsgs_in_cnvs$cn < 2]
# now add to master_df
master_df$tsg_del <- ifelse(master_df$UKB_id %in% indvs_with_tsg_del, "yes", "no")

## DRIVER TSGs
#subset tsgs in cnvs for driver mutations only
drivers_in_cnvs <- inner_join(driver_tsg, genes, by = "gene", multiple = 'all')
## DELETIONS
# get list of individuals with deletion of TSG
indvs_with_dtsg_del <- drivers_in_cnvs$UKB_id[drivers_in_cnvs$cn < 2]
# now add to master_df
master_df$driver_tsg_del <- ifelse(master_df$UKB_id %in% indvs_with_dtsg_del, "yes", "no")


## ONCOGENES
#subset genes in cnvs list to oncs only
oncs_in_cnvs <- merge(oncogenes, genes, by = "gene")
## DELETIONS
# get list of individuals with deletion of onc
indvs_with_onc_del <- oncs_in_cnvs$UKB_id[oncs_in_cnvs$cn < 2]
# now add to master_df
master_df$onc_del <- ifelse(master_df$UKB_id %in% indvs_with_onc_del, "yes", "no")


## DRIVER ONCs
#subset oncs in cnvs for driver mutations only
driver_oncs_in_cnvs <- inner_join(driver_oncogene, genes, by = "gene", multiple = 'all')
## DELETIONS
# get list of individuals with deletion of onc
indvs_with_donc_del <- driver_oncs_in_cnvs$UKB_id[driver_oncs_in_cnvs$cn < 2]
# now add to master_df
master_df$driver_onc_del <- ifelse(master_df$UKB_id %in% indvs_with_donc_del, "yes", "no")


########### TSG DELETIONS ############
# Make contingency table
contingency_table <- table(deletion = factor(master_df$tsg_del, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Run Fisher's exact test
fisher.test(contingency_table, alternative = "greater")
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


########### TSG DRIVER DELETIONS ############
# Make contingency table
contingency_table <- table(deletion = factor(master_df$driver_tsg_del, levels = c("yes","no")), 
                           cancer = factor(master_df$cancer, levels = c("yes","no")))
# Run Fisher's exact test
fisher.test(contingency_table, alternative = "greater")

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
