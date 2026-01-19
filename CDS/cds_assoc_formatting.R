
#setwd
setwd("/data4/smatthews/pheWAS/CDS/")

#libraries
library(dplyr)
library(data.table)

# read in cds data
cds <- fread("intersect_cds_cnvs.txt")
colnames(cds) <- c("chr", "cdsStart", "cdsEnd", "gene", "transcript_id",
                     "exon_num", "V5", "strand", "cnv_chr", "cnvStart",
                     "cnvEnd", "cn", "sample", "UKB_id")

#subset to deletions
cds <- cds %>% filter(cn <= 2)

# remove all IDs from redacted individuals (negative IDs)
cds <- cds %>%
  filter(UKB_id >= 0)

# make new df where each gene has only one row per individual (concatenating individuals where deletion spans multiple CDS)
# add column to say that there is a deletion overlepping that gene cds for this individual
formatted_df <- cds[, .(has_cds_deletion = 1L), 
                    by = .(UKB_id, gene, chr)]

# write table
write.table(formatted_df, "cds_summary.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
