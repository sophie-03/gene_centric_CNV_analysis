library(CNVRanger)
library(GenomicRanges)
library(data.table)
library(dplyr)

#read in cnv calls
cnvs <- read.table("/data4/smatthews/pheWAS/cnv_GWAS/cnv_calls.bed")
colnames(cnvs) <- c("chr","start","end","state","sample_id","UKB_id")

#subset to deletions only
cnvs <- cnvs[cnvs$state < 2, ]

#group the calls by sample ID and convert them to a GRangesList
grl <- GenomicRanges::makeGRangesListFromDataFrame(cnvs, 
                                            split.field="sample_id", keep.extra.columns=TRUE)
#sort
grl <- GenomicRanges::sort(grl)

#summarize cnv calls to cnv regions
cnvrs <- populationRanges(grl, mode = "RO", ro.thresh=0.5)

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

write.table(cnvs_with_regions, "/data4/smatthews/pheWAS/cnv_GWAS/cnv_regions_dels_RO_0.5.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)



## SPLIT CNV REGIONS INTO DIFFERENT FILES FOR PROCESSING
# Extract numeric region index from region_id
cnvs_with_regions[, region_num := as.integer(gsub("region_", "", region_id))]

# Define bins of 10 regions each
cnvs_with_regions[, region_group := ceiling(region_num / 20)]

# Get total number of groups (should be around 71 for ~710 regions)
n_groups <- max(cnvs_with_regions$region_group)

# Split into list of data.tables
cnv_subsets <- split(cnvs_with_regions, by = "region_group", keep.by = FALSE)

# Optionally, assign each subset to a variable in the environment (cnvrs1, cnvrs2, ...)
for (i in seq_along(cnv_subsets)) {
  assign(paste0("cnvrs", i), cnv_subsets[[i]])
}

# Write each subset to a separate file
for (i in seq_along(cnv_subsets)) {
   out_file <- paste0("/data4/smatthews/pheWAS/cnv_GWAS/cnv_regions_RO_0.5_del_split/cnv_regions_RO_0.5_del_part", i, ".txt")
   fwrite(cnv_subsets[[i]], out_file, sep = "\t")
 }
