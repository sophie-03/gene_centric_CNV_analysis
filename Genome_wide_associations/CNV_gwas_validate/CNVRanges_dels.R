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
cnvrs <- populationRanges(grl, density=0.05)

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

write.table(cnvs_with_regions, "/data4/smatthews/pheWAS/cnv_GWAS/cnv_regions_dels_density0.05.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)



## SPLIT CNV REGIONS INTO DIFFERENT FILES FOR PROCESSING
# Get unique region IDs
unique_regions <- unique(cnvs_with_regions$region_id)

# Split into 50 roughly equal chunks
n_chunks <- 50
region_chunks <- split(unique_regions, cut(seq_along(unique_regions), n_chunks, labels = FALSE))

#  Subset the original dataframe for each chunk
df_list <- lapply(region_chunks, function(regs) {
  subset(cnvs_with_regions, region_id %in% regs)
})

# Make sure the output directory exists
out_dir <- "/data4/smatthews/pheWAS/cnv_GWAS/cnv_regions_del_density0.05_split"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Loop through the list and write each dataframe
for(i in seq_along(df_list)) {
  out_file <- file.path(out_dir, paste0("cnv_regions_del_density0.05_part", i, ".txt"))
  fwrite(df_list[[i]], out_file, sep = "\t")
}


