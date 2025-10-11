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

write.table(cnvs_with_regions, "/data4/smatthews/pheWAS/cnv_GWAS/cnv_regions.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
