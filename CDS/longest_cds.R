library(ape)
library(dplyr)

#read in longest isoform gff
longest_isoform <- read.gff("~/Downloads/longest_isoform.gff3")

# filter CDS
cds <- longest_isoform %>% filter(type == "CDS")

# extract IDs
cds <- cds %>%
  mutate(
    gene_id = str_extract(attributes, "(?<=gene_id=)[^;]+"),
    gene_name = str_extract(attributes, "(?<=gene_name=)[^;]+"),
    transcript_id = str_extract(attributes, "(?<=transcript_id=)[^;]+"),
    cds_len = end - start + 1
  )

cds_summary <- cds %>%
  group_by(gene_id, gene_name, transcript_id, seqid, strand) %>%
  summarise(
    cds_start = min(start),
    cds_end   = max(end),
    total_cds_length = sum(cds_len),
    .groups = "drop"
  )

bed <- cds_summary %>%
  select(seqid, cds_start, cds_end, gene_name, total_cds_length, strand)

write.table(
  bed,
  file = "~/Downloads/longest_cds.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
