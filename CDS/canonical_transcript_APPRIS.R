library(ape)
library(dplyr)
library(tidyr)
library(stringr)

#read in gff
gff <- read.gff("/data4/smatthews/pheWAS/CDS/gencode_grch37.gff3")

# keep only transcripts
tx <- gff %>% filter(type == "transcript")

#seperate the attributes column
tx_expanded <- tx %>%
  separate_rows(attributes, sep = ";") %>%
  separate(attributes, into = c("key", "value"), sep = "=", extra = "merge") %>%
  pivot_wider(
    names_from = key, 
    values_from = value,
    values_fn = ~ paste(unique(.x), collapse = ";")  # Combine duplicates
  )

# Now filter for protein-coding genes AND transcripts
protein_coding_tx <- tx_expanded %>%
  filter(gene_type == "protein_coding",
         transcript_type == "protein_coding")  # Double-check with transcript_type

# Get APPRIS principal 1 transcripts 
appris_principal <- protein_coding_tx %>%
  filter(grepl("appris_principal_1", tag))

#do any genes have multiple transcripts tagged at appris_principal_1?
nrow(appris_principal) #20161
length(unique(appris_principal$gene_name)) #14553

# Select one transcript per gene with appris_principal_1, by selecting for basic tag and then length
final_canonical <- appris_principal %>%
  group_by(gene_id) %>%
  arrange(
    desc(grepl("basic", tag)),  # Prefer basic tag
    desc(end - start)           # Then longest transcript
  ) %>%
  slice(1) %>%
  ungroup()

# Identify genes without appris principal 1 tag
genes_without_appris <- setdiff(
  unique(protein_coding_tx$gene_id),
  unique(final_canonical$gene_id)
)


# Use other APPRIS levels for remaining genes
if (length(genes_without_appris) > 0) {
  appris_backup <- protein_coding_tx %>%
    filter(gene_id %in% genes_without_appris) %>%
    # Look for any APPRIS annotation
    filter(grepl("appris_principal_[0-9]", tag)) %>%
    # Extract APPRIS level
    mutate(appris_level = as.numeric(str_extract(tag, "appris_principal_([0-9])", group = 1))) %>%
    group_by(gene_id) %>%
    # Take best APPRIS level, then basic, then longest
    arrange(appris_level, desc(grepl("basic", tag)), desc(end - start)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-appris_level)
  
  # Combine
  all_canonical <- bind_rows(final_canonical, appris_backup)
} else {
  all_canonical <- final_canonical
}

# FINAL ENFORCEMENT: Ensure exactly one per gene
all_canonical <- all_canonical %>%
  group_by(gene_id) %>%
  slice(1) %>%  # Take first if any duplicates remain
  ungroup()

# Verification
cat("\n=== FINAL VERIFICATION ===\n")
cat("Total canonical transcripts:", nrow(all_canonical), "\n")
cat("Unique gene_ids:", length(unique(all_canonical$gene_id)), "\n")
cat("These numbers MUST be equal!\n")

#write file
write.table(all_canonical,"/data4/smatthews/pheWAS/CDS/primary_transcripts_APPRIS.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

#make bed file
bed <- all_canonical %>%
  mutate(
    start_bed = start - 1,  # Convert to 0-based for BED
    ID = transcript_id       # Use transcript_id as ID
  ) %>%
  select(
    seqid,           # Chromosome
    start = start_bed, # 0-based start
    end,             # 1-based end (same as GFF end)
    gene_name,       # Gene name
    ID               # Transcript ID
  )

#write bed file
write.table(bed, "/data4/smatthews/pheWAS/CDS/canonical_transcripts.bed", sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

## GET CDS FOR THESE TRANSCRIPTS
