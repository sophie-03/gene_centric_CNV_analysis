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

# Create scoring system for ALL transcripts at once
all_canonical <- protein_coding_tx %>%
  # Create priority scores
  mutate(
    # Priority 1: MANE Select (highest clinical relevance)
    mane_score = ifelse(grepl("MANE_Select", tag), 1, 2),
    
    # Priority 2: APPRIS level (1 is best)
    appris_score = case_when(
      grepl("appris_principal_1", tag) ~ 1,
      grepl("appris_principal_2", tag) ~ 2,
      grepl("appris_principal_3", tag) ~ 3,
      grepl("appris_principal_4", tag) ~ 4,
      grepl("appris_principal_5", tag) ~ 5,
      TRUE ~ 6
    ),
    
  ) %>%
  
  # Group by gene and select best transcript
  group_by(gene_id) %>%
  arrange(
    mane_score,      # MANE first
    appris_score,    # Then best APPRIS
  ) %>%
  slice(1) %>%
  ungroup() %>%
  
  # Remove scoring columns
  select(-mane_score, -appris_score)


#write file
write.table(all_canonical,"/data4/smatthews/pheWAS/CDS/primary_transcripts_MANE.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

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
write.table(bed, "/data4/smatthews/pheWAS/CDS/canonical_transcripts_MANE.bed", sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

## GET CDS FOR THESE TRANSCRIPTS

# Get the canonical transcript IDs
canonical_tx_ids <- all_canonical$transcript_id

# Filter the original GFF for CDS features of these transcripts
cds_features <- gff %>%
  filter(type == "CDS") %>%
  # Extract transcript_id from attributes
  mutate(
    transcript_id = str_extract(attributes, "transcript_id=([^;]+)", group = 1)
  ) %>%
  filter(transcript_id %in% canonical_tx_ids) %>%
  # Extract additional useful info
  mutate(
    gene_name = str_extract(attributes, "gene_name=([^;]+)", group = 1),
    exon_number = str_extract(attributes, "exon_number=([^;]+)", group = 1),
    protein_id = str_extract(attributes, "protein_id=([^;]+)", group = 1)
  )
