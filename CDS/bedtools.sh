#!/bin/bash
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=s.matthews5@universityofgalway.ie
#SBATCH --job-name=bedtools
#SBATCH -D /data4/smatthews/pheWAS/CDS

module load Anaconda3/2024.02-1
conda activate cnvGWAS

cd /data4/smatthews/pheWAS/CDS/

# get genes where the longest cds is overlapped by a deletion
bedtools intersect -wa -wb -a canonical_transcripts_CDS.bed -b cnv_calls.bed -filenames > intersect_cds_cnvs.txt
