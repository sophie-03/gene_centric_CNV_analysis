#!/bin/bash
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=s.matthews5@universityofgalway.ie
#SBATCH --job-name=longest_cds
#SBATCH -D /data4/smatthews/pheWAS/CDS

module load Anaconda3/2024.02-1
conda activate cnvGWAS

cd /data4/smatthews/pheWAS/CDS

## use AGAT
 agat_sp_keep_longest_isoform.pl --gff gencode_grch37.gff3 -o longest_isoform.gff3
