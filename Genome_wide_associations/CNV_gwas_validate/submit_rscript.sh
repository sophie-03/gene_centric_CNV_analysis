#!/bin/bash
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=s.matthews5@universityofgalway.ie
#SBATCH --job-name=RO_0.3_regions
#SBATCH -D /data4/smatthews/pheWAS/cnv_GWAS
#SBATCH -o logs/cnvrangesRO.out
#SBATCH --cpus-per-task=16

module load Anaconda3/2024.02-1
conda activate cnvGWAS

Rscript /data4/smatthews/pheWAS/github_gene_centric_cnv_analysis/Genome_wide_associations/CNV_gwas_validate/CNVRanges_RO_dels.R
