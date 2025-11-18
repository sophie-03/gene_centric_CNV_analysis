#!/bin/bash
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=s.matthews5@universityofgalway.ie
#SBATCH --job-name=D_0.05_regions
#SBATCH -D /data4/smatthews/pheWAS/cnv_GWAS
#SBATCH -o logs/cnvrangesD_0.05.out
#SBATCH --cpus-per-task=16

module load Anaconda3/2024.02-1
conda activate cnvGWAS

Rscript /data4/smatthews/pheWAS/github_gene_centric_cnv_analysis/Genome_wide_associations/CNV_gwas_validate/CNVRanges_dels.R
