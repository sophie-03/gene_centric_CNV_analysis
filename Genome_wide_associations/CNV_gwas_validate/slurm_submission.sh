#!/bin/bash
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=s.matthews5@universityofgalway.ie
#SBATCH --job-name=C50_dels
#SBATCH -D /data4/smatthews/pheWAS/cnv_GWAS
#SBATCH --array=1-40         # one per CNV part file
#SBATCH -o logs/C50_dels_%A_%a.out
#SBATCH --nice=100

module load Anaconda3/2024.02-1
conda activate cnvGWAS

PHENO="C50"
PART=${SLURM_ARRAY_TASK_ID}

Rscript /data4/smatthews/pheWAS/github_gene_centric_cnv_analysis/Genome_wide_associations/CNV_gwas_validate/CNV_GWAS_dels.R $PHENO $PART
