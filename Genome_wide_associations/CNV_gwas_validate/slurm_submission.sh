#!/bin/bash
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=s.matthews5@universityofgalway.ie
#SBATCH --job-name=C61_dels
#SBATCH -D /data4/smatthews/pheWAS/cnv_GWAS
#SBATCH --array=1-40         # one per CNV part file
#SBATCH -o logs/C61_dels_%A_%a.out
#SBATCH --nice=100
#SBATCH --partition=normal
#SBATCH --exclude=cn070,cn071,cn072,cn073,cn074,cn075,cn076,cn100
#SBATCH --cpus-per-task=16

module load Anaconda3/2024.02-1
conda activate cnvGWAS

PHENO="C61"
PART=${SLURM_ARRAY_TASK_ID}

Rscript /data4/smatthews/pheWAS/github_gene_centric_cnv_analysis/Genome_wide_associations/CNV_gwas_validate/CNV_GWAS_dels.R $PHENO $PART
