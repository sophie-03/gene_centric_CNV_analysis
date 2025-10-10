#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --job-name=C44_dups #job name
#SBATCH -D /data4/smatthews/pheWAS/cnv_GWAS # set working directory to

module load Anaconda3/2024.02-1
conda activate cnvGWAS

Rscript /data4/smatthews/pheWAS/github_gene_centric_cnv_analysis/Genome_wide_associations/CNV_gwas_validate/CNV_GWAS_dups.R C44
