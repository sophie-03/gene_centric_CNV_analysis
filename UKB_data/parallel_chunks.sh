#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --output=out.txt #output file
#SBATCH --error=err.txt #error file
#SBATCH --job-name=par_chunks #job name
#SBATCH -D . # set working directory to
#SBATCH --partition=highmem
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2

cd ~/../../data2/smatthews/UKB/scripts
module load R/R-4.0.2
module load Anaconda3/4.4.0
source activate PennCNV

parallel --jobs 5 Rscript ::: indv_files_1-10.R indv_files_11-20.R indv_files_21-30.R indv_files_31-40.R indv_files_41-49.R
