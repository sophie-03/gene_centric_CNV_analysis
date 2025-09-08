#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --output=PennCNV.out.txt #output file
#SBATCH --error=PennCNV.err.txt #error file
#SBATCH --job-name=PennCNV #job name
#SBATCH -D . # set working directory to
#SBATCH --partition=highmem
#SBATCH --ntasks=5

cd ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2

#create pfb file
perl ~/bin/PennCNV-1.0.5/compile_pfb.pl -listfile ~/../../data2/smatthews/UKB/chr2/individual_files/sample_list.txt -output chr2.pfb

##run PennCNV in batches of 100,000
module load Anaconda3/4.4.0
source activate PennCNV
#need to redo this whithin the conda environment
export PATH=/home/dbennett/bin/perl-5.18.4/:$PATH
#cut sample-list into 5 lists of 100,000
cd ~/../../data2/smatthews/UKB/chr2
split -l 100000 /data2/smatthews/UKB/chr2/individual_files/sample_list.txt
#move to PennCNV directory
cd ~/bin/PennCNV-1.0.5
parallel < /data2/smatthews/UKB/scripts/PennCNV_parallel.sh

#combine all call files
cd ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2
cat calls1.rawcnv calls2.rawcnv calls3.rawcnv calls4.rawcnv calls5.rawcnv > chr2_calls.rawcnv
