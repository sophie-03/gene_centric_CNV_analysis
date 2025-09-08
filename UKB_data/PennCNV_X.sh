#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --output=PennCNV_X.out.txt #output file
#SBATCH --error=PennCNV_X.err.txt #error file
#SBATCH --job-name=PennCNV #job name
#SBATCH -D . # set working directory to
#SBATCH --partition=normal
#SBATCH --ntasks=5

cd ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX

#make snp position file
#echo -e "Name\tChr\tPosition" > "snp_position_chrX2.txt"
#paste <(cut -f2,1,4 ~/../../data2/smatthews/UKB/bim_files/ukb_snp_chrX2_v2.bim) >> "snp_position_chrX2.txt"

#create pfb file
perl ~/bin/PennCNV-1.0.5/compile_pfb.pl -listfile ~/../../data2/smatthews/UKB/chrX/individual_files/sample_list.txt -output chrX.pfb

##run PennCNV in batches of 100,000
module load Anaconda3/4.4.0
source activate PennCNV
#need to redo this whithin the conda environment
export PATH=/home/dbennett/bin/perl-5.18.4/:$PATH
#cut sample-list into 5 lists of 100,000
cd ~/../../data2/smatthews/UKB/chrX
split -l 100000 /data2/smatthews/UKB/chrX/individual_files/sample_list.txt
#move to PennCNV directory
cd ~/bin/PennCNV-1.0.5
parallel < /data2/smatthews/UKB/scripts/PennCNV_parallel_X.sh

#combine all call files
cd ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX
cat calls1.rawcnv calls2.rawcnv calls3.rawcnv calls4.rawcnv calls5.rawcnv > chrX_calls.rawcnv
