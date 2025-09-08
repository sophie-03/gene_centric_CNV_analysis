#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --output=out.txt #output file
#SBATCH --error=err.txt #error file
#SBATCH --job-name=format #job name
#SBATCH -D . # set working directory to
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=2

cd /data2/smatthews/UKB/chrX/individual_files

#replace tab-23-tab with tab-X-tab (this way we are only targeting 23s that are stand alones in their own column)
sed -i 's/\t23\t/\tX\t/g' sample1*.txt
sed -i 's/\t23\t/\tX\t/g' sample2*.txt
sed -i 's/\t23\t/\tX\t/g' sample3*.txt
sed -i 's/\t23\t/\tX\t/g' sample4*.txt
sed -i 's/\t23\t/\tX\t/g' sample5*.txt
sed -i 's/\t23\t/\tX\t/g' sample6*.txt
sed -i 's/\t23\t/\tX\t/g' sample7*.txt
sed -i 's/\t23\t/\tX\t/g' sample8*.txt
sed -i 's/\t23\t/\tX\t/g' sample9*.txt
