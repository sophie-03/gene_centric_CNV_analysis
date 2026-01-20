#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --mail-type=ALL # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --nice=10000
#SBATCH -D /data4/smatthews/pheWAS/CDS

module load Anaconda3/2024.02-1
conda activate cnvGWAS

cd /data4/smatthews/pheWAS/CDS/

# Get the command line argument for the PHENOTYPE variable
PHENOTYPE=$1

#create genes (only needs to be created once)
# awk '{print $2}' /data4/smatthews/pheWAS/CDS/cds_summary.txt | sort -u > /data4/smatthews/pheWAS/CDS/genes.txt

## Define parallel groups
num_groups=100
# Count the total number of lines in the genes.txt file
total_lines=$(wc -l < /data4/smatthews/pheWAS/CDS/genes.txt)
# Calculate the number of lines per group (rounded up)
lines_per_group=$(( (total_lines + num_groups - 1) / num_groups ))
# Create temporary directories for each group
for ((i=1; i<=num_groups; i++))
do
    mkdir -p "group_$i"
done
# Distribute genes into different groups
awk -v num_groups="$num_groups" -v lines_per_group="$lines_per_group" '
    { group_id = int((NR - 1) / lines_per_group) + 1; print > "group_" group_id "/genes.txt" }
' /data4/smatthews/pheWAS/CDS/genes.txt


# run 
while read -r gene
do
    # Create a temp file for the current group
    grep -w "$gene" /data4/smatthews/pheWAS/CDS/cds_summary.txt > "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.rawcnv"

    # Process the boundary with the R script
    Rscript /data4/smatthews/pheWAS/github_gene_centric_cnv_analysis/CDS/cds_assoc_dels.R "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.rawcnv" "$SLURM_ARRAY_TASK_ID" "${PHENOTYPE}"

    # Clean up the temp_cnv.rawcnv file
    rm -f "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.rawcnv"
done < "group_$SLURM_ARRAY_TASK_ID/genes.txt"

#concatenate results into one file
cat group_*/temp_logistic.txt > "${PHENOTYPE}"_dels_logistic_output.txt

# clean up group directories
#rm -r group_*
