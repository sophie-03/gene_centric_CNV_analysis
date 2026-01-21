#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --mail-type=ALL # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --array=1-100

conda activate phewas

## get command line arguments
phenotype_file=$1
overlap_file=$2
covariate_file=$3

## create list of genes
awk '{print $6}' intersect_genes_cnvs.txt | sort -u > genes.txt


## Define the number of parallel groups you want to create
num_groups=100
# Count the total number of lines in the boundaries.txt file
total_lines=$(wc -l < genes.txt)
# Calculate the number of lines per group (rounded up)
lines_per_group=$(( (total_lines + num_groups - 1) / num_groups ))
# Create temporary directories for each group
for ((i=1; i<=num_groups; i++))
do
    mkdir -p "group_$i"
done
# Distribute genes into different groups
awk -v num_groups="$num_groups" -v lines_per_group="$lines_per_group" '
    { group_id = int((NR - 1) / lines_per_group) + 1; print > "group_" group_id "/subset_genes.txt" }
' genes.txt


## loop through genes
while read -r gene
do
    # Create a temp_cnv.rawcnv file for the current boundary in the group
    grep -w "$gene" $overlap_file > "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.txt"

    # Process the boundary with the R script
    Rscript association.R "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.txt" "$SLURM_ARRAY_TASK_ID" "$phenotype_file"

    # Clean up the temp_cnv.rawcnv file
    rm -f "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.txt"
done < "group_$SLURM_ARRAY_TASK_ID/subset_genes.txt"


#concatenate results into one file
cat group_*/temp_logistic.txt > association_results.txt

# clean up group directories
rm -r group_*
