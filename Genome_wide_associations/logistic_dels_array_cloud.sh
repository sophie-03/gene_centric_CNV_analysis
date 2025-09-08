#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --mail-type=ALL # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --array=1-100
#SBATCH --nice=10000

source ~/.bashrc
conda activate phewas

# Get the command line argument for the PHENOTYPE variable
PHENOTYPE=$1

# Set the working directory dynamically based on the PHENOTYPE variable
WORKING_DIR="/define/home/smatthews/pheWAS/binary/pheno_${PHENOTYPE}"
cd "$WORKING_DIR" || exit 1

#create boundaries file (only needs to be created once)
# awk '{print $1}' ../gene_cnvs/any_overlap/protein_coding_geneCN.rawcnv | sort -u > ../gene_cnvs/any_overlap/protein_coding_boundaries.txt

# Define the number of parallel groups you want to create (should be 5 in this case)
num_groups=100
# Count the total number of lines in the boundaries.txt file
total_lines=$(wc -l < /define/home/smatthews/pheWAS/gene_cnvs/protein_coding_boundaries.txt)
# Calculate the number of lines per group (rounded up)
lines_per_group=$(( (total_lines + num_groups - 1) / num_groups ))
# Create temporary directories for each group
for ((i=1; i<=num_groups; i++))
do
    mkdir -p "group_$i"
done
# Distribute boundaries into different groups
awk -v num_groups="$num_groups" -v lines_per_group="$lines_per_group" '
    { group_id = int((NR - 1) / lines_per_group) + 1; print > "group_" group_id "/protein_coding_boundaries.txt" }
' /define/home/smatthews/pheWAS/gene_cnvs/protein_coding_boundaries.txt

while read -r boundary
do
    # Create a temp_cnv.rawcnv file for the current boundary in the group
    grep "$boundary" /define/home/smatthews/pheWAS/gene_cnvs/geneCN.rawcnv > "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.rawcnv"

    # Process the boundary with the R script
    Rscript /define/home/smatthews/pheWAS/scripts/logistic_dels_cloud.R "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.rawcnv" "$SLURM_ARRAY_TASK_ID" "${PHENOTYPE}"

    # Clean up the temp_cnv.rawcnv file
    rm -f "group_$SLURM_ARRAY_TASK_ID/temp_cnvs.rawcnv"
done < "group_$SLURM_ARRAY_TASK_ID/protein_coding_boundaries.txt"

#concatenate results into one file
cat group_*/temp_logistic.txt > "${PHENOTYPE}"_dels_logistic_output.txt

# clean up group directories
#rm -r group_*
