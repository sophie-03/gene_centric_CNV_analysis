## 1. Data prep

Using CNVs called previously from UKB.
Data needs to be formatted in the style of PennCMV output, which is done in the script ```format_geneCN.R```.  

Create covariate file: ```covariates.R```

## 2. Association tests: Logistic regression

For each binary phenotype, I run a logistic regression to see if deletion or duplication status influences phenotype. Separate tests are run for duplications and deletions. The ```logistic_dups_array_cloud.sh``` or ```logistic_dels_array_cloud.sh``` scripts runs an array job that seperates genes into 100 groups, and then loops over each gene in the group, running a logistic regression test. Outputs are then concatenated into an output file. The phenoytpe is specified in the command line.
