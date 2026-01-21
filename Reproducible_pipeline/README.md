Scripts in this directory can be used to carry out a genome wide association test between the copy number of whole genes and a phenotype.  


All required packages are pre-installed in conda environment ```cnvGWAS.yml```.

## Step 1: Overlap
Identify if a CNV overlaps a gene.  
Input must be in bed format - one file containing cnvs and one containing genes.  
  
Run:  
```overlap.sh genes.bed cnvs.bed```

<br/> 
  
>### Example input:  
>  
>#### *genes.bed*:  
>**Columns**: chr, gene start, gene end, ensembl ID, strand, gene name.  
>| chr1 | 13182960 | 13184326 | ENSG00000179412 | -1 | HNRNPCP5 | 
>|-|-|-|-|-|-|
>  
>#### *cnvs.bed*:  
>**Columns:** chr, cnv start, cnv end, strand, sample ID.  
>| chr1 | 13176463 | 13186401 | 1 | *sample ID* | 
>|-|-|-|-|-|  
  
<br/> 
  
>### Example output:  
>**Columns:** chr, gene start, gene end, ensembl ID, strand, gene name, chr, cnv start, cnv end, strand, sample ID.
>|chr1 |13182960   |13184326 |ENSG00000179412| -1 |HNRNPCP5|chr1|13176463  |13186401|1   |  249238   |
>|-----|-----------|---------|-----------|--------|------|------|----------|--------|----|-----------|
  
<br/> 
  
## Step 2: Genome-wide association  

The following scripts carry out a genome-wide association test between gene copy number and phenotype, based on the output generated from ```overlap.sh```. For speed, the script is written to run in an array on a SLURM cluster with 100 parallel jobs.  

Run:  
```sbatch association_array.sh [intersect_genes_cnvs.txt] [phenotype_file] [covariate_file] [type]```

<br/> 

### Input
- *intersect_genes_cnvs.txt*: output file from ```overlap.sh```
- *phenotype_file*: text file containing list of sample IDs for individuals with the phenotype
- *covariate_file*: text file containing covariate data, with column names. MUST contain IID column (sample IDs).
- *type*: either ```deletion``` or ```duplication``` depending on your analysis

<br/> 

> #### Example covariate file:
> |IID |    age   |  sex   |  PC1     |PC2  |   PC3 |    PC4  |   batch |  smoke|
> |-|-|-|-|-|-|-|-|-|
> |id1 |41     | 1   |    -12.1725    |    5.39163 |-1.28103   |     0.841765   |     Batch_b001    |  never|
> |id2 |46     | 2    |   -13.0245   |     6.41514 |-0.183365   |    2.92761| Batch_b001    |  previous|

<br/> 

### Output
**Output file:** association_results.txt  
**Columns:** Gene, estimate, standard error, z-value, p-value.
