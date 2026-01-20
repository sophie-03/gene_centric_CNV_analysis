Scripts in this directory can be used to carry out a genome wide association test between the copy number of whole genes and a phenotype.  


All required packages are pre-installed in conda environment ```cnvGWAS.yml```.

## Step 1: Overlap
Identify if a CNV overlaps a gene.  
Input must be in bed format - one file containing cnvs and one containing genes.  
  
Run:  
```overlap.sh genes.bed cnvs.bed```

### Example input:  
  
#### *genes.bed*:  
**Columns**: chr, gene start, gene end, ensembl ID, strand, gene name.  
| chr1 | 13182960 | 13184326 | ENSG00000179412 | -1 | HNRNPCP5 | 
|-|-|-|-|-|-|

  
#### *cnvs.bed*:  
**Columns:** chr, cnv start, cnv end, strand, sample ID.  
| chr1 | 13176463 | 13186401 | 1 | *sample ID* | 
|-|-|-|-|-|  

  
### Example output:  
**Columns:** chr, gene start, gene end, ensembl ID, strand, gene name, chr, cnv start, cnv end, strand, sample ID.
  
|chr1 |13182960   |13184326 |ENSG00000179412| -1 |HNRNPCP5|chr1|13176463  |13186401|1   |  249238   |
|-----|-----------|---------|-----------|--------|------|------|----------|--------|----|-----------|
  
  
  
## Step 2: Genome-wide association  

The following scripts carry out a genome-wide association test between gene copy number and phenotype, based on the output generated from ```overlap.sh```. For speed, the script is written to run in an array on a SLURM cluster.  

Run:  
