Scripts in this directory can be used to carry out a genome wide association test between the copy number of whole genes and a phenotype.  

## Step 1: Overlap
Identify if a CNV overlaps a gene (+ flanking region).  
Input must be in bed format - one file containing cnvs and one containing genes.  
  
Run:  
```overlap.sh genes.bed cnvs.bed```

### Ouput format:  
