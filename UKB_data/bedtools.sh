#!/bin/bash

## cloud ##

cd ~/UKB/genes_in_cnvs

# get genes that are encompassed completely by CNVs
bedtools intersect -wa -wb -a biomart_grch37_genes.bed -b ../cnv_calls.bed -filenames -f 1.0 > intersect_genes_cnvs.txt

# get genes where >5% of the gene is disrupted by a CNV
bedtools intersect -wa -wb -a biomart_grch37_genes.bed -b ../cnv_calls.bed -filenames -f 0.05 >5_intersect_genes_cnvs.txt
