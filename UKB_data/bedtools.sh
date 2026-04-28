#!/bin/bash

## cloud ##

cd ~/UKB/genes_in_cnvs

# get genes that are encompassed completely by CNVs
bedtools intersect -wa -wb -a biomart_grch37_genes.bed -b ../cnv_calls.bed -filenames -f 1.0 > intersect_genes_cnvs.txt


