#!/bin/bash

echo "Gene bed file: $1"
echo "CNV bed file: $2"

# complete overlap between genes and cnvs
bedtools intersect -wa -wb -a "$1" -b "$2" -filenames -f 1.0 > intersect_genes_cnvs.txt
