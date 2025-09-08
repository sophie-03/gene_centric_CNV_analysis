#!/bin/bash

cd ~/GTEX/bed_files/

## have to run bedtools on subsets of samples as it can only take ~400 files as an input at a time
#for subset1
time bedtools intersect -wa -wb -a ../hg38_genelist_buffer.bed -b subset1/*.bed -filenames -f 1.0 > subset1/subset1_bedtools_intersect.out
#for subset2
time bedtools intersect -wa -wb -a ../hg38_genelist_buffer.bed -b subset2/*.bed -filenames -f 1.0 > subset2/subset2_bedtools_intersect.out
#for subset3
time bedtools intersect -wa -wb -a ../hg38_genelist_buffer.bed -b subset3/*.bed -filenames -f 1.0 > subset3/subset3_bedtools_intersect.out

#combine all outputs into one file
cat subset1/subset1_bedtools_intersect.out subset2/subset2_bedtools_intersect.out subset3/subset3_bedtools_intersect.out > all_samples_buffer_intersect.out
