# Expression analysis

## Step 1: Call CNVs from GTEX WGS data

### 1.1 Download the GTEX data

Download data from GTEx using ```download_gtex-data.sh```

### 1.2 Install CNVpytor

```
source activate GTEX #activate premade conda environment: GTEX
conda install -c bioconda cnvpytor #install CNVpytor
```

### 1.3 Call CNVs using CNVpytor

The ```1by1_gtex_cnvpytor.sh``` script will download a CRAM file from the manifest, call CNVs using CNVpytor, then delete the CRAM file and downoad the next.


### 1.4 Format file names

```remove.sh``` just removes the ```'``` and ```"``` from the call file names 


## Step 2: Identify which CNVs span whole genes

### 2.1 Format files and filter CNV calls
 Both the list of GRCh38.p13 genes and the CNV call files need to be formatted to BED format.

```genelist_format.R```  filters the gtf file downloaded from [Gencode](https://www.gencodegenes.org/human/) (content: comprehensive gene annotation, regions: chr, gtf file) to contain genes only (excludes exons, stop/start codons etc.). It then formats the file to be in BED format, with the 4th column containing the ensembl gene ID and removes the excess columns.  
It also outputs a bed file that has added 1kb to the start and end positions of each gene, this the bed file that will be used in analysis. The added buffer aims to capture promoter regions for each gene.

```format_cnvpytor_calls.R```  reformats the output files to BED format and filters the CNV calls. Filtering is done based on q0 (fraction of reads mapped with zero quality within call region). 


### 2.2 Find the CNVs that span whole genes

Bedtools intersect can only handle ~400 input files, so the CNV call files need to be divided into subsets (3 subsets). ```bedtools_intersect.sh``` can then be run, which runs bedtools intersect on the 3 subsets and concatenates the outputs into one file.

### 2.3 Correlation analysis

Carried out in R: ```correlation_analysis.R```

