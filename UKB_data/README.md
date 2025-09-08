## 1. Calling CNVs

### 1.1 Making signal intensity files
Formatting UKB data for PennCNV:
UKB data gives you one logRratio file and one BAF file. Each file contains the data for every individual. However, the input for PennCNV requires signal intensity files for each individual. Signal intensity files are text files that contain both log and BAF information for each marker. The header for the signal intensity files are: ```Name | Chromosome | Position | LogRratio | BAF```

Therefore, the files from UKB need to be split into individuals (columns) and pasted into a file with the corresponding column from the other file, along with marker information. Due to the number of individuals, this is a time-consuming process. And so I first cut the large data files into chunks of 10,000 individuals, and then split into individual files in parallel.

Splitting the large logRratio and BAF files into smaller files of 10,000 individuals is done using the ```make_chunks.sh``` script
```
./make_chunks.sh #just change to the chr that you want to download
```

Now we've broken up the huge files, we can make the individual files in parallel. The ```paralell_chunks.sh``` file contins the parallel command to run 5 R scripts simultaeneously (one R script for each chunk file i.e ```indv_files-1-10.R``` etc). The R script creates signal intensity files for each individual, containing the following columns: ```marker_name | chr | position | logRration | BAF```

```
./parallel_chunks.sh
```

Now check that you have 488,377 inidividual signal intensity files, and then we can delete the 10,000 chunk files
```
rm -r ~/UKB_data/chr2/10000_chunks/ukb*
```

### 1.2 PennCNv
```PennCNV.sh``` contains the parallel command that runs each line of ```PennCNV_parallel.sh``` in parallel.

>#### Chromosome X
>Calling CNVs from chromosome X requires you to specify -chrx when calling the CNvs and to input a sexfile. Sex information is found in the fam file, which is downloaded via gfetch, the same as downloading the data, however uses the ```-m``` flag. It doesn't matter which chr you use to download the fam file, as they are all the same.
>```
>./gfetch 22431 -c22 -m
>```
>The sexfile is created from the ```create_sexfile.R``` script, which takes the sex information found in the fam file. Sex information is found in the 5th column, where 1 = male, 2 = female and 0 = unknown. The Rscript creates a sexfile containing the sex of the individual and the sample file it corresponds to (these sample files are numbered based on the order they appear in the log and BAF files, and are NOT the UKB IDs).
>```
>Rscript create_sexfile.R
>```
>The UKB data labels chrX as chr23, and so the name of the chromosomes need to be changed in each sample file (It would probably be more efficient to change this before the individual files were all made, but I didn't realise it was a problem until now). The following script does this:
>```
>./format_X_samples.sh
>```
>PennCNV can now be ran on chrX. The ```-chrx``` and ```--sexfile``` flags need to be used, so specific scripts have been made for the X chromosome (```PennCNV_X.sh``` and ```PennCNV_parallel_X.sh```)
>```
>./PennCNV_X.sh
>```

## 2. Identify CNVs that span whole genes

Before and further analysis can be done, the CNV calls from PennCNV need to be formatted into a BED file.  
  
Firstly, the individuals need to be mapped to UKB IDs, so that we can later compare phenotypes to these CNV calls. The order of individuals in the log and BAF files from UKB correspond to the order of IDs in the fam file (downloaded earlier). The ```create_mapfile.R``` script creates a file that maps the UKB ID to the sample number from the log/BAF fies (that have been used so far).  
  
Then I need to make the cnv calls into a bed file, so I can use bedtools. ```cnvs_to_bed.R``` does this. The resulting bed file also contains the sample number used thus far, and the UKB ID.

### 2.2 Get gene list

``` genes_to_bed.R``` converts gene list to bed file
  
 ```buffer_to_genelist.R``` adds a 1kb buffer to both the start and end points of the genes. The output of this script is a bed file containing a list of genes that can be used in bedtools.  
  
  
### 2.3 Identify CNVs that span whole genes  
  

```bedtools.sh``` identifies CNVs that span whole genes 
      
## 3: Tumour suppressor genes and oncogenes  
  
### 3.1 Make TSG and oncogene bed files

A list of cancer genes was obtained from the [cosmic gene census](https://cancer.sanger.ac.uk/census). This list was then subset into a list containing only the tumour suppressor genes and a list containing only oncogenes, and both turned into a bed file: ```oncogenes.bed``` and ```TSGs.bed```. This is done by the script ```Oncogenes_to_bed.R``` and ```TSGs_to_bed_files.R```.  

### 3.2 Analysis

Analysis carried out in TSG_ONC_analysis.qmd
