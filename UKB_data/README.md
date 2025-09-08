## 1. Data Wrangling

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
