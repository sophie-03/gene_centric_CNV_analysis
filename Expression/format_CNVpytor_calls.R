#this script formats the CNVpytor output into bed format and filters CNVs by mapping quality (>50%)

setwd("~/GTEX/CNVpytor_calls")

library(data.table)
library(tidyr)
library(dplyr)

#list of file names to loop through
fileNames <- Sys.glob("GTEX*.tsv")


for (file in fileNames) {

  #read in file from list of files in directory
  calls <- fread(file)
  calls <- as_tibble(calls)

  #filter calls based on associated mapping quality(q0)
  calls <- calls[(calls['V9']<=0.5),]

  #keep just the columns we need
  calls <- calls %>% select(V1, V2, V4)

  #reorder columns
  calls <- calls[,c(2,1,3)]

  #separate out chr, start and end locations
  calls <- separate(calls, V2, into = c("chr", "cnvStart", "cnvEnd"), sep = "[-,:]")

  #keep only chr, remove unplaced scaffolds
    calls <- calls[calls$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                          "chr8", "chr9", "chr10", "chr11", "chr12",
                          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                          "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),]

    #change end positions so they are in line with chr sizes
    #calls$cnvEnd <- ifelse((calls$chr=='chr1')&(calls$cnvEnd>248956422), 248956422, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr2')&(calls$cnvEnd>242193529), 242193529, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr3')&(calls$cnvEnd>198295559), 198295559, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr4')&(calls$cnvEnd>190214555), 190214555, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr5')&(calls$cnvEnd>181538259), 181538259, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr6')&(calls$cnvEnd>170805979), 170805979, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr7')&(calls$cnvEnd>159345973), 159345973, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr8')&(calls$cnvEnd>145138636), 145138636, calls$cnvEnd)
  #calls$cnvEnd <- ifelse((calls$chr=='chr9')&(calls$cnvEnd>138394717), 138394717, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr10')&(calls$cnvEnd>133797422), 133797422, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr11')&(calls$cnvEnd>135086622), 135086622, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr12')&(calls$cnvEnd>133275309), 133275309, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr13')&(calls$cnvEnd>114364328), 114364328, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr14')&(calls$cnvEnd>107043718), 107043718, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr15')&(calls$cnvEnd>101991189), 101991189, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr16')&(calls$cnvEnd>90338345), 90338345, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr17')&(calls$cnvEnd>83257441), 83257441, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr18')&(calls$cnvEnd>80373285), 80373285, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr19')&(calls$cnvEnd>58617616), 58617616, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr20')&(calls$cnvEnd>64444167), 64444167, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr21')&(calls$cnvEnd>46709983), 46709983, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chr22')&(calls$cnvEnd>50818468), 50818468, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chrX')&(calls$cnvEnd>156040895), 156040895, calls$cnvEnd)
    #calls$cnvEnd <- ifelse((calls$chr=='chrY')&(calls$cnvEnd>57227415), 57227415, calls$cnvEnd)

  #write output
  write.table(calls, paste(file,".bed",sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep = "\t")

  #clear memory
  rm(file)
}
