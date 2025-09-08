setwd("~/UKB")

library(tidyr)
library(dplyr)
library(data.table)

#import data
calls <- read.table("all_calls.rawcnv")
calls <- as_tibble(calls)

## format data into bed format
#separate chr, cnv start and cnv end into columns
calls <- separate(calls, V1, into =c("chr","cnvStart","cnvEnd"), sep="[:,-]")
#separate copy number into column
calls <- separate(calls,V4, into = c(NA,"cn"), sep="[=]")
#separate sample number
calls <- separate(calls,V5, into = c(NA,"sample"), sep="sample")
calls <- separate(calls,sample, into = c("sample"), sep=".t")
calls

#keep only the columns I need
calls <- calls %>% select(chr, cnvStart, cnvEnd, cn, sample)
calls

##add column for UKB ID
#read in mapping file
map <- read.table("map_file.txt", header = TRUE)
#map UKB ids to samples using the sample column
calls <- merge(calls, map, by = "sample")

#move sample column to end
calls <- calls[, c("chr","cnvStart","cnvEnd","cn","sample","UKB_id")]
#order by chr number
calls <- calls[order(calls$chr),]

write.table(calls, "cnv_calls.bed", quote=FALSE, col.names=FALSE, row.names=FALSE, sep ="\t")
