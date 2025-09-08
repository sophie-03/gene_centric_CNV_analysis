#!/bin/bash

#this script downloads CRAM files from GTEX 1 by 1 and calls CNVs using CNVpytor. It then deletes the CRAm file and downloads the next.

#for i in {0..1}
for i in 1
do
jq .[$i] ../GTEX/manifest_subset.json > download_manifest.json #extract individual objects and add to new manifest file for download
#format download file
sed -i '1i [' download_manifest.json #add square bracket to start of download file
echo "]" >> download_manifest.json #add square bracket to end of  download file

#download individual CRAM file
../bin/gen3/gen3-client configure --profile=anvil --cred=credentials.json --apiendpoint=https://gen3.theanvil.io #configure
yes | ../bin/gen3/gen3-client download-multiple --profile=anvil --manifest=download_manifest.json --download-path=Crams --protocol=s3 #download

#call CNVs using cnvpytor
samtools index Crams/GTEX-* #index CRAM file
cnvpytor -root GTEX.pytor -rd Crams/*.cram -T Homo_sapiens_assembly38.fasta -j 1 #import read depth signal
cnvpytor -root GTEX.pytor -his 1000 #calculate read depth histograms, GC correction and statistics type 9bin size =1000)
cnvpytor -root GTEX.pytor -partition 1000 #partitioning using mean-shift method
cnvpytor -root GTEX.pytor -call 1000 > CNVpytor_calls/GTEX.calls.1000.tsv #call CNVs

#rename call file so it matches sample ID
cd CNVpytor_calls
old_name=$"GTEX.calls.1000.tsv"
new_name=$(jq '.[0] .file_name' ../download_manifest.json).calls.1000.tsv
mv "$old_name" "$new_name"

#delete CRAMs and .pytor file
cd ../
rm -r Crams/GTEX*
rm -r GTEX.pytor
done
