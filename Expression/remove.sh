#!/bin/bash

cd ~/GTEX/CNVpytor_calls

for file in *GTEX*; do
var=$(echo $file | sed 's/\"//g')
#echo $var
mv $file $var
done
