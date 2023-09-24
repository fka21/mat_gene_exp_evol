#!/bin/bash

conda activate quant

#What are the length read lengths distributions after processing them
#This information will be used for transcriptome indexing

for k in /home/s9/tsai/work/Ferenc/Datasets/*
do

cd "$k"

DIR=$(echo ${k})

 for i in ${DIR}/*.fastq.gz
 do
   SAMPLE=$(echo ${i} | awk -F "/" '{print $11}' | sed "s/.fastq\.gz//")

   LENGTH=$(zcat $i | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}')

   echo ${DIR}/${SAMPLE} ${LENGTH}

 done
done
