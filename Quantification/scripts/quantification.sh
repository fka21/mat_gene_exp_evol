#!/bin/bash

# Script used to pseudoalign the trimmed reads to the indexed transcriptomes
# Indexing was performed previously
# Handles single end reads and paired end reads 

# Using conda environment
conda activate salmon

# Looping through subdirectories
for k in /home/s9/tsai/work/Ferenc/Datasets/*
do

cd "$k"

# Saving filenames to variables
DIR=$(echo ${k})
GTF=$(ls ${k}/*gtf.gz)
GFF=$(ls ${k}/*gff3.gz)
GENOME=$(ls ${k}/*toplevel.fa.gz)
CDNA=$(ls ${k}/*cdna.all.fa.gz)

echo $DIR

   # Looping through files
   for i in ${DIR}/*_1_fastp.fastq.gz
   do
     
     # Getting basenames
     SAMPLE1=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/_1_\.fastp\.fastq\.gz//")
     SAMPLE2=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/_1_fastp\.fastq\.gz/_2_fastp.fastq.gz/")
     
     # Quantification of PE reads
     # Parameters are set following the suggestions from the developers
     salmon quant -i ${DIR}/salmon_index -l A -1 ${DIR}/${SAMPLE1} -2 ${DIR}/${SAMPLE2} -o ${DIR}/${SAMPLE1}_salmon --validateMappings -p 10 --seqBias --gcBias --posBias --numBootstraps 100
   done

   for i in ${DIR}/*_RNA-Seq.fastp.fastq.gz
   do
     # Looping through files
     SAMPLE=$(echo ${i} | awk -F "/" '{print $9}')
     echo $SAMPLE
     
     # Quantification of PE reads
     # Parameters are set following the suggestions from the developers
     salmon quant -i ${DIR}/salmon_index -l A -r ${DIR}/${SAMPLE} -o ${DIR}/${SAMPLE}_salmon -p 2 --seqBias --posBias --validateMappings --numBootstraps 100
   done

 done


done

