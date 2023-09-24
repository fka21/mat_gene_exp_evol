#!/bin/bash

# Script used for trimming raw reads
# Handles single end reads and paired end reads also

# Using conda environment
conda activate quant

# Looping through files paired end files
for k in /home/s9/tsai/work/Ferenc/Datasets/*
do

cd "$k"

DIR=$(echo ${k})

 for i in ${DIR}/*_1.fastq.gz
 do
   # Getting basename of file
   SAMPLE=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/_1\.fastq\.gz//")
   
   echo ${SAMPLE}   
    
   #Trimming
   fastp \
   -i ${DIR}/${SAMPLE}_1.fastq.gz \						# Using basename and extending for PE reads
   -I ${DIR}/${SAMPLE}_2.fastq.gz \
   -o ${DIR}/${SAMPLE}_1_fastp.fastq.gz \
   -O ${DIR}/${SAMPLE}_2_fastp.fastq.gz \
   -q 25 \									# Parameters were iteratively determined using FastQC and universally used
   -5 \
   --cut_front_window_size 5 \
   --cut_front_mean_quality 25 \
   -r \
   --cut_right_window_size 5 \
   --cut_right_mean_quality 25 \
   --correction \
   --low_complexity_filter \
   --detect_adapter_for_pe \
   -w 7 \
   -l 35
 done


# Looping through single end files
  for i in ${DIR}/*.fastq
  do
   # Getting basename of file
   SAMPLE=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/\.fastq//")

   echo ${SAMPLE}

   fastp \
   -i ${DIR}/${SAMPLE}.fastq \
   -o ${DIR}/${SAMPLE}_RNA-Seq.fastp.fastq \
   -q 25 \
   -5 \
   --cut_front_window_size 5 \
   --cut_front_mean_quality 25 \
   -r \
   --cut_right_window_size 5 \
   --cut_right_mean_quality 25 \
   --low_complexity_filter \
   -w 16 \
   -l 35
 done

done
