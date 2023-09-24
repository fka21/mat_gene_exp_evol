#!/bin/bash

#ASSEMBLIES
# This script was used on single end data which already contained processed reads (see Quantification/scripts/fastp.sh) 
# Looping through all reads, assemblies of stage/replicate specific transcripts
# Using an overassembly approach (multiple assemblers with multiple kmers) followed by a redundancy elimination step (EviGene)

for i in /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/Reads/single/*
do

gzip -d "${i}"
conda activate transabyss

READ=$(echo ${i} | awk -F "/" '{print $12}' | sed "s/\.gz//g")
DIR="/home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/Reads/single/"


transabyss \
--se ${DIR}/${READ} \
--outdir /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/transabyss \
--name ${i}_k27.fasta \
-k 27 \
--threads 10

conda activate trinity

Trinity \
--seqType fq \
--max_memory 500G \
--single ${DIR}/${READ} \
--full_cleanup \
--CPU 10 \
--output /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/trinity/${i}_trinity

/home/s9/tsai/work/Ferenc/Tools/SPAdes-3.14.0-Linux/bin/rnaspades.py \
-s ${DIR}/${READ} \
--threads 10 \
--only-assembler \
-k 21,27 \
-o /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/rnaspades/${i}_rnaspades

gzip "${i}"

done

