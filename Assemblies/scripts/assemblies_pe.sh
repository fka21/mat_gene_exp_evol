#!/bin/bash

#Call format: nohup bash -i assemblies.sh

#ASSEMBLIES
# This script was used on paired end data which already contained processed reads (see Quantification/scripts/fastp.sh) 
# Looping through all reads, assemblies of stage/replicate specific transcripts
# Using an overassembly approach (multiple assemblers with multiple kmers) followed by a redundancy elimination step (EviGene)

for i in /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/*_1_fastp.fastq.gz
do

DIR="/home/s9/tsai/work/Ferenc/Datasets/H_dujardini/"
READ=$(echo ${i} | awk -F "/" '{print $9}' | sed "s/_1_fastp\.fastq\.gz//g") 

conda activate transabyss

transabyss --pe \
${DIR}/${READ}_1_fastp.fastq.gz \
${DIR}/${READ}_2_fastp.fastq.gz \
--outdir /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/transabyss --name ${READ}_k40.fasta -k 40 --threads 6

conda activate trinity

Trinity \
--seqType fq \
--max_memory 500G \
--left ${DIR}/${READ}}_1_fastp.fastq.gz \
--right ${DIR}/${READ}_2_fastp.fastq.gz \
--full_cleanup \
--CPU 6 \
--output /home/s9/tsai/work/Ferenc/Datasets/T_transversa/Assembly/trinity/${READ}_trinity


/home/s9/tsai/work/Ferenc/Tools/SPAdes-3.14.0-Linux/bin/rnaspades.py \
-1 ${DIR}/${READ}_1_fastp.fastq.gz \
-2 ${DIR}/${READ}_2_fastp.fastq.gz \
--threads 15 \
--only-assembler \
-k 40,56,64 \
-o /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/rnaspades/${READ}


done

