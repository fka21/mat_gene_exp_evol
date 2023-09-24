#!/bin/bash

# Script used for assessing the completeness of transcriptomes assembled
for i in ./*/*cdna.all.fa
do

OUT=$(echo ${i} | awk -F "/" '{print $2}')

busco -i $i -l /home/s9/tsai/work/Ferenc/Tools/BUSCO/busco_downloads/lineages/metazoa_odb10 -o BUSCO --out_path ${OUT}/Assembly -m transcriptome

done
