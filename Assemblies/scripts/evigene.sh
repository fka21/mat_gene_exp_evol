#!/bin/bash

# Example commands used to eliminate redundancy
# All steps were repeated for each de novo generated transcriptome
# The steps are following the suggestions from the developer of EviGene 

#Formatting header of the concatanated file with all the generated transcripts from all softwares and parameters
/home/s9/tsai/work/Ferenc/Tools/evigene/scripts/rnaseq/trformat.pl \
-prefix=Hduj \
-input /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/Hduj_concat_all.fa > /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/Hduj_trformat.fasta

#Running evigene
mkdir /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/evigene
cd /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/evigene
/home/s9/tsai/work/Ferenc/Tools/evigene/scripts/prot/tr2aacds4.pl -mrnaseq /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/Hduj_trformat.fasta -NCPU 12 -MAXMEM 512000 -logfile Hduj_evigene.log

#Preparing uniprot database for annotation
#In the directory with the assemblies
diamond makedb --in uniprot_sprot.fasta --db uniprot_sprot.fasta

#Blasting against uniprot-sp
env aaset=okayset/Hduj_trformat.okay.aa refaa=/sysdev/s9/tsai/Ferenc/Datasets/H_dujardini/Assembly/evigene/uniprot_sprot.fasta ncpu=4 datad=/home/s9/tsai/work/Ferenc/Datasets/H_dujardini/Assembly/evigene/temp /home/s9/tsai/work/Ferenc/Tools/evigene/scripts/prot/run_evgaablast.sh

#Crosreference table
grep "^>" uniprot_sprot.fasta | awk -F 'OS=' '{print $1}' | perl -pe 's/ /\t/'  | sed "s/^>//g" > crosref_uniprot.tbl

#Naming genes
/home/s9/tsai/work/Ferenc/Tools/evigene/scripts/prot/namegenes.pl -blasttab refaa-Hduj_trformat.okay.aa.btall -names crosref_uniprot.tbl

#Final cleanup
~/work/Ferenc/Tools/evigene/scripts/evgmrna2tsa2.pl \
-trclass /sysdev/s9/tsai/Ferenc/Datasets/H_dujardini/Assembly/Hduj_trformat.trclass \
-mrna /sysdev/s9/tsai/Ferenc/Datasets/H_dujardini/Assembly/evigene/okayset/Hduj_trformat.okay.tr \
-idprefix Hduj \
-organism=Hypsibius_dujardini \
-NCPU=8 \
-novectrim \
-names /sysdev/s9/tsai/Ferenc/Datasets/H_dujardini/Assembly/evigene/refaa-Hduj_trformat.okay.aa.btall.named

