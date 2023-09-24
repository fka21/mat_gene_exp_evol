#!/bin/bash

conda activate quant

# Creating decoy sequences as suggested by developers of salmon
for k in /home/s9/tsai/work/Ferenc/Datasets/*
do

cd "$k"

DIR=$(echo ${k})
GTF=$(ls ${k}/*gtf.gz)
GFF=$(ls ${k}/*gff3.gz)
GENOME=$(ls ${k}/*toplevel.fa.gz)
CDNA=$(ls ${k}/*cdna.all.fa.gz)

echo $DIR
 
grep "^>" <(gunzip -c $GENOME ) | cut -d " " -f1 > ${DIR}/decoys.txt
sed -i.bak -e 's/>//g' ${DIR}/decoys.txt
cat $CDNA $GENOME > ${DIR}/gentrome.fa.gz

done

# Indexing each dataset with appropriate k-mer sizes to the read distributions
# Kmer sizes were determined using read_lengths.sh
#  Approx 0.6*read_length was used

salmon index -t /home/s9/tsai/work/Ferenc/Datasets/A_stephensi/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/A_stephensi/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/A_stephensi/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/A_suum/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/A_suum/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/A_suum/salmon_index -k19
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/B_germanica/Bger_cdna.all.fa.gz -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/B_germanica/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/B_lancelatum_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/B_lancelatum_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/B_lancelatum_2018/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/B_mori/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/B_mori/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/B_mori/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/B_taurus/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/B_taurus/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/B_taurus/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/C_angaria_2017/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/C_angaria_2017/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/C_angaria_2017/salmon_index -k19
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/C_elegans_2017/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/C_elegans_2017/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/C_elegans_2017/salmon_index -k19
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/C_gigas_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/C_gigas_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/C_gigas_2018/salmon_index -k23
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/C_hircus/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/C_hircus/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/C_hircus/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/C_intestinalis/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/C_intestinalis/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/C_intestinalis/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_ananassae_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_ananassae_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_ananassae_2018/salmon_index -k23
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_erecta_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_erecta_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_erecta_2018/salmon_index -k23
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_melanogaster_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_melanogaster_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_melanogaster_2018/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_mojavensis_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_mojavensis_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_mojavensis_2018/salmon_index -k23
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_persimilis_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_persimilis_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_persimilis_2018/salmon_index -k23
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_rerio_2017/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_rerio_2017/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_rerio_2017/salmon_index -k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_simulans_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_simulans_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_simulans_2018/salmon_index -k23
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_virilis_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_virilis_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_virilis_2018/salmon_index -k23
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_willistoni_2018/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_willistoni_2018/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_willistoni_2018/salmon_index -k23
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/D_yakuba_2015/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/D_yakuba_2015/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/D_yakuba_2015/salmon_index	-k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/G_gallus/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/G_gallus/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/G_gallus/salmon_index	-k31
salmon index -t	/home/s9/tsai/work/Ferenc/Datasets/H_dujardini/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/H_dujardini/salmon_index	-k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/H_sapiens/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/H_sapiens/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/H_sapiens/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/M_capitata/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/M_capitata/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/M_capitata/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/M_franciscanus/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/M_franciscanus/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/M_franciscanus/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/M_leidyi/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/M_leidyi/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/M_leidyi/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/M_mulatta/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/M_mulatta/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/M_mulatta/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/M_musculus/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/M_musculus/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/M_musculus/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/N_vectensis/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/N_vectensis/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/N_vectensis/salmon_index -k23
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/P_caudatus/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/P_caudatus/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/P_caudatus/salmon_index -k23
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/P_dumerili/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/P_dumerili/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/P_dumerili/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/P_miniata/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/P_miniata/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/P_miniata/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/S_carpocapsae_2017/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/S_carpocapsae_2017/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/S_carpocapsae_2017/salmon_index -k19
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/S_feltiae_2017/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/S_feltiae_2017/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/S_feltiae_2017/salmon_index -k19
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/S_scrofa/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/S_scrofa/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/S_scrofa/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/T_castenum_2019/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/T_castenum_2019/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/T_castenum_2019/salmon_index -k31
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/T_transversa/Ttra_cdna.all.fa -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/T_transversa/salmon_index -k23
salmon index -t /home/s9/tsai/work/Ferenc/Datasets/X_laevis_2016/gentrome.fa.gz -d /home/s9/tsai/work/Ferenc/Datasets/X_laevis_2016/decoys.txt -p 15 -i /home/s9/tsai/work/Ferenc/Datasets/X_laevis_2016/salmon_index -k31

#Manually did these
#salmon index -t Lvar_cdna.all.fa -p 8 -i ./salmon_index -k 31
#salmon index -t Htub_cdna.all.fa -p 8 -i ./salmon_index -k 31
#
