#Script for building data frame from orthogroup data with expression data
#Check for orthogroup similarities across species
#Ferenc Kagan
#19.12.2020


##### LIBRARIES #####
library(adephylo)
library(phylobase)
library(ape)
library(phylotools)
library(phytools)
library(geiger)
library(pheatmap)
library(tidyverse)
library(tximport)
library(edgeR)
library(motmot)
library(castor)
library(MuMIn)
library(tidyverse)
library(geiger)
library(zFPKM)

##### READ IN DATA ######
maternal <- readRDS( "~/Documents/Gene_expr_evol/Intermediate_files/downregulated_IDs_NFE-T.RDS")
maternal_tpm <- readRDS("~/Documents/Gene_expr_evol/Intermediate_files/maternal_IDs.RDS")
ogroups <- read_tsv("~/Documents/Orthology/N0.tsv", na = c("", NA)) %>%
  select(-c(HOG,`Gene Tree Parent Clade`))


##### FUNCTIONS #####
source("~/Documents/Gene_expr_evol/Scripts/functions.R")

##### SAVE & LOAD #####
#load("/Users/ferenckagan/Documents/Bioinformatic_analysis/laptop/Bioinformatic_analysis/Quantification/Intermediate_files/OG_tbl.RData")
#save.image('/Users/ferenckagan/Documents/Bioinformatic_analysis/laptop/Bioinformatic_analysis/Quantification/Intermediate_files/OG_tbl.RData')

##### LOAD SALMON COUNTS #####

variables <- list.files("~/Documents/Gene_expr_evol/Intermediate_files/", pattern = "txi_.....RDS", full.names = T)
variable_names <- str_extract_all(variables, "txi_....", simplify = T)

for(i in 1:length(variables)){
  temp <- readRDS(variables[i])$abundance
  
  #Across species normalization with TPM10K approach
  temp <- (temp * dim(temp)[1]) / 10^4
  
  temp_name <- variable_names[i]
  temp_name <- str_replace_all(temp_name, "txi", "tpm")
  assign(temp_name, temp)
}

#Some datasets gene names need adjustments to combine with OrthoFinder results
rownames(tpm_sfel) <- str_remove_all(rownames(tpm_sfel), "gene:")

##### DATA DESCRIPTION #####
coldata_aste <- data.frame(row.names=colnames(tpm_aste), Stages_aste = factor(str_remove_all(colnames(tpm_aste), "_[0-9]+")))
coldata_asum <- data.frame(row.names=colnames(tpm_asum), Stages_asum = factor(str_remove_all(colnames(tpm_asum), "_[0-9]+")))
coldata_bger <- data.frame(row.names=colnames(tpm_bger), Stages_bger = factor(str_remove_all(colnames(tpm_bger), "_[0-9]+")))
coldata_btau <- data.frame(row.names=colnames(tpm_btau), Stages_btau = factor(str_remove_all(colnames(tpm_btau), "_[0-9]+")))
coldata_chir <- data.frame(row.names=colnames(tpm_chir), Stages_chir = factor(str_remove_all(colnames(tpm_chir), "_[0-9]+")))
coldata_blan <- data.frame(row.names=colnames(tpm_blan), Stages_blan = factor(str_remove_all(colnames(tpm_blan), "_[0-9]+")))
coldata_bmor <- data.frame(row.names=colnames(tpm_bmor), Stages_bmor = factor(str_remove_all(colnames(tpm_bmor), "_[0-9]+")))
coldata_cana <- data.frame(row.names=colnames(tpm_cana), Stages_cana = factor(str_remove_all(colnames(tpm_cana), "_[0-9]+")))
coldata_cele <- data.frame(row.names=colnames(tpm_cele), Stages_cele = factor(str_remove_all(colnames(tpm_cele), "_[0-9]+")))
coldata_cgig <- data.frame(row.names=colnames(tpm_cgig), Stages_cgig = factor(str_remove_all(colnames(tpm_cgig), "_[0-9]+")))
coldata_dana <- data.frame(row.names=colnames(tpm_dana), Stages_dana = factor(str_remove_all(colnames(tpm_dana), "_[0-9]+")))
coldata_dere <- data.frame(row.names=colnames(tpm_dere), Stages_dere = factor(str_remove_all(colnames(tpm_dere), "_[0-9]+")))
coldata_dmel <- data.frame(row.names=colnames(tpm_dmel), Stages_dmel = factor(str_remove_all(colnames(tpm_dmel), "_[0-9]+")))
coldata_dmoj <- data.frame(row.names=colnames(tpm_dmoj), Stages_dmoj = factor(str_remove_all(colnames(tpm_dmoj), "_[0-9]+")))
coldata_dper <- data.frame(row.names=colnames(tpm_dper), Stages_dper = factor(str_remove_all(colnames(tpm_dper), "_[0-9]+")))
coldata_drer <- data.frame(row.names=colnames(tpm_drer), Stages_drer = factor(str_remove_all(colnames(tpm_drer), "_[0-9]+")))
coldata_dsim <- data.frame(row.names=colnames(tpm_dsim), Stages_dsim = factor(str_remove_all(colnames(tpm_dsim), "_[0-9]+")))
coldata_dvir <- data.frame(row.names=colnames(tpm_dvir), Stages_dvir = factor(str_remove_all(colnames(tpm_dvir), "_[0-9]+")))
coldata_dwil <- data.frame(row.names=colnames(tpm_dwil), Stages_dwil = factor(str_remove_all(colnames(tpm_dwil), "_[0-9]+")))
coldata_dyak <- data.frame(row.names=colnames(tpm_dyak), Stages_dyak = factor(str_remove_all(colnames(tpm_dyak), "_[0-9]+")))
coldata_ggal <- data.frame(row.names=colnames(tpm_ggal), Stages_ggal = factor(str_remove_all(colnames(tpm_ggal), "_[0-9]+")))
coldata_hsap <- data.frame(row.names=colnames(tpm_hsap), Stages_hsap = factor(str_remove_all(colnames(tpm_hsap), "_[0-9]+")))
coldata_mmul <- data.frame(row.names=colnames(tpm_mmul), Stages_mmul = factor(str_remove_all(colnames(tpm_mmul), "_[0-9]+")))
coldata_mmus <- data.frame(row.names=colnames(tpm_mmus), Stages_mmus = factor(str_remove_all(colnames(tpm_mmus), "_[0-9]+")))
coldata_nvec <- data.frame(row.names=colnames(tpm_nvec), Stages_nvec = factor(str_remove_all(colnames(tpm_nvec), "_[0-9]+")))
coldata_pcau <- data.frame(row.names=colnames(tpm_pcau), Stages_pcau = factor(str_remove_all(colnames(tpm_pcau), "_[0-9]+")))
coldata_scar <- data.frame(row.names=colnames(tpm_scar), Stages_scar = factor(str_remove_all(colnames(tpm_scar), "_[0-9]+")))
coldata_sfel <- data.frame(row.names=colnames(tpm_sfel), Stages_sfel = factor(str_remove_all(colnames(tpm_sfel), "_[0-9]+")))
coldata_tcas <- data.frame(row.names=colnames(tpm_tcas), Stages_tcas = factor(str_remove_all(colnames(tpm_tcas), "_[0-9]+")))
coldata_xtro <- data.frame(row.names=colnames(tpm_xtro), Stages_xtro = factor(str_remove_all(colnames(tpm_xtro), "_[0-9]+")))
coldata_hduj <- data.frame(row.names=colnames(tpm_hduj), Stages_hduj = factor(str_remove_all(colnames(tpm_hduj), "_[0-9]+")))
coldata_cint <- data.frame(row.names=colnames(tpm_cint), Stages_cint = factor(str_remove_all(colnames(tpm_cint), "_[0-9]+")))
coldata_hery <- data.frame(row.names=colnames(tpm_hery), Stages_hery = factor(str_remove_all(colnames(tpm_hery), "_[0-9]+")))
coldata_htub <- data.frame(row.names=colnames(tpm_htub), Stages_htub = factor(str_remove_all(colnames(tpm_htub), "_[0-9]+")))
coldata_lvar <- data.frame(row.names=colnames(tpm_lvar), Stages_lvar = factor(str_remove_all(colnames(tpm_lvar), "_[0-9]+")))
coldata_mcap <- data.frame(row.names=colnames(tpm_mcap), Stages_mcap = factor(str_remove_all(colnames(tpm_mcap), "_[0-9]+")))
coldata_mfra <- data.frame(row.names=colnames(tpm_mfra), Stages_mfra = factor(str_remove_all(colnames(tpm_mfra), "_[0-9]+")))
coldata_mlei <- data.frame(row.names=colnames(tpm_mlei), Stages_mlei = factor(str_remove_all(colnames(tpm_mlei), "_[0-9]+")))
coldata_pdum <- data.frame(row.names=colnames(tpm_pdum), Stages_pdum = factor(str_remove_all(colnames(tpm_pdum), "_[0-9]+")))
coldata_pliv <- data.frame(row.names=colnames(tpm_pliv), Stages_pliv = factor(str_remove_all(colnames(tpm_pliv), "_[0-9]+")))
coldata_pmin <- data.frame(row.names=colnames(tpm_pmin), Stages_pmin = factor(str_remove_all(colnames(tpm_pmin), "_[0-9]+")))
coldata_sscr <- data.frame(row.names=colnames(tpm_sscr), Stages_sscr = factor(str_remove_all(colnames(tpm_sscr), "_[0-9]+")))
coldata_ttra <- data.frame(row.names=colnames(tpm_ttra), Stages_ttra = factor(str_remove_all(colnames(tpm_ttra), "_[0-9]+")))


##### BUILD DF WITH ORTHOGROUPS AND EXPRESSION VALUES FOR PERSISTENTLY EXPRESSED #####
variables <- ls(pattern = "tpm_....$")


for(i in 1:length(variables)){
  temp <- as.data.frame(t(get(variables[i]))) #Extracting and transforming abundance table
  temp_metadata <- get(paste("coldata_", str_remove_all(variables[i], "tpm_"), sep = "")) #Get metadata for counts
  temp$replicate <- as.factor(str_extract(rownames(temp_metadata), "[0-9]+$")) #Add replicate metadata to tbl
  temp$stage <- as.factor(temp_metadata[, 1]) #Add stage metadata to tbl
  
  temp_name <- paste("tpm_tbl_", str_remove_all(variables[i], "tpm_"), sep = "")
  
  temp <- tibble(temp) %>%
    filter(temp$stage == "stage1") %>%
    dplyr::select(-stage, -replicate)
  
  temp <- as.data.frame(t(temp[, -1])); colnames(temp) <- paste(rep( str_remove_all(variables[i], "tpm_"), dim(temp)[2]), 1:dim(temp)[2], sep = "_")
  # temp <- subset(temp, (rownames(temp) %in% maternal_tpm) & !(rownames(temp) %in% maternal)) #using maternal genes which are not in degraded category
  
  if(!(temp_name == "tpm_tbl_sfel")){
    rownames(temp) <- str_remove_all(rownames(temp), "\\.[0-9]+$") 
  }
  
  print(temp_name)
  assign(temp_name, temp)
  
}

#For C.gigas batch effect was noticed between the repeated experiments, using only NET1 experiment derived samples only
tpm_tbl_cgig <- tpm_tbl_cgig[, 3:4]

rownames(tpm_tbl_asum) <- str_remove_all(rownames(tpm_tbl_asum), "_t[0-9]+")

tpm_tbl_aste$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_aste), ogroups$Anopheles_stephensi)]
tpm_tbl_btau$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_btau), ogroups$Bos_taurus)]
tpm_tbl_bger$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_bger), ogroups$Blattella_germanica)]
tpm_tbl_chir$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_chir), ogroups$Capra_hircus)]
tpm_tbl_ggal$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_ggal), ogroups$Gallus_gallus)]
tpm_tbl_hsap$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_hsap), ogroups$Homo_sapiens)]
tpm_tbl_mmul$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_mmul), ogroups$Macaca_mulatta)]
tpm_tbl_mmus$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_mmus), ogroups$Mus_musculus)]
tpm_tbl_nvec$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_nvec), ogroups$Nematostella_vectensis)]
tpm_tbl_asum$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_asum), ogroups$Ascaris_suum)]
tpm_tbl_bmor$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_bmor), ogroups$Bombyx_mori)]
tpm_tbl_blan$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_blan), ogroups$Branchiostoma_lanceolatum)]
tpm_tbl_cana$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_cana), ogroups$Caenorhabditis_angaria)]
tpm_tbl_cele$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_cele), ogroups$Caenorhabditis_elegans)]
tpm_tbl_cgig$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_cgig), ogroups$Crassostrea_gigas)]
tpm_tbl_dana$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dana), ogroups$Drosophila_ananassae)]
tpm_tbl_dere$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dere), ogroups$Drosophila_erecta)]
tpm_tbl_dmel$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dmel), ogroups$Drosophila_melanogaster)]
tpm_tbl_dmoj$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dmoj), ogroups$Drosophila_mojavensis)]
tpm_tbl_dper$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dper), ogroups$Drosophila_persimilis)]
tpm_tbl_drer$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_drer), ogroups$Danio_rerio)]
tpm_tbl_dsim$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dsim), ogroups$Drosophila_simulans)]
tpm_tbl_dvir$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dvir), ogroups$Drosophila_virilis)]
tpm_tbl_dwil$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dwil), ogroups$Drosophila_willistoni)]
tpm_tbl_dyak$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_dyak), ogroups$Drosophila_yakuba)]
tpm_tbl_tcas$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_tcas), ogroups$Tribolium_castaneum)]
tpm_tbl_pcau$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_pcau), ogroups$Priapulus_caudatus)]
tpm_tbl_scar$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_scar), ogroups$Steinernema_carpocapsae)]
tpm_tbl_sfel$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_sfel), ogroups$Steinernema_feltiae)]
tpm_tbl_xtro$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_xtro), ogroups$Xenopus_tropicalis)]
tpm_tbl_hduj$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_hduj), ogroups$Hypsibius_dujardini)]
tpm_tbl_cint$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_cint), ogroups$Ciona_intestinalis)]
tpm_tbl_hery$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_hery), ogroups$Heliocidaris_erythrogramma)]
tpm_tbl_htub$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_htub), ogroups$Heliocidaris_tuberculata)]
tpm_tbl_lvar$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_lvar), ogroups$Lytechinus_variegatus)]
tpm_tbl_mfra$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_mfra), ogroups$Mesocentrotus_franciscanus)]
tpm_tbl_mlei$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_mlei), ogroups$Mnemiopsis_leidyi)]
tpm_tbl_mcap$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_mcap), ogroups$Montipora_capitata)]
tpm_tbl_pliv$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_pliv), ogroups$Paracentrotus_lividus)]
tpm_tbl_pmin$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_pmin), ogroups$Patiria_miniata)]
tpm_tbl_pdum$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_pdum), ogroups$Platynereis_dumerilii)]
tpm_tbl_sscr$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_sscr), ogroups$Sus_scrofa)]
tpm_tbl_ttra$ID <- ogroups$Orthogroups[match(rownames(tpm_tbl_ttra), ogroups$Terebratalia_transversa)]

#Change colnames
colnames(tpm_tbl_asum) <- c(paste(rep("Ascaris_suum", length(tpm_tbl_asum)-1), 1:(length(tpm_tbl_asum)-1), sep = "_"), "ID")
colnames(tpm_tbl_aste) <- c(paste(rep("Anopheles_stephensi", length(tpm_tbl_aste)-1), 1:(length(tpm_tbl_aste)-1), sep = "_"), "ID")
colnames(tpm_tbl_btau) <- c(paste(rep("Bos_taurus", length(tpm_tbl_btau)-1), 1:(length(tpm_tbl_btau)-1), sep = "_"), "ID")
colnames(tpm_tbl_bger) <- c(paste(rep("Blattella_germanica", length(tpm_tbl_bger)-1), 1:(length(tpm_tbl_bger)-1), sep = "_"), "ID")
colnames(tpm_tbl_bmor) <- c(paste(rep("Bombyx_mori", length(tpm_tbl_bmor)-1), 1:(length(tpm_tbl_bmor)-1), sep = "_"), "ID")
colnames(tpm_tbl_cana) <- c(paste(rep("Caenorhabditis_angaria", length(tpm_tbl_cana)-1), 1:(length(tpm_tbl_cana)-1), sep = "_"), "ID")
colnames(tpm_tbl_cele) <- c(paste(rep("Caenorhabditis_elegans", length(tpm_tbl_cele)-1), 1:(length(tpm_tbl_cele)-1), sep = "_"), "ID")
colnames(tpm_tbl_cgig) <- c(paste(rep("Crassostrea_gigas", length(tpm_tbl_cgig)-1), 1:(length(tpm_tbl_cgig)-1), sep = "_"), "ID")
colnames(tpm_tbl_chir) <- c(paste(rep("Capra_hircus", length(tpm_tbl_chir)-1), 1:(length(tpm_tbl_chir)-1), sep = "_"), "ID")
colnames(tpm_tbl_dana) <- c(paste(rep("Drosophila_ananassae", length(tpm_tbl_dana)-1), 1:(length(tpm_tbl_dana)-1), sep = "_"), "ID")
colnames(tpm_tbl_dere) <- c(paste(rep("Drosophila_erecta", length(tpm_tbl_dere)-1), 1:(length(tpm_tbl_dere)-1), sep = "_"), "ID")
colnames(tpm_tbl_dmel) <- c(paste(rep("Drosophila_melanogaster", length(tpm_tbl_dmel)-1), 1:(length(tpm_tbl_dmel)-1), sep = "_"), "ID")
colnames(tpm_tbl_dmoj) <- c(paste(rep("Drosophila_mojavensis", length(tpm_tbl_dmoj)-1), 1:(length(tpm_tbl_dmoj)-1), sep = "_"), "ID")
colnames(tpm_tbl_dper) <- c(paste(rep("Drosophila_persimilis", length(tpm_tbl_dper)-1), 1:(length(tpm_tbl_dper)-1), sep = "_"), "ID")
colnames(tpm_tbl_dsim) <- c(paste(rep("Drosophila_simulans", length(tpm_tbl_dsim)-1), 1:(length(tpm_tbl_dsim)-1), sep = "_"), "ID")
colnames(tpm_tbl_dvir) <- c(paste(rep("Drosophila_virilis", length(tpm_tbl_dvir)-1), 1:(length(tpm_tbl_dvir)-1), sep = "_"), "ID")
colnames(tpm_tbl_dyak) <- c(paste(rep("Drosophila_yakuba", length(tpm_tbl_dyak)-1), 1:(length(tpm_tbl_dyak)-1), sep = "_"), "ID")
colnames(tpm_tbl_dwil) <- c(paste(rep("Drosophila_willistoni", length(tpm_tbl_dwil)-1), 1:(length(tpm_tbl_dwil)-1), sep = "_"), "ID")
colnames(tpm_tbl_ggal) <- c(paste(rep("Gallus_gallus", length(tpm_tbl_ggal)-1), 1:(length(tpm_tbl_ggal)-1), sep = "_"), "ID")
colnames(tpm_tbl_hsap) <- c(paste(rep("Homo_sapiens", length(tpm_tbl_hsap)-1), 1:(length(tpm_tbl_hsap)-1), sep = "_"), "ID")
colnames(tpm_tbl_mmul) <- c(paste(rep("Macaca_mulatta", length(tpm_tbl_mmul)-1), 1:(length(tpm_tbl_mmul)-1), sep = "_"), "ID")
colnames(tpm_tbl_mmus) <- c(paste(rep("Mus_musculus", length(tpm_tbl_mmus)-1), 1:(length(tpm_tbl_mmus)-1), sep = "_"), "ID")
colnames(tpm_tbl_nvec) <- c(paste(rep("Nematostella_vectensis", length(tpm_tbl_nvec)-1), 1:(length(tpm_tbl_nvec)-1), sep = "_"), "ID")
colnames(tpm_tbl_sfel) <- c(paste(rep("Steinernema_feltiae", length(tpm_tbl_sfel)-1), 1:(length(tpm_tbl_sfel)-1), sep = "_"), "ID")
colnames(tpm_tbl_scar) <- c(paste(rep("Steinernema_carpocapsae", length(tpm_tbl_scar)-1), 1:(length(tpm_tbl_scar)-1), sep = "_"), "ID")
colnames(tpm_tbl_pcau) <- c(paste(rep("Priapulus_caudatus", length(tpm_tbl_pcau)-1), 1:(length(tpm_tbl_pcau)-1), sep = "_"), "ID")
colnames(tpm_tbl_tcas) <- c(paste(rep("Tribolium_castaneum", length(tpm_tbl_tcas)-1), 1:(length(tpm_tbl_tcas)-1), sep = "_"), "ID")
colnames(tpm_tbl_xtro) <- c(paste(rep("Xenopus_tropicalis", length(tpm_tbl_xtro)-1), 1:(length(tpm_tbl_xtro)-1), sep = "_"), "ID")
colnames(tpm_tbl_blan) <- c(paste(rep("Branchiostoma_lanceolatum", length(tpm_tbl_blan)-1), 1:(length(tpm_tbl_blan)-1), sep = "_"), "ID")
colnames(tpm_tbl_drer) <- c(paste(rep("Danio_rerio", length(tpm_tbl_drer)-1), 1:(length(tpm_tbl_drer)-1), sep = "_"), "ID")
colnames(tpm_tbl_hduj) <- c(paste(rep("Hypsibius_dujardini", length(tpm_tbl_hduj)-1), 1:(length(tpm_tbl_hduj)-1), sep = "_"), "ID")
colnames(tpm_tbl_cint) <- c(paste(rep("Ciona_intestinalis", length(tpm_tbl_cint)-1), 1:(length(tpm_tbl_cint)-1), sep = "_"), "ID")
colnames(tpm_tbl_hery) <- c(paste(rep("Heliocidaris_erythrogramma", length(tpm_tbl_hery)-1), 1:(length(tpm_tbl_hery)-1), sep = "_"), "ID")
colnames(tpm_tbl_htub) <- c(paste(rep("Heliocidaris_tuberculata", length(tpm_tbl_htub)-1), 1:(length(tpm_tbl_htub)-1), sep = "_"), "ID")
colnames(tpm_tbl_lvar) <- c(paste(rep("Lytechinus_variegatus", length(tpm_tbl_lvar)-1), 1:(length(tpm_tbl_lvar)-1), sep = "_"), "ID")
colnames(tpm_tbl_mfra) <- c(paste(rep("Mesocentrotus_franciscanus", length(tpm_tbl_mfra)-1), 1:(length(tpm_tbl_mfra)-1), sep = "_"), "ID")
colnames(tpm_tbl_mlei) <- c(paste(rep("Mnemiopsis_leidyi", length(tpm_tbl_mlei)-1), 1:(length(tpm_tbl_mlei)-1), sep = "_"), "ID")
colnames(tpm_tbl_mcap) <- c(paste(rep("Montipora_capitata", length(tpm_tbl_mcap)-1), 1:(length(tpm_tbl_mcap)-1), sep = "_"), "ID")
colnames(tpm_tbl_pliv) <- c(paste(rep("Paracentrotus_lividus", length(tpm_tbl_pliv)-1), 1:(length(tpm_tbl_pliv)-1), sep = "_"), "ID")
colnames(tpm_tbl_pmin) <- c(paste(rep("Patiria_miniata", length(tpm_tbl_pmin)-1), 1:(length(tpm_tbl_pmin)-1), sep = "_"), "ID")
colnames(tpm_tbl_pdum) <- c(paste(rep("Platynereis_dumerilii", length(tpm_tbl_pdum)-1), 1:(length(tpm_tbl_pdum)-1), sep = "_"), "ID")
colnames(tpm_tbl_sscr) <- c(paste(rep("Sus_scrofa", length(tpm_tbl_sscr)-1), 1:(length(tpm_tbl_sscr)-1), sep = "_"), "ID")
colnames(tpm_tbl_ttra) <- c(paste(rep("Terebratalia_transversa", length(tpm_tbl_ttra)-1), 1:(length(tpm_tbl_ttra)-1), sep = "_"), "ID")

#Extracting unassigned genes
unassigned_lst <- list(asum = tpm_tbl_asum[is.na(tpm_tbl_asum$ID),],
                       aste = tpm_tbl_aste[is.na(tpm_tbl_aste$ID),],
                       bmor = tpm_tbl_bmor[is.na(tpm_tbl_bmor$ID),],
                       chir = tpm_tbl_chir[is.na(tpm_tbl_chir$ID),],
                       ggal = tpm_tbl_ggal[is.na(tpm_tbl_ggal$ID),],
                       hsap = tpm_tbl_hsap[is.na(tpm_tbl_hsap$ID),],
                       mmul = tpm_tbl_mmul[is.na(tpm_tbl_mmul$ID),],
                       mmus = tpm_tbl_mmus[is.na(tpm_tbl_mmus$ID),],
                       nvec = tpm_tbl_nvec[is.na(tpm_tbl_nvec$ID),],
                       blan = tpm_tbl_blan[is.na(tpm_tbl_blan$ID),],
                       btau = tpm_tbl_btau[is.na(tpm_tbl_btau$ID),],
                       cana = tpm_tbl_cana[is.na(tpm_tbl_cana$ID),],
                       bger = tpm_tbl_bger[is.na(tpm_tbl_bger$ID),],
                       cele = tpm_tbl_cele[is.na(tpm_tbl_cele$ID),],
                       cgig = tpm_tbl_cgig[is.na(tpm_tbl_cgig$ID),],
                       dana = tpm_tbl_dana[is.na(tpm_tbl_dana$ID),],
                       dere = tpm_tbl_dere[is.na(tpm_tbl_dere$ID),],
                       dmel = tpm_tbl_dmel[is.na(tpm_tbl_dmel$ID),],
                       dmoj = tpm_tbl_dmoj[is.na(tpm_tbl_dmoj$ID),],
                       dper = tpm_tbl_dper[is.na(tpm_tbl_dper$ID),],
                       dsim = tpm_tbl_dsim[is.na(tpm_tbl_dsim$ID),],
                       drer = tpm_tbl_drer[is.na(tpm_tbl_drer$ID),],
                       dvir = tpm_tbl_dvir[is.na(tpm_tbl_dvir$ID),],
                       dwil = tpm_tbl_dwil[is.na(tpm_tbl_dwil$ID),],
                       dyak = tpm_tbl_dyak[is.na(tpm_tbl_dyak$ID),],
                       scar = tpm_tbl_scar[is.na(tpm_tbl_scar$ID),],
                       sfel = tpm_tbl_sfel[is.na(tpm_tbl_sfel$ID),],
                       tcas = tpm_tbl_tcas[is.na(tpm_tbl_tcas$ID),],
                       pcau = tpm_tbl_pcau[is.na(tpm_tbl_pcau$ID),],
                       xtro= tpm_tbl_xtro[is.na(tpm_tbl_xtro$ID),],
                       hduj= tpm_tbl_hduj[is.na(tpm_tbl_hduj$ID),])


X <- vector(mode="list", length=length(variables))

for(i in 1:length(variables)){
  X[[i]] <- get(variables[[i]])$ID
  names(X)[i] <- str_remove_all(variables[[i]], "tpm_tbl_")
}

upset_obj <- fromList(X)

write.table(upset_obj, "~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/OG_presence.tbl",
            quote = F, col.names = T, row.names = F, sep = "\t")

#Discarding unassigned OGs and normalizing within species
variables <- ls(pattern = "tpm_tbl_....$")

for(i in 1:length(variables)){
  temp <- get(variables[i])
  libsize <- colSums(temp[, 1:2])[1]
  temp_name <- paste(variables[i], "clean", sep = "_")
  
  print(c(dim(temp[!is.na(temp$ID),])[1] / dim(temp)[1], i))
  temp <- temp[!is.na(temp$ID),]
  temp[, 1:(dim(temp)[2]-1)] <- applyNormFactorsMat(temp[, 1:dim(temp)[2]-1], 
                                                    calcNormFactors.default(temp[, 1:dim(temp)[2]-1],
                                                                            method = "TMM", 
                                                                            lib.size = rep(libsize, dim(temp)[2]-1))) 
  assign(temp_name, temp)
}



#Collapsing all tables into one data frame 
#If multiple paralogs in one orthogroup I take the median of all genes in the respective group
#Using median because within orthogroups there might be outliers (f.ex. a small proportion of paralogs might be not expressed at all, while most are)
OG_df <- full_join(tpm_tbl_asum_clean, tpm_tbl_aste_clean, by="ID") %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Anopheles_stephensi_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_blan_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Branchiostoma_lanceolatum_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_btau_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Bos_taurus_4, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_bmor_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Bombyx_mori_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_bger_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Blattella_germanica_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_cana_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Caenorhabditis_angaria_4, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_cele_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Caenorhabditis_elegans_4, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_cgig_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Crassostrea_gigas_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_chir_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Capra_hircus_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dana_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_ananassae_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dere_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_erecta_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dmel_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_melanogaster_5, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dmoj_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_mojavensis_6, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dper_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_persimilis_6, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dsim_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_simulans_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dvir_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_virilis_7, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dwil_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_willistoni_6, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_dyak_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Drosophila_yakuba_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_drer_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Danio_rerio_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_hsap_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Homo_sapiens_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_mmul_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Macaca_mulatta_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_mmus_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Mus_musculus_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_nvec_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Nematostella_vectensis_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_ggal_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Gallus_gallus_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_tcas_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Tribolium_castaneum_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_scar_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Steinernema_carpocapsae_4, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_sfel_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Steinernema_feltiae_4, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_pcau_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Priapulus_caudatus_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_xtro_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Xenopus_tropicalis_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_hduj_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Hypsibius_dujardini_4, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_cint_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Ciona_intestinalis_2, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_hery_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Heliocidaris_erythrogramma_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_htub_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Heliocidaris_tuberculata_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_lvar_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Lytechinus_variegatus_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_mfra_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Mesocentrotus_franciscanus_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_mlei_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Mnemiopsis_leidyi_4, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_mcap_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Montipora_capitata_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_pliv_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Paracentrotus_lividus_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_pmin_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Patiria_miniata_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_pdum_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Platynereis_dumerilii_7, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_sscr_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Sus_scrofa_3, median, na.rm = T))
OG_df <- full_join(OG_df, tpm_tbl_ttra_clean, by="ID")  %>% tibble() %>% group_by(ID) %>% dplyr::summarize(across(Ascaris_suum_1:Terebratalia_transversa_2, median, na.rm = T))


ids <- OG_df[,1]
OG_df <- OG_df[!(is.na(ids$ID)), ]
ids <- ids[!is.na(ids)]
OG_df <- as.data.frame(OG_df[, -1])
rownames(OG_df) <- ids

intDF <- collapseReplicate(OG_df)
discrDF <- apply(intDF, 2, function(x){
  x <- as.matrix(do.call(cbind, tibble(dplyr::case_when((x < 2) ~ "1", 
                                                        (x >= 2) & (x < 10) ~ "2", 
                                                        (x > 10) & (x < 1000) ~ "3", 
                                                        (x >= 1000) ~ "4"))))
})
rownames(discrDF) <- rownames(intDF)

binaryDF <- intDF
binaryDF[!(binaryDF == "NaN")] <- 2
binaryDF[binaryDF == "NaN"] <- 1


#Remove orthogroups where less than 10 species have expression values and log transform
exprsDF <- collapseReplicate(OG_df)
OG_df_se <- SEcollapseReplicate(OG_df)
keep <- rownames(exprsDF[apply(exprsDF, 1, function(x){sum(!is.na(x))}) >= 10, ])
OG_df <- log(OG_df + 0.01)
discrDF <- discrDF[keep, ]

#Batch correction
mod <- read_tsv('~/batch_stage1.txt', col_names = F)
mod$X1 <- str_remove_all(mod$X1, "_[0-9]+")

mod1 <- model.matrix(~X1, data = mod)
bat.1 <- limma::removeBatchEffect(x = OG_df, batch = mod$X2, design = mod1)
bat.1 <- limma::removeBatchEffect(x = bat.1, batch = mod$X3, design = mod1)
bat.1 <- limma::removeBatchEffect(x = bat.1, batch = mod$X4, design = mod1)
bat.1 <- limma::removeBatchEffect(x = bat.1, batch = mod$X5, design = mod1)


write.table(bat.1, "~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/batch_OG_DF_norm-final.tsv", 
            quote=FALSE, sep='\t',  col.names = NA)
write.table(discrDF, "~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/discr_OG_DF_norm-final.tsv", 
            quote=FALSE, sep='\t',  col.names = NA)
write.table(binaryDF, "~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/binary_OG_DF_norm-final.tsv", 
            quote=FALSE, sep='\t',  col.names = NA)
write.table(OG_df_se, "~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/OG_DF_se_norm-final.tsv", 
            quote=FALSE, sep='\t',  col.names = NA)
saveRDS(unassigned_lst, "~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/unassigned_lst.RDS")

##### CATEGORIZING #######

rm(list = ls(pattern = "bat.1"))
rm(list = ls(pattern = "*DF*"))
rm(list = ls(pattern = "*OG_df*"))
rm(list = ls(pattern = "tpm_tbl_...._clean$"))
rm(list = ls(pattern = "unassigned*"))
rm(list = ls(pattern = "coldata*"))

cat_df <- data.frame(ID = as.character(),
                     Category = as.character(),
                     species = as.character())
#Looping through species specific tables to add categories
for(i in 1:length(variables)){
  temp <- get(variables[i])
  temp_name <- paste(variables[i], "cat", sep = "_")
  print(paste0("Working on: ", temp_name))
  
  temp$gID <- rownames(temp)
  temp_sp <- str_remove_all(colnames(temp)[1], "_[0-9]+$")
  temp <- tibble(temp) %>% 
    dplyr::select(gID, ID) %>% 
    dplyr::group_by(gID) %>% 
    mutate(category = dplyr::case_when((gID %in% maternal) ~ "D", 
                                     (gID %in% maternal_tpm) ~ "M", 
                                       !(gID %in% maternal) & !(gID %in% maternal_tpm) ~ "NA")) %>%
    ungroup() %>%
    dplyr::select(-gID)
  colnames(temp)[colnames(temp) == "category"] <- temp_sp
  
  #saveRDS(temp, paste0("~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/Intermediate_files/", paste0(temp_sp, ".RDS")))
  #temp <- readRDS(paste0("~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/Intermediate_files/", paste0(temp_sp, ".RDS")))
  temp$species <- temp_sp
  colnames(temp)[2] <- "Category"
  temp[temp[, 2] == "NA", 2] <- NA
  temp <- temp[!(is.na(temp$ID)), ]
  cat_df <- rbind(cat_df, temp)
  
  print(paste0("Done: ", temp_name))
}

cat_df <- unique(cat_df)


write.table(cat_df, "~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/OG_categories.tsv")

save.image("~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/OG_tbl.RData")
#load("~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/OG_tbl.RData")
