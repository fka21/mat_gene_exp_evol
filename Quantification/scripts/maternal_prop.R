##### LIBRARIES #####
library(tidyverse)
library(RColorBrewer)
library(extrafont)
font_import()
loadfonts(device="win") 

##### FUNCTIONS #####
source("~/Documents/Gene_expr_evol/Scripts/functions.R")
cols<-brewer.pal(4,"YlGnBu")

##### READ IN DATA #####
#Iteratively read in txiimport RDS files and assign to new variables
files <- list.files("~/Documents/Gene_expr_evol/Intermediate_files", 
                    pattern = "txi_*", recursive = T, full.names = T)
files_ext <- list.files("~/Documents/Gene_expr_evol/Intermediate_files", 
                        pattern = "txi_*", recursive = T, full.names = F)

for(i in 1:length(files)){
  temp <- readRDS(files[i])
  temp_name <- str_remove_all(files_ext[i], "txi_")
  temp_name <- str_remove_all(temp_name, ".RDS")
  
  assign(paste0("tpm_maternal_", temp_name), temp$abundance)
}


#Read in ID list containing all degraded gene IDS for all speciesspecies
files <- list.files("~/Documents/Gene_expr_evol/Intermediate_files", 
                    pattern = "*_maternal_NFE-T.RDS", recursive = T, full.names = T)
files_ext <- list.files("~/Documents/Gene_expr_evol/Intermediate_files", 
                        pattern = "*_maternal_NFE-T.RDS", recursive = T, full.names = F)

for(i in 1:length(files)){
  temp <- readRDS(files[i])
  temp_name <- str_remove_all(files_ext[i], "_maternal_NFE-T.RDS")
  
  assign(paste0("dge_maternal_", temp_name), temp)
}

##### TPM DEFINITIONS #####
#This section is used to subset abundance files for maternal genes and downregulated genes respectively

#Initizalize variables
variables_tpm <- ls(pattern = "tpm_maternal_....")
variables_dge <- ls(pattern = "dge_maternal_....")
tpm_maternal_df <- data.frame(gene = character(),
                              exp = double(),
                              species = character(),
                              category = character())
dge_maternal_df <- data.frame(gene = character(),
                              exp = double(),
                              species = character(),
                              category = character())

#For loop to go over all species and collapse replicates
#Compiling all results into "maternal_df"

for(i in 1:length(variables_tpm)){
  
  print(variables_tpm[i])
  temp_name <- str_remove_all(variables_tpm[i], "tpm_maternal_")
  
  temp_df <- get(variables_tpm[i])
  
  #For these species header adjustments are required
  if(temp_name %in% c("aste", "btau", "chir", "cint","drer", "ggal", "hsap", "mmul", "mmus", "sscr", "xtro")){
    rownames(temp_df) <- str_remove_all(rownames(temp_df), "\\.[0-9]+$")
  }
  
  #Collapsing replicates into single expression values and keeping only oocyte stages
  temp <- collapseReplicate(temp_df)
  temp <- temp[, "stage1"]
  
  #Saving into appropriate result tables
  temp_df <- data.frame(gene = names(temp), exp = temp, species = str_remove_all(variables_tpm[i], "tpm_maternal_"), category = "maternal") 
  temp_df <- subset(temp_df, !(temp_df$gene %in% get(variables_dge[i])))
  tpm_maternal_df <- rbind(tpm_maternal_df, temp_df)
  
  temp_df <- data.frame(gene = names(temp), exp = temp, species = str_remove_all(variables_tpm[i], "tpm_maternal_"), category = "down_regulated")
  temp_df <- subset(temp_df, temp_df$gene %in% get(variables_dge[i]))
  dge_maternal_df <- rbind(dge_maternal_df, temp_df)

}

#Defining ranges of TPM according to https://www.ebi.ac.uk/gxa/FAQ.html (instead of TPM 1 using TPM 2 for cut-off, 
#following this article https://link.springer.com/article/10.1007/s12064-013-0178-3)

save.image("~/Documents/Gene_expr_evol/Intermediate_files/maternal_props.RData")

##### PLOT PREPARATIONS #####
library(phytools)
library(RColorBrewer)
library(evobiR)

#Read in tree and format a species name
species_tree <- read.tree("~/tree_calibration/congr_sp_tree.dated.tre")
species_tree$tip.label[species_tree$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"
species_tree$tip.label[species_tree$tip.label == "Hypsibius_dujardini"] <- "Hypsibius_exemplaris"

#Indexing table for species name abreviations and full length species names
crosref <- data.frame(code = str_remove_all(str_remove_all(variables_tpm, "tpm_maternal_"), ".RDS"),
                      names = c("Anopheles_stephensi","Ascaris_suum","Blattella_germanica","Branchiostoma_lanceolatum","Bombyx_mori","Bos_taurus",      
                                "Caenorhabditis_angaria","Caenorhabditis_elegans","Crassostrea_gigas","Capra_hircus","Ciona_intestinalis","Drosophila_ananassae",           
                                "Drosophila_erecta","Drosophila_melanogaster","Drosophila_mojavensis","Drosophila_persimilis","Danio_rerio","Drosophila_simulans","Drosophila_virilis",             
                                "Drosophila_willistoni","Drosophila_yakuba","Gallus_gallus","Hypsibius_exemplaris","Heliocidaris_erythrogramma","Homo_sapiens","Heliocidaris_tuberculata",           
                                "Lytechinus_variegatus","Montipora_capitata","Mesocentrotus_franciscanus","Mnemiopsis_leidyi","Macaca_mulatta","Mus_musculus","Nematostella_vectensis",         
                                "Priapulus_caudatus","Platynereis_dumerilii","Paracentrotus_lividus","Patiria_miniata","Steinernema_carpocapsae","Steinernema_feltiae","Sus_scrofa",                     
                                "Tribolium_castaneum","Terebratalia_transversa","Xenopus_tropicalis"))

tpm_maternal_df$species <- str_remove_all(tpm_maternal_df$species, ".RDS")
dge_maternal_df$species <- str_remove_all(dge_maternal_df$species, ".RDS")
species_tree$tip.label <- str_replace_all(species_tree$tip.label, " ", "_")

#Test if there is difference in the proportions of gene expression categories

tpm <- tpm_maternal_df %>% 
  mutate(category = dplyr::case_when((exp < 2) ~ "no_expression",
                                     (exp >= 2) & (exp < 10) ~ "low",
                                     (exp > 10) & (exp < 1000) ~ "medium",
                                     (exp >= 1000) ~ "high")) %>% 
  na.omit() %>% 
  group_by(species, category) %>% 
  dplyr:: summarise(n = n()) %>%
  ungroup() %>% 
  dplyr::select(species, category, n) %>%
  group_by(species) %>%
  mutate(total = sum(n)) %>%
  mutate(gene_category = "mat")

dge <- dge_maternal_df %>% 
  mutate(category = dplyr::case_when((exp < 2) ~ "no_expression",
                                     (exp >= 2) & (exp < 10) ~ "low",
                                     (exp > 10) & (exp < 1000) ~ "medium",
                                     (exp >= 1000) ~ "high")) %>% 
  na.omit() %>% 
  group_by(species, category) %>% 
  dplyr:: summarise(n = n()) %>%
  ungroup() %>% 
  dplyr::select(species, category, n) %>%
  group_by(species) %>%
  mutate(total = sum(n)) %>%
  mutate(gene_category = "degr")

test_df <- bind_rows(tpm, dge) %>% 
  spread(key = category, value = n, fill = 0) 

test_res <- data.frame(speces = as.character(),
                       high = as.double(),
                       medium = as.double(),
                       low = as.double(),
                       no_expression = as.double())


#Test if degraded gene categories are more abundant in highly expressed genes
for(i in 1:length(unique(test_df$species))){
  species_temp <- unique(test_df$species)[i]
  
  int_df <- test_df %>% filter(species == species_temp)
  
  test_res[i, 1] <- species_temp
  test_res[i, 2] <- prop.test(x = unlist(c(int_df[int_df$gene_category == "degr", "high"], 
                                    int_df[int_df$gene_category == "mat", "high"])),
                              n = unlist(c(int_df[int_df$gene_category == "degr", "total"],
                                    int_df[int_df$gene_category == "mat", "total"])),
                              alternative = "greater")$p.value
  test_res[i, 3] <- prop.test(x = unlist(c(int_df[int_df$gene_category == "degr", "medium"], 
                                           int_df[int_df$gene_category == "mat", "medium"])),
                              n = unlist(c(int_df[int_df$gene_category == "degr", "total"],
                                           int_df[int_df$gene_category == "mat", "total"])),
                              alternative = "greater")$p.value
  test_res[i, 4] <- prop.test(x = unlist(c(int_df[int_df$gene_category == "degr", "low"], 
                                           int_df[int_df$gene_category == "mat", "low"])),
                              n = unlist(c(int_df[int_df$gene_category == "degr", "total"],
                                           int_df[int_df$gene_category == "mat", "total"])),
                              alternative = "greater")$p.value
  test_res[i, 5] <- prop.test(x = unlist(c(int_df[int_df$gene_category == "degr", "no_expression"], 
                                           int_df[int_df$gene_category == "mat", "no_expression"])),
                              n = unlist(c(int_df[int_df$gene_category == "degr", "total"],
                                           int_df[int_df$gene_category == "mat", "total"])),
                              alternative = "greater")$p.value
  
}

#View results on prop.test
tibble(test_res[, 2:5] <= 0.05) %>% mutate(sp = test_res$speces) %>% View()

#Categorize expression values in tables and calculate proportion af given categories for maternal genes
tpm <- tpm_maternal_df %>% 
  mutate(category = dplyr::case_when((exp < 2) ~ "No expression",
                                     (exp >= 2) & (exp < 10) ~ "Low",
                                     (exp > 10) & (exp < 1000) ~ "Medium",
                                     (exp >= 1000) ~ "High")) %>% 
  na.omit() %>% 
  group_by(species, category) %>% 
  dplyr:: summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>% 
  dplyr::select(species, category, freq) %>% 
  spread(key = category, value = freq, fill = 0)

#Categorize expression values in tables and calculate proportion af given categories for degraded genes
dge <- dge_maternal_df %>% 
  mutate(category = dplyr::case_when((exp < 2) ~ "No expression",
                                     (exp >= 2) & (exp < 10) ~ "Low",
                                     (exp > 10) & (exp < 1000) ~ "Medium",
                                     (exp >= 1000) ~ "High")) %>% 
  na.omit() %>% 
  group_by(species, category) %>% 
  dplyr:: summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>% 
  dplyr::select(species, category, freq) %>% 
  spread(key = category, value = freq, fill = 0)


#Add full species name to proportion tables
dge$species <- crosref$names[match(crosref$code, dge$species)]
tpm$species <- crosref$names[match(crosref$code, tpm$species)]

#Format and order tables into matrixes required by plotting function
X <- as.matrix(tpm[,2:5])
rownames(X) <- tpm$species
X <- ReorderData(species_tree, X)
Y <- as.matrix(dge[,2:5])
rownames(Y) <- dge$species
Y <- ReorderData(species_tree, Y)

#Plotting
svg("~/Desktop/Publication_plots/maternal_proportions.svg",
    width = 8, height = 8)
par(mfrow=c(1,3), family = "Arial") ## we can also use layout
plotTree.barplot(species_tree, X, 
                 add=TRUE, 
                 pt.cex = 12, 
                 args.barplot=list(xlab="Maternal proportions",
                                   col=cols,
                                   mar=c(5.1,0,2.1,2.1)))

legend("topleft", legend = colnames(X), pch=22,
       pt.cex = 2, pt.bg=cols, cex=0.75, bty="n", horiz=T)
plotTree.barplot(species_tree, Y,args.barplot=list(xlab="Down-regulated proportions",mar=c(5.1,0,2.1,2.1),
                                                   col=cols),args.plotTree=list(plot=FALSE),add=TRUE)
legend("topleft", legend = colnames(X), pch=22,
       pt.cex = 2, pt.bg=cols, cex=0.75, bty="n", horiz=T)

dev.off()



##### SAVING DATA ######
#save.image("~/Documents/Gene_expr_evol/Intermediate_files/maternal_props.RData")
load("~/Documents/Gene_expr_evol/Intermediate_files/maternal_props.RData")
