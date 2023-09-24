#Script for estimating evolutionary model parameters in fold change data
#As input it requires the same inputs as the model fitting script with the addition of the summary table output from that script
#As output it will output a table containing the parameters of the best fitting model for each orthogroup


##### LIBRARYIES ######
library(motmot)
library(dplyr)
library(stringr)
library(geiger)
library(phytools)
library(OUwie)
library(foreach)
library(doParallel)

##### FUNCTION ######
setwd("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/multiregime_sandbox_fc/")
source("./functions.R")

##### READ TREE #####
#For reproducibility
set.seed(7)

#Read in data
species_tree <- read.tree("./congr_sp_tree.dated.tre")
species_tree$tip.label[species_tree$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"
species_tree$tip.label[species_tree$tip.label == "Hypsibius_dujardini"] <- "Hypsibius_exemplaris"


regimes <- setNames(factor(c(rep("oviparity", 19),
                             rep("ovuliparity", 4), 
                             rep("hemotrophic_viviparous", 6), 
                             "oviparity", 
                             rep("ovuliparity", 13))),
                    species_tree$tip.label)
species_tree <- make.simmap(species_tree, regimes, model = "ER")

#Reading in data and preparing it
fcMat <- t(read.table("./fc.tsv", sep = "\t", header = T, row.names = 1))
mods <- read.table("./exprs_fc_fits.tsv", sep = "\t", header = T) %>%
  na.omit() %>%
  filter(aic.weight >= 0.5 & Perm_test <= 0.05 & !(Model %in% "fitWhite")) %>%
  mutate(Model = str_remove_all(Model, "fit"))

fcMat <- fcMat[,colnames(fcMat) %in% mods$OG]

#Defining output data frames for parameters
res <- data.frame(OG = character(),
                  Model = character(),
                  mean_expr = double(),
                  species_nr = double(),
                  sigma2 = double(),
                  alpha = double(),
                  ovuliparity.sigma2 = double(),
                  oviparity.sigma2 = double(),
                  hemotrophic_viviparous.sigma2 = double(),
                  ovuliparity.alpha = double(),
                  oviparity.alpha = double(),
                  hemotrophic_viviparous.alpha = double(),
                  ovuliparity.theta = double(),
                  oviparity.theta = double(),
                  hemotrophic_viviparous.theta = double(),
                  ovuliparity.theta.se = double(),
                  oviparity.theta.se = double(),
                  hemotrophic_viviparous.theta.se = double(),
                  theta = double(),
                  theta.anc = double(),
                  theta.se = double(),
                  theta.anc.se = double())

#Looping through fitted models
for(i in 1:length(mods$OG)){
  
  #Subsetting for fitted model
  temp_og <- mods[i, 1]
  temp_mod <- mods[i, 2]
  
  print(temp_og)
  
  temp <- fcMat[, temp_og, drop = F]
  temp <- na.omit(temp)

  #Preparing data
  species_tree_temp <- keep.tip(species_tree, species_tree$tip.label[species_tree$tip.label %in% rownames(temp)])
  regimes_temp <- regimes[names(regimes) %in% rownames(temp)]
  regimes_temp <- regimes_temp[order(match(names(regimes_temp), rownames(temp)))]
  
  temp_df <- data.frame(species = rownames(temp),
                        regime = regimes_temp,
                        expr = temp[, 1])

  #Fitting the best winning model again without the need to perform permutation test
  fit <- OUwie(species_tree_temp, temp_df, model=temp_mod, simmap.tree=TRUE)
  
  saveRDS(fit, paste0("./output/", paste0(temp_og, "_fc.RDS")))

  if(temp_mod == "BM1"){
    
    res[i, "OG"] <- temp_og
    res[i, "Model"] <- temp_mod
    res[i, "mean_expr"] <- mean(temp_df$expr)
    res[i, "species_nr"] <- length(temp_df$species)
    res[i ,"sigma2"] <- fit$solution[2, 1]
    res[i, "theta.anc"] <- fit$theta[1, 1]
    res[i, "theta"] <- fit$theta[2, 1]
    res[i, "theta.anc.se"] <- fit$theta[1, 2]
    res[i, "theta.se"] <- fit$theta[2, 2]
    
  } else if(temp_mod == "BMS"){
    
    rownames(fit$theta)[2:dim(fit$theta)[1]] <- colnames(fit$solution)
    
    if(length(unique(temp_df$regime)) == 3){

      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i, "theta.anc"] <- fit$theta[1, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "theta.anc.se"] <- fit$theta[1, 2]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "ovuliparity")) == 2){

      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i, "theta.anc"] <- fit$theta[1, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "theta.anc.se"] <- fit$theta[1, 2]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "hemotrphic viviparous")) == 2){

      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i, "theta.anc"] <- fit$theta[1, 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "theta.anc.se"] <- fit$theta[1, 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
    } else if(sum(unique(temp_df$regime) %in% c("ovuliparity", "hemotrphic viviparous")) == 2){

      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i, "theta.anc"] <- fit$theta[1, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "theta.anc.se"] <- fit$theta[1, 2]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
    }
  } else if(temp_mod == "OU1"){
    
    res[i, "OG"] <- temp_og
    res[i, "Model"] <- temp_mod
    res[i, "mean_expr"] <- mean(temp_df$expr)
    res[i, "species_nr"] <- length(temp_df$species)
    res[i ,"alpha"] <- fit$solution[1, 1]
    res[i ,"sigma2"] <- fit$solution[2, 1]
    res[i, "theta"] <- fit$theta[1, 1]
    res[i, "theta.se"] <- fit$theta[1, 2]
    
  } else if(temp_mod == "OUM"){
    
    rownames(fit$theta)[1:dim(fit$theta)[1]] <- colnames(fit$solution)
    
    if(length(unique(temp_df$regime)) == 3){

      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"alpha"] <- fit$solution[1, 1]
      res[i ,"sigma2"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]

    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "ovuliparity")) == 2){

      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"alpha"] <- fit$solution[1, 1]
      res[i ,"sigma2"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("hemotrophic_viviparous", "ovuliparity")) == 2){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"alpha"] <- fit$solution[1, 1]
      res[i ,"sigma2"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "hemotrophic_viviparous")) == 2){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"alpha"] <- fit$solution[1, 1]
      res[i ,"sigma2"] <- fit$solution[2, 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    }
    
    
  } else if(temp_mod == "OUMA"){
    
    rownames(fit$theta)[1:dim(fit$theta)[1]] <- colnames(fit$solution)
    
    if(length(unique(temp_df$regime)) == 3){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.alpha"] <- fit$solution[1, "ovuliparity"]
      res[i ,"oviparity.alpha"] <- fit$solution[1, "oviparity"]
      res[i ,"hemotrophic_viviparous.alpha"] <- fit$solution[1, "hemotrophic_viviparous"]
      res[i ,"sigma2"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "ovuliparity")) == 2){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.alpha"] <- fit$solution[1, "ovuliparity"]
      res[i ,"oviparity.alpha"] <- fit$solution[1, "oviparity"]
      res[i ,"sigma2"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]

    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "hemotrphic viviparous")) == 2){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"hemotrophic_viviparous.alpha"] <- fit$solution[1, "hemotrophic_viviparous"]
      res[i ,"oviparity.alpha"] <- fit$solution[1, "oviparity"]
      res[i ,"sigma2"] <- fit$solution[2, 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("ovuliparity", "hemotrphic viviparous")) == 2){
     
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"hemotrophic_viviparous.alpha"] <- fit$solution[1, "hemotrophic_viviparous"]
      res[i ,"ovuliparity.alpha"] <- fit$solution[1, "ovuliparity"]
      res[i ,"sigma2"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    }
    
  } else if(temp_mod == "OUMV"){
    
    rownames(fit$theta)[1:dim(fit$theta)[1]] <- colnames(fit$solution)
    
    if(length(unique(temp_df$regime)) == 3){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i ,"alpha"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "ovuliparity")) == 2){
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"alpha"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "hemotrphic viviparous")) == 2){
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"alpha"] <- fit$solution[2, 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("ovuliparity", "hemotrphic viviparous")) == 2){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"alpha"] <- fit$solution[2, 1]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    }
  } else if(temp_mod == "OUMVA"){
    
    rownames(fit$theta)[1:dim(fit$theta)[1]] <- colnames(fit$solution)
    
    if(length(unique(temp_df$regime)) == 3){
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i ,"ovuliparity.alpha"] <- fit$solution[1, "ovuliparity"]
      res[i ,"oviparity.alpha"] <- fit$solution[1, "oviparity"]
      res[i ,"hemotrophic_viviparous.alpha"] <- fit$solution[1, "hemotrophic_viviparous"]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "ovuliparity")) == 2){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"ovuliparity.alpha"] <- fit$solution[1, "ovuliparity"]
      res[i ,"oviparity.alpha"] <- fit$solution[1, "oviparity"]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("oviparity", "hemotrphic viviparous")) == 2){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i ,"oviparity.sigma2"] <- fit$solution[2, "oviparity"]
      res[i ,"oviparity.alpha"] <- fit$solution[1, "oviparity"]
      res[i ,"hemotrophic_viviparous.alpha"] <- fit$solution[1, "hemotrophic_viviparous"]
      res[i, "oviparity.theta"] <- fit$theta["oviparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "oviparity.theta.se"] <- fit$theta["oviparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    } else if(sum(unique(temp_df$regime) %in% c("ovuliparity", "hemotrphic viviparous")) == 2){
      
      res[i, "OG"] <- temp_og
      res[i, "Model"] <- temp_mod
      res[i, "mean_expr"] <- mean(temp_df$expr)
      res[i, "species_nr"] <- length(temp_df$species)
      res[i ,"hemotrophic_viviparous.sigma2"] <- fit$solution[2, "hemotrophic_viviparous"]
      res[i ,"ovuliparity.sigma2"] <- fit$solution[2, "ovuliparity"]
      res[i ,"ovuliparity.alpha"] <- fit$solution[1, "ovuliparity"]
      res[i ,"hemotrophic_viviparous.alpha"] <- fit$solution[1, "hemotrophic_viviparous"]
      res[i, "ovuliparity.theta"] <- fit$theta["ovuliparity", 1]
      res[i, "hemotrophic_viviparous.theta"] <- fit$theta["hemotrophic_viviparous", 1]
      res[i, "ovuliparity.theta.se"] <- fit$theta["ovuliparity", 2]
      res[i, "hemotrophic_viviparous.theta.se"] <- fit$theta["hemotrophic_viviparous", 2]
      
    }
    
  }
  
}

write.table(res, "./output/fc_param_output.tsv",
            sep = "\t", col.names = T, row.names = F,
            quote = F)
