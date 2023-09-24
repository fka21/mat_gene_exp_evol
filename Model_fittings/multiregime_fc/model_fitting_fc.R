#Script for estimating evolutionary models in fold change data
#As input it requires the dated species tree, the gene expression data, and a table containing the IDs of maternally expressed genes
#It will output the models themselves and a summary table of the model fitting in ./output/ directory


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
setwd("/mnt/SDB/Feri/PhD/multiregime_sandbox_fc")
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
                             rep("hemotrphic viviparous", 6), 
                             "oviparity", 
                             rep("ovuliparity", 13))),
                    species_tree$tip.label)
species_tree <- make.simmap(species_tree, regimes, model = "ER")

fcMat <- t(read.table("./fc.tsv", sep = "\t", header = T, row.names = 1))
geneCat <- read.table("./OG_categories.tsv", header = T, sep = " ", row.names = 1)

#Initialize output table
res_df_fc <- data.frame(OG = character(), 
                           Model = character(),
                           aic.weight  = double(),
                           Tree_size = double(),
                           Mean_trait = double(),
                           Perm_test = double())

#Using parallel computing
registerDoParallel(5) 

#Looping through each orthogroup
for(k in 1:dim(fcMat)[2]){
  
  #Extract orthogroup with its standard errors
  temp <- fcMat[, k, drop = F]
  temp <- na.omit(temp)
  
  print(colnames(temp))
  
  #Pointless to use small trees as models will get misspecified
  #Using a tree size of 10 as cutoff
  if(length(temp[, 1]) < 15){next}
  
  #Preparing data
  species_tree_temp <- ape::keep.tip(species_tree, species_tree$tip.label[species_tree$tip.label %in% rownames(temp)])
  
  regimes_temp <- regimes[names(regimes) %in% rownames(temp)]
  regimes_temp <- regimes_temp[order(match(names(regimes_temp), rownames(temp)))]
  
  temp_df <- data.frame(rownames(temp),
                        regimes_temp,
                        temp[, 1])
  
  #Fitting different models
  fitOUM <- OUwie(species_tree_temp, temp_df, model="OUM", simmap.tree=TRUE)
  fitBM1 <- OUwie(species_tree_temp, temp_df, model="BM1", simmap.tree=TRUE)
  fitBMS <- OUwie(species_tree_temp, temp_df, model="BMS", simmap.tree=TRUE)
  fitOUMV <- OUwie(species_tree_temp, temp_df, model="OUMV", simmap.tree=TRUE)
  fitOUMA <- OUwie(species_tree_temp, temp_df, model="OUMA", simmap.tree=TRUE)
  fitOUMVA <- OUwie(species_tree_temp, temp_df, model="OUMVA", simmap.tree=TRUE)
  fitOU1 <- OUwie(species_tree_temp, temp_df, model="OU1", simmap.tree=TRUE)
  fitWhite <- fitContinuous(species_tree_temp, temp[,1], model = "white", ncores = 5)
  
  #Extracting AICc weights
  ouwie_aicc <- setNames(c(fitBM1$AICc, fitBMS$AICc, fitOUM$AICc, fitOUMV$AICc, fitOUMA$AICc, fitOUMVA$AICc, fitWhite$opt$aicc, fitOU1$AICc),
                         c("fitBM1", "fitBMS", "fitOUM", "fitOUMV", "fitOUMA", "fitOUMVA", "fitWhite", "fitOU1"))
  aic.weight <- aic.w(ouwie_aicc)

  #Permutation test  
  perm <- c()
  
  perm <- foreach (i=1:250, .combine = 'c') %dopar% {
    asd <- temp_df
    asd[, 3] <- sample(asd[, 3])
    
    if(names(which.max(aic.weight)) == "fitWhite"){
      
      asd[, 3] <- rnorm(length(asd[, 3]), mean(asd[, 3]), sd(asd[, 3]))
      return(fitContinuous(species_tree_temp, asd[,3, drop = F], model = "white", ncores = 1)$opt$aicc)
      
      saveRDS(fitWhite,
              paste0("./output/", paste0(colnames(temp), "_fc.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitBM1") {
      
      return(OUwie(species_tree_temp, asd, model="BM1", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitBM1,
              paste0("./output/", paste0(colnames(temp), "_fc.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitBMS") {
      
      return(OUwie(species_tree_temp, asd, model="BMS", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitBMS,
              paste0("./output/", paste0(colnames(temp), "_fc.RDS")))      
      
      
    } else if(names(which.max(aic.weight)) == "fitOUM") {
      
      return(OUwie(species_tree_temp, asd, model="OUM", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOUM,
              paste0("./output/", paste0(colnames(temp), "_fc.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitOUMA") {
      
      return(OUwie(species_tree_temp, asd, model="OUMA", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOUMA,
              paste0("./output/", paste0(colnames(temp), "_fc.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitOUMV") {
      
      return(OUwie(species_tree_temp, asd, model="OUMV", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOUMV,
              paste0("./output/", paste0(colnames(temp), "_fc.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitOUMVA") {
      
      return(OUwie(species_tree_temp, asd, model="OUMVA", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOUMVA,
              paste0("./output/", paste0(colnames(temp), "_fc.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitOU1") {
      
      return(OUwie(species_tree_temp, asd, model="OU1", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOU1,
              paste0("./output/", paste0(colnames(temp), "_fc.RDS")))
      
    }
  }


  #Saving the model fitting statistics to the output table  
  if(names(which.max(aic.weight)) == "fitWhite"){
    temp_aicc <- get(names(which.max(aic.weight)))$opt$aicc
  } else {
    temp_aicc <- get(names(which.max(aic.weight)))$AICc
  }  
  
  res_df_fc <- rbind(res_df_fc,
                    data.frame(OG = colnames(temp),
                               Model = names(which.max(aic.weight)),
                               aic.weight  = max(aic.weight),
                               Tree_size = dim(temp)[1],
                               Mean_trait = mean(temp[, 1]),
                               Perm_test = ks.test(temp_aicc, perm, alternative = "greater")$p.value))
  #Saving output
  write.table(res_df_fc,
              "./exprs_fc_fits.tsv",
              sep = "\t",
              quote = F,
              col.names = T,
              row.names = F)

 print(paste("Finished", k, sep = " "))
  
}



