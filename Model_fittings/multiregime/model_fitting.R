#Script for estimating evolutionary models in gene expression data
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
setwd("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/Model_fittings/multiregime_sandbox/")
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
                             rep("hemotrophic viviparous", 6), 
                             "oviparity", 
                             rep("ovuliparity", 13))),
                    species_tree$tip.label)
species_tree <- make.simmap(species_tree, regimes, model = "ER")

#Export image of painted regimes
svg("~/Desktop/Publication_plots/species_tree_regimes.svg")
 
plotSimmap(species_tree, 
     color = setNames(brewer.pal(3, "Set2"), c("hemotrophic viviparous", "oviparity", "ovuliparity")))

dev.off()

#Read in data
exprsMat <- read.table("./batch_OG_DF_norm-final.tsv", sep = "\t", header = T, row.names = 1)
geneCat <- read.table("./OG_categories.tsv", header = T, sep = " ", row.names = 1)

#Collapse replicates before fittings
exprsMatSE <- t(SEcollapseReplicate(exprsMat))
exprsMat <- t(collapseReplicate(exprsMat))
exprsMat[is.nan(exprsMat)] <- NA
exprsMatSE[is.nan(exprsMatSE)] <- NA
rownames(exprsMat)[rownames(exprsMat) == "Hypsibius_dujardini"] <- "Hypsibius_exemplaris"
rownames(exprsMatSE)[rownames(exprsMatSE) == "Hypsibius_dujardini"] <- "Hypsibius_exemplaris"

#Initialize output table
res_df_exprs <- data.frame(OG = character(), 
                           Model = character(),
                           aic.weight  = double(),
                           Tree_size = double(),
                           Mean_trait = double(),
                            Perm_test = double())


#Performing parallel computing
registerDoParallel(5) 

#Loop through each orthogroup
for(k in 1:dim(exprsMat)[2]){
  
  #Extract orthogroup with its standard errors
    temp <- exprsMat[, k, drop = F]
    temp_se <- exprsMatSE[, k, drop = F]
  
  #Looping through maternal gene categories to prune trees and fit models
    index <- tibble(geneCat) %>%
      filter(ID == colnames(temp) & Category == "M") %>% 
      unique()
  
  #Subsetting for OGs found in maternal gene category
    temp <- as.matrix(subset(temp, rownames(temp) %in% index$species))
    temp_se_int <- as.vector(subset(temp_se, rownames(temp_se) %in% index$species))
    names(temp_se_int) <- rownames(subset(temp_se, rownames(temp_se) %in% index$species))
  
  #Some species have lost expression values, keeping these will artificially shift model fittings towards OU models
  #No standard errors are attributable to lost expression so using this as a subsetting criteria
    temp <- subset(temp, !(rownames(temp) %in% names(temp_se_int[temp_se_int == 0])))
    temp_se_int <- temp_se_int[temp_se_int != 0]
  
  #Pointless to use small trees as models will get misspecified
  #Using a tree size of 15 as cutoff
    if(length(temp[, 1]) < 15){next}
  
  #Preparing data
    species_tree_temp <- keep.tip(species_tree, species_tree$tip.label[species_tree$tip.label %in% rownames(temp)])
    regimes_temp <- regimes[names(regimes) %in% rownames(temp)]
    regimes_temp <- regimes_temp[order(match(names(regimes_temp), rownames(temp)))]
  
    temp_df <- data.frame(rownames(temp),
                          regimes_temp,
                          temp[, 1],
                          temp_se_int)
    
    #Fitting evolutionary models
    fitOUM <- OUwie(species_tree_temp, temp_df, model="OUM", simmap.tree=TRUE, mserr = "known")
    fitBM1 <- OUwie(species_tree_temp, temp_df, model="BM1", simmap.tree=TRUE, mserr = "known")
    fitBMS <- OUwie(species_tree_temp, temp_df, model="BMS", simmap.tree=TRUE, mserr = "known")
    fitOUMV <- OUwie(species_tree_temp, temp_df, model="OUMV", simmap.tree=TRUE, mserr = "known")
    fitOUMA <- OUwie(species_tree_temp, temp_df, model="OUMA", simmap.tree=TRUE, mserr = "known")
    fitOUMVA <- OUwie(species_tree_temp, temp_df, model="OUMVA", simmap.tree=TRUE, mserr = "known")
    fitOU1 <- OUwie(species_tree_temp, temp_df, model="OU1", simmap.tree=TRUE, mserr = "known")
    fitWhite <- fitContinuous(species_tree_temp, temp[,1], model = "white", SE = temp_se_int)
  
  #Extracting AICc of each model fit
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
              paste0("./output/", paste0(colnames(temp), "_exprs.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitBM1") {
      
      return(OUwie(species_tree_temp, asd, model="BM1", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitBM1,
              paste0("./output/", paste0(colnames(temp), "_exprs.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitBMS") {
      
      return(OUwie(species_tree_temp, asd, model="BMS", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitBMS,
              paste0("./output/", paste0(colnames(temp), "_exprs.RDS")))      
      
      
    } else if(names(which.max(aic.weight)) == "fitOUM") {
      
      return(OUwie(species_tree_temp, asd, model="OUM", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOUM,
              paste0("./output/", paste0(colnames(temp), "_exprs.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitOUMA") {
      
      return(OUwie(species_tree_temp, asd, model="OUMA", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOUMA,
              paste0("./output/", paste0(colnames(temp), "_exprs.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitOUMV") {
      
      return(OUwie(species_tree_temp, asd, model="OUMV", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOUMV,
              paste0("./output/", paste0(colnames(temp), "_exprs.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitOUMVA") {
      
      return(OUwie(species_tree_temp, asd, model="OUMVA", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOUMVA,
              paste0("./output/", paste0(colnames(temp), "_exprs.RDS")))
      
    } else if(names(which.max(aic.weight)) == "fitOU1") {
      
      return(OUwie(species_tree_temp, asd, model="OU1", simmap.tree=TRUE)$AICc)
      
      saveRDS(fitOU1,
              paste0("./output/", paste0(colnames(temp), "_exprs.RDS")))
      
    }
  }
   

  #Saving best winning model to output table
  if(names(which.max(aic.weight)) == "fitWhite"){
    temp_aicc <- get(names(which.max(aic.weight)))$opt$aicc
  } else {
    temp_aicc <- get(names(which.max(aic.weight)))$AICc
  }  
  
  res_df_exprs <- rbind(res_df_exprs,
                    data.frame(OG = colnames(temp),
                               Model = names(which.max(aic.weight)),
                               aic.weight  = max(aic.weight),
                               Tree_size = dim(temp)[1],
                               Mean_trait = mean(temp[, 1]),
                               Perm_test = ks.test(temp_aicc, perm, alternative = "greater")$p.value))
   
  #Exporting output table
  write.table(res_df_exprs,
              "./exprs_cont_fits.tsv",
              sep = "\t",
              quote = F,
              col.names = T,
              row.names = F)

 print(paste("Finished", k, sep = " "))
  
}



