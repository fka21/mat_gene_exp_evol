#Script used for fitting phylogenetic logistic models to fold change values
#As input it requires the fold change values in the form of a matrix, the dated species tree and
# for the output table it will also require the best fitting model from the fitting of continuous trait models
#As output it will output each fitted model and an associated regression plot in the output/ folder


##### LIBRARYIES ######
library(motmot)
library(dplyr)
library(stringr)
library(phytools)
library(phylolm)
library(future)

setwd("/mnt/data1/Feri/pglm")

##### FUNCTION ######
source("./functions.R")

##### READ DATA #####
#For reproducibility
set.seed(7)

#Phylogeny and expression data
tree <- read.tree('./congr_sp_tree.dated.tre')
tree$tip.label[tree$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"

fc <- t(read.table("./fc.tsv", sep = "\t"))

mod_fc <- read.table("./exprs_fc_fits.tsv", header = T, sep = "\t", fill = T) %>%
  as_tibble() %>%
  filter(aic.weight >= 0.5 & Perm_test <= 0.05)

fc <- fc[, colnames(fc) %in% mod_fc$OG]


res_fc <- data.frame(OG = character(),
                     Model = character(),
                     wAIC = double(),
                     LRT = double(),
                     Coeff_inter = double(),
                     Coeff_slope = double(),
                     phyCorrelation = double(),
                     Category = character(),
                     binOne = integer(),
                     binZero = integer(),
                     stringsAsFactors = FALSE)

# Define representative species for plotting
repr_full <- data.frame(vivipar = factor(c(rep(0, 19),
                                           rep(0, 4), 
                                           rep(1, 6), 
                                           0,
                                           rep(0, 13))),
                        ovulipar = factor(c(rep(0, 19),
                                            rep(1, 4), 
                                            rep(0, 6), 
                                            0, 
                                            rep(1, 13))),
                        ovipar = factor(c(rep(1, 19),
                                          rep(0, 4), 
                                          rep(0, 6), 
                                          1, 
                                          rep(0, 13))))
rownames(repr_full) <- tree$tip.label

# Set up parallel processing
plan(multicore, workers = 3)

#Iterating over reproductive modes and for each reproductive mode through each gene
for(k in 1:dim(repr_full)[2]){
  
  repr <- repr_full[, k, drop = F]
  name <- paste0("_", colnames(repr_full)[k])
  
  for(i in 1:dim(fc)[2]){
    
    temp_name <- colnames(fc[, i, drop = F])
    
    #Skipping ones without fitted model
    if(!(temp_name %in% mod_fc$OG)){next}
    
    #Preparing inputs
    sortedDat <- sortTraitData(phy = tree, y = fc[, i, drop = F], log.trait = F, pass.ultrametric = T)
    print(paste(colnames(fc[, i, drop = F]), colnames(repr_full)[k], sep = ": "))
    
    #Building data table
    df <- repr
    df$expr <- sortedDat$trait[match(rownames(df), rownames(sortedDat$trait))]
    df <- na.omit(df)
    
    #Modelling only where there is multiple reproductive modes
    if(dim(distinct(select(df, 1)))[1] < 2){next}
    
    colnames(df) <- c("repr", "expr")
    
    #Fitting null model and full model
    model.full <- phyloglm(repr ~ expr, data = df, phy = tree, 
                           method = "logistic_IG10", boot = 100)
    model.null <- phyloglm(repr ~ 1, data = df, phy = tree, 
                           method = "logistic_IG10")
    
    selection <-  MuMIn::Weights(data.frame(AICc = c(model.full$aic, model.null$aic)))
    pval <- pchisq(-2 * (model.null$logLik - model.full$logLik), df = 1, lower.tail = F)
    
    res_fc[i, 1] <- temp_name
    res_fc[i, 2] <- mod_fc$Model[match(temp_name, mod_fc$OG)]
    res_fc[i, 3] <- selection[1]
    res_fc[i, 4] <- pval
    res_fc[i, 5] <- model.full$bootmean[1]
    res_fc[i, 6] <- model.full$bootmean[2]
    res_fc[i, 7] <- model.full$bootmean[3]
    res_fc[i, 8] <- colnames(repr_full)[k]
    res_fc[i, 9] <- sum(df$repr == 1)
    res_fc[i, 10] <- sum(df$repr == 0)
    
    saveRDS(model.full, paste0("./output_fc/",
                               paste0(temp_name), paste0(name, ".RDS")))
    
    write.table(res_fc, 
                paste0("./output_fc/pglm_res_fc", paste0(name, ".tsv")),
                sep = "\t",
                col.names = T,
                row.names = F,
                quote = F)
    
    base_name <- paste0("./output_fc/",
                        temp_name)
    
    svg(paste0(base_name,  paste0(name,".svg")))
    
    par(las=1, bty="l", family = "Helvetica") ## cosmetic
    
    plot(df$expr, jitter(as.numeric(as.character(df$repr)), factor=0, amount=0.02),
         xlab = "log(fold change)", ylab = "Reproductive mode",
         cex = 1.5,     
         lwd = 2,
         cex.lab = 1.5,    # X-axis and Y-axis labels size
         cex.axis = 1.2)
    
    curve(plogis(model.full$bootmean[1] + model.full$bootmean[2]*x), 
          col="red", 
          add=TRUE, 
          lwd = 3)
    
    mtext(paste0("LRT: ", round(pval, 3)),
          adj = 1,
          side = 3,
          cex = 1.5)
    
    dev.off()
    
    
  }
  
}

save.image("./output_fc/PGLS_fc.RData")
