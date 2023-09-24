

##### LOAD PACKAGES #####
Libraries <- c("DESeq2", "pheatmap", "vsn", "hexbin", "cowplot", "apeglm", "tximport", 
               "viridis", "PoiClaClu", "genefilter", "vidger", "reshape2", "ggplot2", "stringr", 
               "gplots", "sva", "tidyverse", "RColorBrewer")
lapply(Libraries, library, character.only = TRUE)

##### FUNCTION FOR MZT #######
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
source("~/Documents/Gene_expr_evol/Scripts/functions.R")


##### LOAD QUANTIFICATION ######
#Ordered species codes for naming
species <- c("aste", "asum", "bger", "blan", "bmor", "btau", "cana", "cele", "cgig", "chir", "cint", "dana", "dere", "dmel",
             "dmoj", "dper", "drer", "dsim", "dvir", "dwil", "dyak", "ggal", "hduj", "hery", "hsap", "htub", "lvar", "mcap", 
             "mfra", "mlei", "mmul", "mmus", "nvec", "pcau", "pdum", "pliv", "pmin", "scar", "sscr", "sfel", "tcas", "xtro")

##### IMPORT SALMON COUNTS #####
files <- list.files("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files", 
                    pattern = "txi_*", recursive = T, full.names = T)
files_ext <- list.files("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files", 
                    pattern = "txi_*", recursive = T, full.names = F)

for(i in 1:length(files)){
  temp <- readRDS(files[i])
  temp_name <- str_remove_all(files_ext[i], "txi_")
  temp_name <- str_remove_all(temp_name, ".RDS$")
  
  assign(paste0("txi_", temp_name), temp)
}

##### FORMAT COUNT TABLES #####

variables <- ls(pattern = "txi_*")

#Iteratively creating metadata for all species
for(i in 1:length(variables)){
  temp <- get(variables[i])
  temp_name <- str_remove_all(variables[i], "txi_")
  
  temp_df <- data.frame(row.names = colnames(temp$abundance), V1 = factor(str_remove_all(colnames(temp$abundance), "_[0-9]+")))
  colnames(temp_df) <- paste0("Stages_", temp_name)
  
  assign(paste0("coldata_", temp_name), temp_df)
}

#Adding further metadata where known
coldata_aste$var_aste <- c(rep("single", 4), rep("paired", 5))
coldata_bmor$var_bmor <- c(rep(c("female", "male"), 4), rep("unknown", 3))
coldata_cgig$var_cgig <- c(rep("NET3_exp", 4), rep("NET1_exp", 4))
coldata_dyak$var_dyak <- c(rep("male", 3), rep("female", 3), rep("unknown", 2))
coldata_mmul$var_mmul <- c("multiple", "multiple", "female1", "female1", "female2", "female1", "female1", "female3", "female1", "female1", "female1", "female2",
                           "female2", "female2", "female2", "female4")
coldata_xtro$var_xtro <- c(rep("polyA", 6), rep("ribozero", 5), rep("polyA", 6))
coldata_mlei$var_mlei <- c("0001", "0002", "0003", "0004", "0042", "0043", "0044")

##### CREATE DDS OBJECtS ####

dds_asum <- DESeqDataSetFromTximport(txi_asum, colData = coldata_asum, design = ~Stages_asum)
dds_aste <- DESeqDataSetFromTximport(txi_aste, colData = coldata_aste, design = ~Stages_aste + var_aste)
dds_bger <- DESeqDataSetFromTximport(txi_bger, colData = coldata_bger, design = ~Stages_bger)
dds_btau <- DESeqDataSetFromTximport(txi_btau, colData = coldata_btau, design = ~Stages_btau)
dds_chir <- DESeqDataSetFromTximport(txi_chir, colData = coldata_chir, design = ~Stages_chir)
dds_ggal <- DESeqDataSetFromTximport(txi_ggal, colData = coldata_ggal, design = ~Stages_ggal)
dds_hsap <- DESeqDataSetFromTximport(txi_hsap, colData = coldata_hsap, design = ~Stages_hsap)
dds_mmul <- DESeqDataSetFromTximport(txi_mmul, colData = coldata_mmul, design = ~Stages_mmul + var_mmul)
dds_mmus <- DESeqDataSetFromTximport(txi_mmus, colData = coldata_mmus, design = ~Stages_mmus)
dds_nvec <- DESeqDataSetFromTximport(txi_nvec, colData = coldata_nvec, design = ~Stages_nvec)
dds_tcas <- DESeqDataSetFromTximport(txi_tcas, colData = coldata_tcas, design = ~Stages_tcas)
dds_blan <- DESeqDataSetFromTximport(txi_blan, colData = coldata_blan, design = ~Stages_blan)
dds_cgig <- DESeqDataSetFromTximport(txi_cgig, colData = coldata_cgig, design = ~Stages_cgig + var_cgig)
dds_drer <- DESeqDataSetFromTximport(txi_drer, colData = coldata_drer, design = ~Stages_drer)
dds_xtro <- DESeqDataSetFromTximport(txi_xtro, colData = coldata_xtro, design = ~Stages_xtro + var_xtro)
dds_bmor <- DESeqDataSetFromTximport(txi_bmor, colData = coldata_bmor, design = ~Stages_bmor)
dds_cana <- DESeqDataSetFromTximport(txi_cana, colData = coldata_cana, design = ~Stages_cana)
dds_cele <- DESeqDataSetFromTximport(txi_cele, colData = coldata_cele, design = ~Stages_cele)
dds_dana <- DESeqDataSetFromTximport(txi_dana, colData = coldata_dana, design = ~Stages_dana)
dds_dere <- DESeqDataSetFromTximport(txi_dere, colData = coldata_dere, design = ~Stages_dere)
dds_dmel <- DESeqDataSetFromTximport(txi_dmel, colData = coldata_dmel, design = ~Stages_dmel)
dds_dmoj <- DESeqDataSetFromTximport(txi_dmoj, colData = coldata_dmoj, design = ~Stages_dmoj)
dds_dper <- DESeqDataSetFromTximport(txi_dper, colData = coldata_dper, design = ~Stages_dper)
dds_dsim <- DESeqDataSetFromTximport(txi_dsim, colData = coldata_dsim, design = ~Stages_dsim)
dds_dvir <- DESeqDataSetFromTximport(txi_dvir, colData = coldata_dvir, design = ~Stages_dvir)
dds_dyak <- DESeqDataSetFromTximport(txi_dyak, colData = coldata_dyak, design = ~Stages_dyak)
dds_dwil <- DESeqDataSetFromTximport(txi_dwil, colData = coldata_dwil, design = ~Stages_dwil)
dds_scar <- DESeqDataSetFromTximport(txi_scar, colData = coldata_scar, design = ~Stages_scar)
dds_sfel <- DESeqDataSetFromTximport(txi_sfel, colData = coldata_sfel, design = ~Stages_sfel)
dds_pcau <- DESeqDataSetFromTximport(txi_pcau, colData = coldata_pcau, design = ~Stages_pcau)
dds_hduj <- DESeqDataSetFromTximport(txi_hduj, colData = coldata_hduj, design = ~Stages_hduj)
dds_cint <- DESeqDataSetFromTximport(txi_cint, colData = coldata_cint, design = ~Stages_cint)
dds_pdum <- DESeqDataSetFromTximport(txi_pdum, colData = coldata_pdum, design = ~Stages_pdum)
dds_pliv <- DESeqDataSetFromTximport(txi_pliv, colData = coldata_pliv, design = ~Stages_pliv)
dds_mfra <- DESeqDataSetFromTximport(txi_mfra, colData = coldata_mfra, design = ~Stages_mfra)
dds_mcap <- DESeqDataSetFromTximport(txi_mcap, colData = coldata_mcap, design = ~Stages_mcap)
dds_hery <- DESeqDataSetFromTximport(txi_hery, colData = coldata_hery, design = ~Stages_hery)
dds_htub <- DESeqDataSetFromTximport(txi_htub, colData = coldata_htub, design = ~Stages_htub)
dds_lvar <- DESeqDataSetFromTximport(txi_lvar, colData = coldata_lvar, design = ~Stages_lvar)
dds_sscr <- DESeqDataSetFromTximport(txi_sscr, colData = coldata_sscr, design = ~Stages_sscr)
dds_mlei <- DESeqDataSetFromTximport(txi_mlei, colData = coldata_mlei, design = ~Stages_mlei)
dds_ttra <- DESeqDataSetFromTximport(txi_ttra, colData = coldata_ttra, design = ~Stages_ttra)
dds_pmin <- DESeqDataSetFromTximport(txi_pmin, colData = coldata_pmin, design = ~Stages_pmin)

#Freeing up space
rm(list = ls(pattern = "txi_....$"))

##### FILTERING LOW COUNTS OUT #####

variables <- ls(pattern = "dds_")

for(i in 1:length(variables)){
  temp <- get(variables[i])
  temp_name <- variables[i]
  keep <- rowSums(counts(temp)) >= 10
  temp <- temp[keep,]
  assign(temp_name, temp)
}

##### TRANSFORMATIONS #####

variables <- ls(pattern = "dds_....$")
variables_rld <- str_replace_all(variables, "dds", "rld")
variables_vst <- str_replace_all(variables, "dds", "vst")

for(i in 1:length(variables)){
  temp <- get(variables[i])
  temp_name <- variables_rld[i]
  temp <- rlogTransformation(temp)
  assign(temp_name, temp)
}

for(i in 1:length(variables)){
  temp <- get(variables[i])
  temp_name <- variables_vst[i]
  temp <- vst(temp)
  assign(temp_name, temp)
}

#Exporting plots of PCA and meanSD
for(i in 1:length(variables_rld)){
  
  p <- plotPCA(get(variables_rld[i]), intgroup=paste0("Stages_", str_remove_all(variables_rld[i], "rld_")))
  ggsave(plot = p, filename = paste0("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Results/Plot/PCA/", paste0(str_remove_all(variables_rld[i], "rld_"), "_rld_pca.png")))
  p <- plotPCA(get(variables_vst[i]), intgroup=paste0("Stages_", str_remove_all(variables_rld[i], "rld_")))
  ggsave(plot = p, filename = paste0("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Results/Plot/PCA/", paste0(str_remove_all(variables_rld[i], "rld_"), "_vst_pca.png")))
  p <- meanSdPlot(assay(get(variables_rld[i])))
  ggsave(filename = paste0("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Results/Plot/MeanSD/", paste0(str_remove_all(variables_rld[i], "rld_"), "_meansd.png")))
  
}

##### MZT TIMING #####
#Estimating MZT timting (if unknown and multiple timepoint samples available) for setting up contrasts

mzt.dist(rld_mcap, 4) #latest stage is stage4 to contrast for
mzt.dist(rld_mmul, 2) #latest stage is stage4 to contrast for
mzt.dist(vst_pdum, 2) #latest stage is stage2 to contrast for, using vst as many samples available
mzt.dist(rld_pmin, 2) #latest stage is stage3 to contrast for
mzt.dist(vst_sscr, 2) #latest stage is stage3 to contrast for, using vst as many samples available
mzt.dist(vst_pliv, 3) #latest stage is stage3 to contrast for, using vst as many samples available
mzt.dist(vst_ttra, 3) #latest stage is stage3 to contrast for
mzt.dist(rld_pcau, 2) #latest stage is stage3 to contrast for
mzt.dist(rld_hduj, 2)

#Droping ambiguous samples based on sample distances
dds_btau <- dds_btau[, -(colnames(dds_btau) %in% "stage1_1")] #Clustered with stage2 instead of stage1
dds_mlei <- dds_mlei[, -c(colnames(dds_sscr) %in%  c('stage1_1', 'stage1_2'))] #Clustered way closer to stage2, instead of stage1
dds_sscr <- dds_sscr[, -(colnames(dds_sscr) %in%  c('stage2_5', 'stage2_6'))] #Removing these as they group with stage3 in distance plots


##### DGE & SVA #####

variables <- ls(pattern = "dds_....")

#Need to perform DESeq before SVA batch correction
for(i in 1:length(variables)){
  temp <- get(variables[i])
  temp_name <- variables[i]
  temp <- DESeq(temp)
  print(temp_name)
  assign(temp_name, temp)
}

#SVA batch correction
for(i in 1:length(variables)){
  #Initializing
  temp <- get(variables[i])
  temp_name <- variables[i]
  temp_name_out <- str_replace_all(temp_name, "dds_", "ddssva_")
  print(temp_name)
  
    dat <- counts(temp, normalized=TRUE)
    idx <- rowMeans(dat) > 1
    dat <- dat[idx,]
    mod <- model.matrix(as.formula(paste("~Stages_", str_remove_all(temp_name, "dds_"), sep = "")), colData(temp))
    mod0 <- model.matrix(~1, colData(temp))
    svseq <- svaseq(dat, mod, mod0)
    
    temp_list <- vector(mode = "list", dim(svseq$sv)[2])
    
    #Getting number of surrogate variables and saving them
    for(i in 1:dim(svseq$sv)[2]){
      temp_list[[i]] <- as.data.frame(svseq$sv[, i])
      colnames(temp_list[[i]]) <- paste("SV", i, sep = "")
    }
    
    ddssva <- temp
    
    #Updating colData with surrogate variables
    for(i in 1:dim(svseq$sv)[2]){
      colData(ddssva) <- cbind(colData(ddssva), temp_list[[i]])
    }
    
    #Using surrogate variables repeat DESeq
    if(dim(as.data.frame(colData(ddssva)))[2] == 2){
      design(ddssva) <- as.formula(paste(paste("~Stages_", str_remove_all(temp_name, "dds_"), sep = ""), 1:dim(svseq$sv)[2], sep = " + "))
      ddssva <- DESeq(ddssva)
      assign(temp_name_out, ddssva)} else{
        formula_temp <- paste(paste("~Stages_", str_remove_all(temp_name, "dds_"), sep = ""), paste("var_", str_remove_all(temp_name, "dds_"), sep = ""), sep = " + ")
        design(ddssva) <- as.formula(paste(paste("~Stages_", str_remove_all(temp_name, "dds_"), sep = ""), 1:dim(svseq$sv)[2], sep = " + "))
        ddssva <- DESeq(ddssva)
        assign(temp_name_out, ddssva)
      }
  }
  

##### EXTRACT RESULTS #####
contrast_aste_1 = c("Stages_aste", "stage2", "stage1")
contrast_aste_2 = c("Stages_aste", "stage3", "stage1")
contrast_aste_3 = c("Stages_aste", "stage3", "stage2")
res_aste_1 = results(ddssva_aste, contrast_aste_1, alpha = 0.05)
res_aste_2 = results(ddssva_aste, contrast_aste_2, alpha = 0.05)
res_aste_3 = results(ddssva_aste, contrast_aste_3, alpha = 0.05)
aste_maternal <- unique(c(rownames(subset(res_aste_1, res_aste_1$padj < 0.05 & res_aste_1$log2FoldChange < -2)),
                          rownames(subset(res_aste_2, res_aste_2$padj < 0.05 & res_aste_2$log2FoldChange < -2)),
                          rownames(subset(res_aste_3, res_aste_3$padj < 0.05 & res_aste_3$log2FoldChange < -2))))
aste_zygotic <- unique(c(rownames(subset(res_aste_1, res_aste_1$padj < 0.05 & res_aste_1$log2FoldChange >2)),
                         rownames(subset(res_aste_2, res_aste_2$padj < 0.05 & res_aste_2$log2FoldChange >2)),
                         rownames(subset(res_aste_3, res_aste_3$padj < 0.05 & res_aste_3$log2FoldChange >2))))
aste_maternal <- str_remove_all(aste_maternal, "\\.[0-9]+$")
aste_zygotic <- str_remove_all(aste_zygotic, "\\.[0-9]+$")

contrast_asum_1 = c("Stages_asum", "stage2", "stage1")
contrast_asum_2 = c("Stages_asum", "stage3", "stage1")
contrast_asum_3 = c("Stages_asum", "stage3", "stage2")
res_asum_1 = results(ddssva_asum, contrast_asum_1, alpha = 0.05)
res_asum_2 = results(ddssva_asum, contrast_asum_2, alpha = 0.05)
res_asum_3 = results(ddssva_asum, contrast_asum_3, alpha = 0.05)
asum_maternal <- unique(c(rownames(subset(res_asum_1, res_asum_1$padj < 0.05 & res_asum_1$log2FoldChange < -2)),
                          rownames(subset(res_asum_2, res_asum_2$padj < 0.05 & res_asum_2$log2FoldChange < -2)),
                          rownames(subset(res_asum_3, res_asum_3$padj < 0.05 & res_asum_3$log2FoldChange < -2))))
asum_zygotic <- unique(c(rownames(subset(res_asum_1, res_asum_1$padj < 0.05 & res_asum_1$log2FoldChange >2)),
                         rownames(subset(res_asum_2, res_asum_2$padj < 0.05 & res_asum_2$log2FoldChange >2)),
                         rownames(subset(res_asum_3, res_asum_3$padj < 0.05 & res_asum_3$log2FoldChange >2))))

contrast_bger_1 = c("Stages_bger", "stage3", "stage2")
contrast_bger_2 = c("Stages_bger", "stage4", "stage2")
contrast_bger_3 = c("Stages_bger", "stage4", "stage3")
res_bger_1 = results(ddssva_bger, contrast_bger_1, alpha = 0.05)
res_bger_2 = results(ddssva_bger, contrast_bger_2, alpha = 0.05)
res_bger_3 = results(ddssva_bger, contrast_bger_3, alpha = 0.05)
bger_maternal <- unique(c(rownames(subset(res_bger_1, res_bger_1$padj < 0.05 & res_bger_1$log2FoldChange < -2)),
                          rownames(subset(res_bger_2, res_bger_2$padj < 0.05 & res_bger_2$log2FoldChange < -2)),
                          rownames(subset(res_bger_3, res_bger_3$padj < 0.05 & res_bger_3$log2FoldChange < -2))))
bger_zygotic <- unique(c(rownames(subset(res_bger_1, res_bger_1$padj < 0.05 & res_bger_1$log2FoldChange >2)),
                         rownames(subset(res_bger_2, res_bger_2$padj < 0.05 & res_bger_2$log2FoldChange >2)),
                         rownames(subset(res_bger_3, res_bger_3$padj < 0.05 & res_bger_3$log2FoldChange >2))))

contrast_btau_1 = c("Stages_btau", "stage2", "stage1")
res_btau_1 = results(ddssva_btau, contrast_btau_1, alpha = 0.05)
btau_maternal <- unique(c(rownames(subset(res_btau_1, res_btau_1$padj < 0.05 & res_btau_1$log2FoldChange < -2))))
btau_zygotic <- unique(c(rownames(subset(res_btau_1, res_btau_1$padj < 0.05 & res_btau_1$log2FoldChange >2))))
btau_maternal <- str_remove_all(btau_maternal, "\\.[0-9]+$")
btau_zygotic <- str_remove_all(btau_zygotic, "\\.[0-9]+$")

contrast_chir_1 = c("Stages_chir", "stage2", "stage1")
contrast_chir_2 = c("Stages_chir", "stage3", "stage1")
contrast_chir_2 = c("Stages_chir", "stage3", "stage2")
res_chir_1 = results(ddssva_chir, contrast_chir_1, alpha = 0.05)
res_chir_2 = results(ddssva_chir, contrast_chir_2, alpha = 0.05)
res_chir_3 = results(ddssva_chir, contrast_chir_2, alpha = 0.05)
chir_maternal <- unique(c(rownames(subset(res_chir_1, res_chir_1$padj < 0.05 & res_chir_1$log2FoldChange < -2)),
                          rownames(subset(res_chir_2, res_chir_2$padj < 0.05 & res_chir_2$log2FoldChange < -2)),
                          rownames(subset(res_chir_3, res_chir_3$padj < 0.05 & res_chir_3$log2FoldChange < -2))))
chir_zygotic <- unique(c(rownames(subset(res_chir_1, res_chir_1$padj < 0.05 & res_chir_1$log2FoldChange >2)),
                         rownames(subset(res_chir_2, res_chir_2$padj < 0.05 & res_chir_2$log2FoldChange >2)),
                         rownames(subset(res_chir_3, res_chir_3$padj < 0.05 & res_chir_3$log2FoldChange < -2))))
chir_maternal <- str_remove_all(chir_maternal, "\\.[0-9]+$")
chir_zygotic <- str_remove_all(chir_zygotic, "\\.[0-9]+$")

contrast_hsap_1 = c("Stages_hsap", "stage2", "stage1")
contrast_hsap_2 = c("Stages_hsap", "stage3", "stage1")
contrast_hsap_3 = c("Stages_hsap", "stage4", "stage1")
contrast_hsap_4 = c("Stages_hsap", "stage3", "stage2")
contrast_hsap_5 = c("Stages_hsap", "stage4", "stage3")
res_hsap_1 = results(ddssva_hsap, contrast_hsap_1, alpha = 0.05)
res_hsap_2 = results(ddssva_hsap, contrast_hsap_2, alpha = 0.05)
res_hsap_3 = results(ddssva_hsap, contrast_hsap_3, alpha = 0.05)
res_hsap_4 = results(ddssva_hsap, contrast_hsap_4, alpha = 0.05)
res_hsap_5 = results(ddssva_hsap, contrast_hsap_5, alpha = 0.05)
hsap_maternal <- unique(c(rownames(subset(res_hsap_1, res_hsap_1$padj < 0.05 & res_hsap_1$log2FoldChange < -2)),
                          rownames(subset(res_hsap_2, res_hsap_2$padj < 0.05 & res_hsap_2$log2FoldChange < -2)),
                          rownames(subset(res_hsap_3, res_hsap_3$padj < 0.05 & res_hsap_3$log2FoldChange < -2)),
                          rownames(subset(res_hsap_4, res_hsap_4$padj < 0.05 & res_hsap_3$log2FoldChange < -2)),
                          rownames(subset(res_hsap_5, res_hsap_5$padj < 0.05 & res_hsap_4$log2FoldChange < -2))))
hsap_zygotic <- unique(c(rownames(subset(res_hsap_1, res_hsap_1$padj < 0.05 & res_hsap_1$log2FoldChange >2)),
                         rownames(subset(res_hsap_2, res_hsap_2$padj < 0.05 & res_hsap_2$log2FoldChange >2)),
                         rownames(subset(res_hsap_3, res_hsap_3$padj < 0.05 & res_hsap_3$log2FoldChange >2)),
                         rownames(subset(res_hsap_4, res_hsap_3$padj < 0.05 & res_hsap_4$log2FoldChange >2)),
                         rownames(subset(res_hsap_5, res_hsap_3$padj < 0.05 & res_hsap_5$log2FoldChange >2))))
hsap_maternal <- str_remove_all(hsap_maternal, "\\.[0-9]+$")
hsap_zygotic <- str_remove_all(hsap_zygotic, "\\.[0-9]+$")

contrast_mmul_1 = c("Stages_mmul", "stage2", "stage1")
contrast_mmul_2 = c("Stages_mmul", "stage3", "stage1")
contrast_mmul_3 = c("Stages_mmul", "stage4", "stage1")
contrast_mmul_4 = c("Stages_mmul", "stage3", "stage2")
contrast_mmul_5 = c("Stages_mmul", "stage4", "stage3")
res_mmul_1 = results(ddssva_mmul, contrast_mmul_1, alpha = 0.05)
res_mmul_2 = results(ddssva_mmul, contrast_mmul_2, alpha = 0.05)
res_mmul_3 = results(ddssva_mmul, contrast_mmul_3, alpha = 0.05)
res_mmul_4 = results(ddssva_mmul, contrast_mmul_4, alpha = 0.05)
res_mmul_5 = results(ddssva_mmul, contrast_mmul_5, alpha = 0.05)
mmul_maternal <- unique(c(rownames(subset(res_mmul_1, res_mmul_1$padj < 0.05 & res_mmul_1$log2FoldChange < -2)),
                          rownames(subset(res_mmul_2, res_mmul_2$padj < 0.05 & res_mmul_2$log2FoldChange < -2)),
                          rownames(subset(res_mmul_3, res_mmul_3$padj < 0.05 & res_mmul_3$log2FoldChange < -2)),
                          rownames(subset(res_mmul_4, res_mmul_4$padj < 0.05 & res_mmul_4$log2FoldChange < -2)),
                          rownames(subset(res_mmul_5, res_mmul_5$padj < 0.05 & res_mmul_5$log2FoldChange < -2))))
mmul_zygotic <- unique(c(rownames(subset(res_mmul_1, res_mmul_1$padj < 0.05 & res_mmul_1$log2FoldChange >2)),
                        rownames(subset(res_mmul_2, res_mmul_2$padj < 0.05 & res_mmul_2$log2FoldChange >2)),
                         rownames(subset(res_mmul_3, res_mmul_3$padj < 0.05 & res_mmul_3$log2FoldChange >2)),
                         rownames(subset(res_mmul_4, res_mmul_4$padj < 0.05 & res_mmul_4$log2FoldChange >2)),
                         rownames(subset(res_mmul_5, res_mmul_5$padj < 0.05 & res_mmul_5$log2FoldChange >2))))
mmul_maternal <- str_remove_all(mmul_maternal, "\\.[0-9]+$")
mmul_zygotic <- str_remove_all(mmul_zygotic, "\\.[0-9]+$")

contrast_ggal_1 = c("Stages_ggal", "stage2", "stage1")
contrast_ggal_2 = c("Stages_ggal", "stage3", "stage1")
contrast_ggal_3 = c("Stages_ggal", "stage4", "stage1")
contrast_ggal_4 = c("Stages_ggal", "stage3", "stage2")
contrast_ggal_5 = c("Stages_ggal", "stage4", "stage3")
res_ggal_1 = results(ddssva_ggal, contrast_ggal_1, alpha = 0.05)
res_ggal_2 = results(ddssva_ggal, contrast_ggal_2, alpha = 0.05)
res_ggal_3 = results(ddssva_ggal, contrast_ggal_3, alpha = 0.05)
res_ggal_4 = results(ddssva_ggal, contrast_ggal_4, alpha = 0.05)
res_ggal_5 = results(ddssva_ggal, contrast_ggal_5, alpha = 0.05)
ggal_maternal <- unique(c(rownames(subset(res_ggal_1, res_ggal_1$padj < 0.05 & res_ggal_1$log2FoldChange < -2)),
                          rownames(subset(res_ggal_2, res_ggal_2$padj < 0.05 & res_ggal_2$log2FoldChange < -2)),
                          rownames(subset(res_ggal_3, res_ggal_3$padj < 0.05 & res_ggal_3$log2FoldChange < -2)),
                          rownames(subset(res_ggal_4, res_ggal_4$padj < 0.05 & res_ggal_4$log2FoldChange < -2)),
                          rownames(subset(res_ggal_5, res_ggal_5$padj < 0.05 & res_ggal_5$log2FoldChange < -2))))
ggal_zygotic <- unique(c(rownames(subset(res_ggal_1, res_ggal_1$padj < 0.05 & res_ggal_1$log2FoldChange >2)),
                         rownames(subset(res_ggal_2, res_ggal_2$padj < 0.05 & res_ggal_2$log2FoldChange >2)),
                         rownames(subset(res_ggal_3, res_ggal_3$padj < 0.05 & res_ggal_3$log2FoldChange >2)),
                         rownames(subset(res_ggal_4, res_ggal_4$padj < 0.05 & res_ggal_4$log2FoldChange >2)),
                         rownames(subset(res_ggal_5, res_ggal_5$padj < 0.05 & res_ggal_5$log2FoldChange >2))))
ggal_maternal <- str_remove_all(ggal_maternal, "\\.[0-9]+$")
ggal_zygotic <- str_remove_all(ggal_zygotic, "\\.[0-9]+$")

contrast_mmus_1 = c("Stages_mmus", "stage2", "stage1")
contrast_mmus_2 = c("Stages_mmus", "stage3", "stage1")
contrast_mmus_3 = c("Stages_mmus", "stage3", "stage2")
res_mmus_1 = results(ddssva_mmus, contrast_mmus_1, alpha = 0.05)
res_mmus_2 = results(ddssva_mmus, contrast_mmus_2, alpha = 0.05)
res_mmus_3 = results(ddssva_mmus, contrast_mmus_3, alpha = 0.05)
mmus_maternal <- unique(c(rownames(subset(res_mmus_1, res_mmus_1$padj < 0.05 & res_mmus_1$log2FoldChange < -2)),
                          rownames(subset(res_mmus_2, res_mmus_2$padj < 0.05 & res_mmus_2$log2FoldChange < -2)),
                          rownames(subset(res_mmus_3, res_mmus_3$padj < 0.05 & res_mmus_3$log2FoldChange < -2))))
mmus_zygotic <- unique(c(rownames(subset(res_mmus_1, res_mmus_1$padj < 0.05 & res_mmus_1$log2FoldChange >2)),
                         rownames(subset(res_mmus_2, res_mmus_2$padj < 0.05 & res_mmus_2$log2FoldChange >2)),
                         rownames(subset(res_mmus_3, res_mmus_3$padj < 0.05 & res_mmus_3$log2FoldChange >2))))
mmus_maternal <- str_remove_all(mmus_maternal, "\\.[0-9]+$")
mmus_zygotic <- str_remove_all(mmus_zygotic, "\\.[0-9]+$")

contrast_nvec_1 = c("Stages_nvec", "stage2", "stage1")
contrast_nvec_2 = c("Stages_nvec", "stage3", "stage1")
contrast_nvec_3 = c("Stages_nvec", "stage3", "stage2")
res_nvec_1 = results(dds_nvec, contrast_nvec_1, alpha = 0.05)
res_nvec_2 = results(dds_nvec, contrast_nvec_2, alpha = 0.05)
res_nvec_3 = results(dds_nvec, contrast_nvec_3, alpha = 0.05)
nvec_maternal <- unique(c(rownames(subset(res_nvec_1, res_nvec_1$padj < 0.05 & res_nvec_1$log2FoldChange < -2)),
                          rownames(subset(res_nvec_2, res_nvec_2$padj < 0.05 & res_nvec_2$log2FoldChange < -2)),
                          rownames(subset(res_nvec_3, res_nvec_3$padj < 0.05 & res_nvec_3$log2FoldChange < -2))))
nvec_zygotic <- unique(c(rownames(subset(res_nvec_1, res_nvec_1$padj < 0.05 & res_nvec_1$log2FoldChange >2)),
                         rownames(subset(res_nvec_2, res_nvec_2$padj < 0.05 & res_nvec_2$log2FoldChange >2)),
                         rownames(subset(res_nvec_3, res_nvec_3$padj < 0.05 & res_nvec_3$log2FoldChange >2))))

contrast_xtro_1 = c("Stages_xtro", "stage2", "stage1")
contrast_xtro_2 = c("Stages_xtro", "stage3", "stage1")
contrast_xtro_3 = c("Stages_xtro", "stage4", "stage1")
contrast_xtro_4 = c("Stages_xtro", "stage5", "stage1")
contrast_xtro_5 = c("Stages_xtro", "stage3", "stage2")
contrast_xtro_6 = c("Stages_xtro", "stage4", "stage3")
contrast_xtro_7 = c("Stages_xtro", "stage5", "stage4")
res_xtro_1 = results(ddssva_xtro, contrast_xtro_1, alpha = 0.05)
res_xtro_2 = results(ddssva_xtro, contrast_xtro_2, alpha = 0.05)
res_xtro_3 = results(ddssva_xtro, contrast_xtro_3, alpha = 0.05)
res_xtro_4 = results(ddssva_xtro, contrast_xtro_4, alpha = 0.05)
res_xtro_5 = results(ddssva_xtro, contrast_xtro_5, alpha = 0.05)
res_xtro_6 = results(ddssva_xtro, contrast_xtro_6, alpha = 0.05)
res_xtro_7 = results(ddssva_xtro, contrast_xtro_7, alpha = 0.05)
xtro_maternal <- unique(c(rownames(subset(res_xtro_1, res_xtro_1$padj < 0.05 & res_xtro_1$log2FoldChange < -2)),
                          rownames(subset(res_xtro_2, res_xtro_2$padj < 0.05 & res_xtro_2$log2FoldChange < -2)),
                          rownames(subset(res_xtro_3, res_xtro_3$padj < 0.05 & res_xtro_3$log2FoldChange < -2)),
                          rownames(subset(res_xtro_4, res_xtro_4$padj < 0.05 & res_xtro_4$log2FoldChange < -2)),
                          rownames(subset(res_xtro_5, res_xtro_5$padj < 0.05 & res_xtro_5$log2FoldChange < -2)),
                          rownames(subset(res_xtro_6, res_xtro_6$padj < 0.05 & res_xtro_6$log2FoldChange < -2)),
                          rownames(subset(res_xtro_7, res_xtro_7$padj < 0.05 & res_xtro_7$log2FoldChange < -2))))
xtro_zygotic <- unique(c(rownames(subset(res_xtro_1, res_xtro_1$padj < 0.05 & res_xtro_1$log2FoldChange >2)),
                         rownames(subset(res_xtro_2, res_xtro_2$padj < 0.05 & res_xtro_2$log2FoldChange >2)),
                         rownames(subset(res_xtro_3, res_xtro_3$padj < 0.05 & res_xtro_3$log2FoldChange >2)),
                         rownames(subset(res_xtro_4, res_xtro_4$padj < 0.05 & res_xtro_4$log2FoldChange >2)),
                         rownames(subset(res_xtro_5, res_xtro_5$padj < 0.05 & res_xtro_5$log2FoldChange >2)),
                         rownames(subset(res_xtro_6, res_xtro_6$padj < 0.05 & res_xtro_6$log2FoldChange >2)),
                         rownames(subset(res_xtro_7, res_xtro_7$padj < 0.05 & res_xtro_7$log2FoldChange >2))))
xtro_maternal <- str_remove_all(xtro_maternal, "\\.[0-9]+$")
xtro_zygotic <- str_remove_all(xtro_zygotic, "\\.[0-9]+$")

contrast_drer_1 = c("Stages_drer", "stage2", "stage1")
contrast_drer_2 = c("Stages_drer", "stage3", "stage1")
contrast_drer_3 = c("Stages_drer", "stage4", "stage1")
contrast_drer_4 = c("Stages_drer", "stage5", "stage1")
contrast_drer_5 = c("Stages_drer", "stage3", "stage2")
contrast_drer_6 = c("Stages_drer", "stage4", "stage3")
contrast_drer_7 = c("Stages_drer", "stage5", "stage4")
res_drer_1 = results(ddssva_drer, contrast_drer_1, alpha = 0.05)
res_drer_2 = results(ddssva_drer, contrast_drer_2, alpha = 0.05)
res_drer_3 = results(ddssva_drer, contrast_drer_3, alpha = 0.05)
res_drer_4 = results(ddssva_drer, contrast_drer_4, alpha = 0.05)
res_drer_5 = results(ddssva_drer, contrast_drer_5, alpha = 0.05)
res_drer_6 = results(ddssva_drer, contrast_drer_6, alpha = 0.05)
res_drer_7 = results(ddssva_drer, contrast_drer_7, alpha = 0.05)
drer_maternal <- unique(c(rownames(subset(res_drer_1, res_drer_1$padj < 0.05 & res_drer_1$log2FoldChange < -2)),
                          rownames(subset(res_drer_2, res_drer_2$padj < 0.05 & res_drer_2$log2FoldChange < -2)),
                          rownames(subset(res_drer_3, res_drer_3$padj < 0.05 & res_drer_3$log2FoldChange < -2)),
                          rownames(subset(res_drer_4, res_drer_4$padj < 0.05 & res_drer_4$log2FoldChange < -2)),
rownames(subset(res_drer_5, res_drer_5$padj < 0.05 & res_drer_5$log2FoldChange < -2)),
rownames(subset(res_drer_6, res_drer_6$padj < 0.05 & res_drer_6$log2FoldChange < -2)),
rownames(subset(res_drer_7, res_drer_7$padj < 0.05 & res_drer_7$log2FoldChange < -2))))
drer_zygotic <- unique(c(rownames(subset(res_drer_1, res_drer_1$padj < 0.05 & res_drer_1$log2FoldChange >2)),
                         rownames(subset(res_drer_2, res_drer_2$padj < 0.05 & res_drer_2$log2FoldChange >2)),
                         rownames(subset(res_drer_3, res_drer_3$padj < 0.05 & res_drer_3$log2FoldChange >2)),
                         rownames(subset(res_drer_4, res_drer_4$padj < 0.05 & res_drer_4$log2FoldChange >2)),
rownames(subset(res_drer_5, res_drer_5$padj < 0.05 & res_drer_5$log2FoldChange >2)),
rownames(subset(res_drer_6, res_drer_6$padj < 0.05 & res_drer_6$log2FoldChange >2)),
rownames(subset(res_drer_7, res_drer_7$padj < 0.05 & res_drer_7$log2FoldChange >2))))
drer_maternal <- str_remove_all(drer_maternal, "\\.[0-9]+$")
drer_zygotic <- str_remove_all(drer_zygotic, "\\.[0-9]+$")

contrast_tcas_1 = c("Stages_tcas", "stage2", "stage1")
contrast_tcas_2 = c("Stages_tcas", "stage3", "stage1")
contrast_tcas_3 = c("Stages_tcas", "stage4", "stage1")
contrast_tcas_4 = c("Stages_tcas", "stage5", "stage1")
contrast_tcas_5 = c("Stages_tcas", "stage3", "stage2")
contrast_tcas_6 = c("Stages_tcas", "stage4", "stage3")
contrast_tcas_7 = c("Stages_tcas", "stage5", "stage4")
res_tcas_1 = results(ddssva_tcas, contrast_tcas_1, alpha = 0.05)
res_tcas_2 = results(ddssva_tcas, contrast_tcas_2, alpha = 0.05)
res_tcas_3 = results(ddssva_tcas, contrast_tcas_3, alpha = 0.05)
res_tcas_4 = results(ddssva_tcas, contrast_tcas_4, alpha = 0.05)
res_tcas_5 = results(ddssva_tcas, contrast_tcas_5, alpha = 0.05)
res_tcas_6 = results(ddssva_tcas, contrast_tcas_6, alpha = 0.05)
res_tcas_7 = results(ddssva_tcas, contrast_tcas_7, alpha = 0.05)
res_tcas_maternal_1 <- subset(res_tcas_1, padj < 0.05 & log2FoldChange < -2)
res_tcas_maternal_2 <- subset(res_tcas_2, padj < 0.05 & log2FoldChange < -2)
res_tcas_maternal_3 <- subset(res_tcas_3, padj < 0.05 & log2FoldChange < -2)
res_tcas_maternal_4 <- subset(res_tcas_4, padj < 0.05 & log2FoldChange < -2)
res_tcas_maternal_5 <- subset(res_tcas_5, padj < 0.05 & log2FoldChange < -2)
res_tcas_maternal_6 <- subset(res_tcas_6, padj < 0.05 & log2FoldChange < -2)
res_tcas_maternal_7 <- subset(res_tcas_7, padj < 0.05 & log2FoldChange < -2)
tcas_maternal <- unique(c(rownames(res_tcas_maternal_1),
                          rownames(res_tcas_maternal_2),
                          rownames(res_tcas_maternal_3),
                          rownames(res_tcas_maternal_4),
                          rownames(res_tcas_maternal_5),
                          rownames(res_tcas_maternal_6),
                          rownames(res_tcas_maternal_7)))
tcas_zygotic <- unique(c(rownames(subset(res_tcas_1, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_tcas_2, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_tcas_3, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_tcas_4, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_tcas_5, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_tcas_6, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_tcas_7, padj < 0.05 & log2FoldChange >2))))

contrast_blan_1 = c("Stages_blan", "stage2", "stage1")
contrast_blan_2 = c("Stages_blan", "stage3", "stage1")
contrast_blan_3 = c("Stages_blan", "stage3", "stage2")
res_blan_1 = results(ddssva_blan, contrast_blan_1, alpha = 0.05)
res_blan_2 = results(ddssva_blan, contrast_blan_2, alpha = 0.05)
res_blan_3 = results(ddssva_blan, contrast_blan_3, alpha = 0.05)
blan_maternal <- unique(c(rownames(subset(res_blan_1, res_blan_1$padj < 0.05 & res_blan_1$log2FoldChange < -2)),
                          rownames(subset(res_blan_2, res_blan_2$padj < 0.05 & res_blan_2$log2FoldChange < -2)),
                          rownames(subset(res_blan_3, res_blan_3$padj < 0.05 & res_blan_3$log2FoldChange < -2))))
blan_zygotic <- unique(c(rownames(subset(res_blan_1, res_blan_1$padj < 0.05 & res_blan_1$log2FoldChange >2)),
                         rownames(subset(res_blan_2, res_blan_2$padj < 0.05 & res_blan_2$log2FoldChange >2)),
                         rownames(subset(res_blan_3, res_blan_3$padj < 0.05 & res_blan_3$log2FoldChange >2))))

contrast_cgig_1 = c("Stages_cgig", "stage2", "stage1")
res_cgig_1 = results(ddssva_cgig, contrast_cgig_1, alpha = 0.05)
cgig_maternal <- unique(c(rownames(subset(res_cgig_1, res_cgig_1$padj < 0.05 & res_cgig_1$log2FoldChange < -2))))
cgig_zygotic <- unique(c(rownames(subset(res_cgig_1, res_cgig_1$padj < 0.05 & res_cgig_1$log2FoldChange >2))))

contrast_bmor_1 = c("Stages_bmor", "stage2", "stage1")
res_bmor_1 = results(dds_bmor, contrast_bmor_1, alpha = 0.05)
bmor_maternal <- unique(c(rownames(subset(res_bmor_1, padj < 0.05 & log2FoldChange < -2))))
bmor_zygotic <- unique(c(rownames(subset(res_bmor_1, padj < 0.05 & log2FoldChange >2))))


contrast_cana_1 = c("Stages_cana", "stage2", "stage1")
contrast_cana_2 = c("Stages_cana", "stage3", "stage1")
contrast_cana_3 = c("Stages_cana", "stage4", "stage1")
contrast_cana_4 = c("Stages_cana", "stage5", "stage1")
contrast_cana_5 = c("Stages_cana", "stage3", "stage2")
contrast_cana_6 = c("Stages_cana", "stage4", "stage3")
contrast_cana_7 = c("Stages_cana", "stage5", "stage4")
res_cana_1 = results(ddssva_cana, contrast_cana_1, alpha = 0.05)
res_cana_2 = results(ddssva_cana, contrast_cana_2, alpha = 0.05)
res_cana_3 = results(ddssva_cana, contrast_cana_3, alpha = 0.05)
res_cana_4 = results(ddssva_cana, contrast_cana_4, alpha = 0.05)
res_cana_5 = results(ddssva_cana, contrast_cana_5, alpha = 0.05)
res_cana_6 = results(ddssva_cana, contrast_cana_6, alpha = 0.05)
res_cana_7 = results(ddssva_cana, contrast_cana_7, alpha = 0.05)
res_cana_maternal_1 <- subset(res_cana_1, padj < 0.05 & log2FoldChange < -2)
res_cana_maternal_2 <- subset(res_cana_2, padj < 0.05 & log2FoldChange < -2)
res_cana_maternal_3 <- subset(res_cana_3, padj < 0.05 & log2FoldChange < -2)
res_cana_maternal_4 <- subset(res_cana_4, padj < 0.05 & log2FoldChange < -2)
res_cana_maternal_5 <- subset(res_cana_5, padj < 0.05 & log2FoldChange < -2)
res_cana_maternal_6 <- subset(res_cana_6, padj < 0.05 & log2FoldChange < -2)
res_cana_maternal_7 <- subset(res_cana_7, padj < 0.05 & log2FoldChange < -2)
cana_maternal <- unique(c(rownames(res_cana_maternal_1),
                          rownames(res_cana_maternal_2),
                          rownames(res_cana_maternal_3),
                          rownames(res_cana_maternal_4),
                          rownames(res_cana_maternal_5),
                          rownames(res_cana_maternal_6),
                          rownames(res_cana_maternal_7)))
cana_zygotic <- unique(c(rownames(subset(res_cana_1, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cana_2, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cana_3, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cana_4, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cana_5, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cana_6, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cana_7, padj < 0.05 & log2FoldChange >2))))
contrast_cele_1 = c("Stages_cele", "stage2", "stage1")
contrast_cele_2 = c("Stages_cele", "stage3", "stage1")
contrast_cele_3 = c("Stages_cele", "stage4", "stage1")
contrast_cele_4 = c("Stages_cele", "stage5", "stage1")
contrast_cele_5 = c("Stages_cele", "stage3", "stage2")
contrast_cele_6 = c("Stages_cele", "stage4", "stage3")
contrast_cele_7 = c("Stages_cele", "stage5", "stage4")
res_cele_1 = results(ddssva_cele, contrast_cele_1, alpha = 0.05)
res_cele_2 = results(ddssva_cele, contrast_cele_2, alpha = 0.05)
res_cele_3 = results(ddssva_cele, contrast_cele_3, alpha = 0.05)
res_cele_4 = results(ddssva_cele, contrast_cele_4, alpha = 0.05)
res_cele_5 = results(ddssva_cele, contrast_cele_5, alpha = 0.05)
res_cele_6 = results(ddssva_cele, contrast_cele_6, alpha = 0.05)
res_cele_7 = results(ddssva_cele, contrast_cele_7, alpha = 0.05)
res_cele_maternal_1 <- subset(res_cele_1, padj < 0.05 & log2FoldChange < -2)
res_cele_maternal_2 <- subset(res_cele_2, padj < 0.05 & log2FoldChange < -2)
res_cele_maternal_3 <- subset(res_cele_3, padj < 0.05 & log2FoldChange < -2)
res_cele_maternal_4 <- subset(res_cele_4, padj < 0.05 & log2FoldChange < -2)
res_cele_maternal_5 <- subset(res_cele_5, padj < 0.05 & log2FoldChange < -2)
res_cele_maternal_6 <- subset(res_cele_6, padj < 0.05 & log2FoldChange < -2)
res_cele_maternal_7 <- subset(res_cele_7, padj < 0.05 & log2FoldChange < -2)
cele_maternal <- unique(c(rownames(res_cele_maternal_1),
                          rownames(res_cele_maternal_2),
                          rownames(res_cele_maternal_3),
                          rownames(res_cele_maternal_4),
                          rownames(res_cele_maternal_5),
                          rownames(res_cele_maternal_6),
                          rownames(res_cele_maternal_7)))
cele_zygotic <- unique(c(rownames(subset(res_cele_1, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cele_2, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cele_3, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cele_4, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cele_5, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cele_6, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_cele_7, padj < 0.05 & log2FoldChange >2))))


contrast_dana = c("Stages_dana", "stage2", "stage1")
res_dana = results(ddssva_dana, contrast_dana, alpha = 0.05)
dana_maternal <- unique(c(as.character(rownames(subset(res_dana, padj < 0.05 & log2FoldChange < -2)))))
dana_zygotic <- unique(c(as.character(rownames(subset(res_dana, padj < 0.05 & log2FoldChange >2)))))

contrast_dere = c("Stages_dere", "stage2", "stage1")
res_dere = results(ddssva_dere, contrast_dere, alpha = 0.05)
dere_zygotic <- unique(c(as.character(rownames(subset(res_dere, padj < 0.05 & log2FoldChange >2)))))
dere_maternal <- unique(c(as.character(rownames(subset(res_dere, padj < 0.05 & log2FoldChange < -2)))))

contrast_dmoj = c("Stages_dmoj", "stage2", "stage1")
res_dmoj = results(ddssva_dmoj, contrast_dmoj, alpha = 0.05)
dmoj_zygotic <- unique(c(as.character(rownames(subset(res_dmoj, padj < 0.05 & log2FoldChange >2)))))
dmoj_maternal <- unique(c(as.character(rownames(subset(res_dmoj, padj < 0.05 & log2FoldChange < -2)))))

contrast_dper = c("Stages_dper", "stage2", "stage1")
res_dper = results(ddssva_dper, contrast_dper, alpha = 0.05)
dper_maternal <- unique(c(as.character(rownames(subset(res_dper, padj < 0.05 & log2FoldChange < -2)))))
dper_zygotic <- unique(c(as.character(rownames(subset(res_dper, padj < 0.05 & log2FoldChange >2)))))

contrast_dyak = c("Stages_dyak", "stage2", "stage1")
res_dyak = results(ddssva_dyak, contrast_dyak, alpha = 0.05)
dyak_maternal <- unique(c(as.character(rownames(subset(res_dyak, padj < 0.05 & log2FoldChange < -2)))))
dyak_zygotic <- unique(c(as.character(rownames(subset(res_dyak, padj < 0.05 & log2FoldChange >2)))))

contrast_dsim = c("Stages_dsim", "stage2", "stage1")
res_dsim = results(ddssva_dsim, contrast_dsim, alpha = 0.05)
dsim_zygotic <- unique(c(as.character(rownames(subset(res_dsim, padj < 0.05 & log2FoldChange >2)))))
dsim_maternal <- unique(c(as.character(rownames(subset(res_dsim, padj < 0.05 & log2FoldChange < -2)))))

contrast_dvir = c("Stages_dvir", "stage2", "stage1")
res_dvir = results(ddssva_dvir, contrast_dvir, alpha = 0.05)
dvir_maternal <- unique(c(as.character(rownames(subset(res_dvir, padj < 0.05 & log2FoldChange < -2)))))
dvir_zygotic <- unique(c(as.character(rownames(subset(res_dvir, padj < 0.05 & log2FoldChange >2)))))

contrast_dwil = c("Stages_dwil", "stage2", "stage1")
res_dwil = results(ddssva_dwil, contrast_dwil, alpha = 0.05)
dwil_maternal <- unique(c(as.character(rownames(subset(res_dwil, padj < 0.05 & log2FoldChange < -2)))))
dwil_zygotic <- unique(c(as.character(rownames(subset(res_dwil, padj < 0.05 & log2FoldChange >2)))))

contrast_dmel_1 = c("Stages_dmel", "stage2", "stage1")
res_dmel_1 = results(ddssva_dmel, contrast_dmel_1, alpha = 0.05)
dmel_maternal <- unique(c(rownames(subset(res_dmel_1, padj < 0.05 & log2FoldChange < -2))))
dmel_zygotic <- unique(c(rownames(subset(res_dmel_1, padj < 0.05 & log2FoldChange >2))))

contrast_hduj_1 = c("Stages_hduj", "stage2", "stage1")
res_hduj_1 = results(ddssva_hduj, contrast_hduj_1, alpha = 0.05)
hduj_maternal <- unique(c(rownames(subset(res_hduj_1, padj < 0.05 & log2FoldChange < -2))))
hduj_zygotic <- unique(c(rownames(subset(res_hduj_1, padj < 0.05 & log2FoldChange >2))))

contrast_pcau_1 = c("Stages_pcau", "stage2", "stage1")
contrast_pcau_2 = c("Stages_pcau", "stage3", "stage1")
contrast_pcau_3 = c("Stages_pcau", "stage3", "stage2")
res_pcau_1 = results(ddssva_pcau, contrast_pcau_1, alpha = 0.05)
res_pcau_2 = results(ddssva_pcau, contrast_pcau_2, alpha = 0.05)
res_pcau_3 = results(ddssva_pcau, contrast_pcau_3, alpha = 0.05)
pcau_maternal <- unique(c(rownames(subset(res_pcau_1, res_pcau_1$padj < 0.05 & res_pcau_1$log2FoldChange < -2)),
                          rownames(subset(res_pcau_2, res_pcau_2$padj < 0.05 & res_pcau_2$log2FoldChange < -2)),
                          rownames(subset(res_pcau_3, res_pcau_3$padj < 0.05 & res_pcau_3$log2FoldChange < -2))))
pcau_zygotic <- unique(c(rownames(subset(res_pcau_1, res_pcau_1$padj < 0.05 & res_pcau_1$log2FoldChange >2)),
                         rownames(subset(res_pcau_2, res_pcau_2$padj < 0.05 & res_pcau_2$log2FoldChange >2)),
                         rownames(subset(res_pcau_3, res_pcau_3$padj < 0.05 & res_pcau_3$log2FoldChange >2))))

contrast_scar_1 = c("Stages_scar", "stage2", "stage1")
contrast_scar_2 = c("Stages_scar", "stage3", "stage1")
contrast_scar_3 = c("Stages_scar", "stage4", "stage1")
contrast_scar_4 = c("Stages_scar", "stage5", "stage1")
contrast_scar_5 = c("Stages_scar", "stage3", "stage2")
contrast_scar_6 = c("Stages_scar", "stage4", "stage3")
contrast_scar_7 = c("Stages_scar", "stage5", "stage4")
res_scar_1 = results(ddssva_scar, contrast_scar_1, alpha = 0.05)
res_scar_2 = results(ddssva_scar, contrast_scar_2, alpha = 0.05)
res_scar_3 = results(ddssva_scar, contrast_scar_3, alpha = 0.05)
res_scar_4 = results(ddssva_scar, contrast_scar_4, alpha = 0.05)
res_scar_5 = results(ddssva_scar, contrast_scar_5, alpha = 0.05)
res_scar_6 = results(ddssva_scar, contrast_scar_6, alpha = 0.05)
res_scar_7 = results(ddssva_scar, contrast_scar_7, alpha = 0.05)
res_scar_maternal_1 <- subset(res_scar_1, padj < 0.05 & log2FoldChange < -2)
res_scar_maternal_2 <- subset(res_scar_2, padj < 0.05 & log2FoldChange < -2)
res_scar_maternal_3 <- subset(res_scar_3, padj < 0.05 & log2FoldChange < -2)
res_scar_maternal_4 <- subset(res_scar_4, padj < 0.05 & log2FoldChange < -2)
res_scar_maternal_5 <- subset(res_scar_5, padj < 0.05 & log2FoldChange < -2)
res_scar_maternal_6 <- subset(res_scar_6, padj < 0.05 & log2FoldChange < -2)
res_scar_maternal_7 <- subset(res_scar_7, padj < 0.05 & log2FoldChange < -2)
scar_maternal <- unique(c(rownames(res_scar_maternal_1),
                          rownames(res_scar_maternal_2),
                          rownames(res_scar_maternal_3),
                          rownames(res_scar_maternal_4),
                          rownames(res_scar_maternal_5),
                          rownames(res_scar_maternal_6),
                          rownames(res_scar_maternal_7)))
scar_zygotic <- unique(c(rownames(subset(res_scar_1, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_scar_2, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_scar_3, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_scar_4, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_scar_5, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_scar_6, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_scar_7, padj < 0.05 & log2FoldChange >2))))
contrast_sfel_1 = c("Stages_sfel", "stage2", "stage1")
contrast_sfel_2 = c("Stages_sfel", "stage3", "stage1")
contrast_sfel_3 = c("Stages_sfel", "stage4", "stage1")
contrast_sfel_4 = c("Stages_sfel", "stage5", "stage1")
contrast_sfel_5 = c("Stages_sfel", "stage3", "stage2")
contrast_sfel_6 = c("Stages_sfel", "stage4", "stage3")
contrast_sfel_7 = c("Stages_sfel", "stage5", "stage4")
res_sfel_1 = results(ddssva_sfel, contrast_sfel_1, alpha = 0.05)
res_sfel_2 = results(ddssva_sfel, contrast_sfel_2, alpha = 0.05)
res_sfel_3 = results(ddssva_sfel, contrast_sfel_3, alpha = 0.05)
res_sfel_4 = results(ddssva_sfel, contrast_sfel_4, alpha = 0.05)
res_sfel_5 = results(ddssva_sfel, contrast_sfel_5, alpha = 0.05)
res_sfel_6 = results(ddssva_sfel, contrast_sfel_6, alpha = 0.05)
res_sfel_7 = results(ddssva_sfel, contrast_sfel_7, alpha = 0.05)
res_sfel_maternal_1 <- subset(res_sfel_1, padj < 0.05 & log2FoldChange < -2)
res_sfel_maternal_2 <- subset(res_sfel_2, padj < 0.05 & log2FoldChange < -2)
res_sfel_maternal_3 <- subset(res_sfel_3, padj < 0.05 & log2FoldChange < -2)
res_sfel_maternal_4 <- subset(res_sfel_4, padj < 0.05 & log2FoldChange < -2)
res_sfel_maternal_5 <- subset(res_sfel_5, padj < 0.05 & log2FoldChange < -2)
res_sfel_maternal_6 <- subset(res_sfel_6, padj < 0.05 & log2FoldChange < -2)
res_sfel_maternal_7 <- subset(res_sfel_7, padj < 0.05 & log2FoldChange < -2)
sfel_maternal <- unique(c(rownames(res_sfel_maternal_1),
                          rownames(res_sfel_maternal_2),
                          rownames(res_sfel_maternal_3),
                          rownames(res_sfel_maternal_4),
                          rownames(res_sfel_maternal_5),
                          rownames(res_sfel_maternal_6),
                          rownames(res_sfel_maternal_7)))
sfel_zygotic <- unique(c(rownames(subset(res_sfel_1, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_sfel_2, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_sfel_3, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_sfel_4, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_sfel_5, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_sfel_6, padj < 0.05 & log2FoldChange >2)),
                         rownames(subset(res_sfel_7, padj < 0.05 & log2FoldChange >2))))

contrast_cint = c("Stages_cint", "stage2", "stage1")
res_cint = results(ddssva_cint, contrast_cint, alpha = 0.05)
cint_maternal <- unique(c(as.character(rownames(subset(res_cint, padj < 0.05 & log2FoldChange < -2)))))
cint_zygotic <- unique(c(as.character(rownames(subset(res_cint, padj < 0.05 & log2FoldChange >2)))))
cint_maternal <- str_remove_all(cint_maternal, "\\.[0-9]+$")
cint_zygotic <- str_remove_all(cint_zygotic, "\\.[0-9]+$")


contrast_hery_1 = c("Stages_hery", "stage2", "stage1")
contrast_hery_2 = c("Stages_hery", "stage3", "stage1")
contrast_hery_3 = c("Stages_hery", "stage4", "stage1")
contrast_hery_4 = c("Stages_hery", "stage3", "stage2")
contrast_hery_5 = c("Stages_hery", "stage4", "stage3")
res_hery_1 = results(ddssva_hery, contrast_hery_1, alpha = 0.05)
res_hery_2 = results(ddssva_hery, contrast_hery_2, alpha = 0.05)
res_hery_3 = results(ddssva_hery, contrast_hery_3, alpha = 0.05)
res_hery_4 = results(ddssva_hery, contrast_hery_4, alpha = 0.05)
res_hery_5 = results(ddssva_hery, contrast_hery_5, alpha = 0.05)
hery_maternal <- unique(c(rownames(subset(res_hery_1, res_hery_1$padj < 0.05 & res_hery_1$log2FoldChange < -2)),
                          rownames(subset(res_hery_2, res_hery_2$padj < 0.05 & res_hery_2$log2FoldChange < -2)),
                          rownames(subset(res_hery_3, res_hery_3$padj < 0.05 & res_hery_3$log2FoldChange < -2)),
                          rownames(subset(res_hery_4, res_hery_4$padj < 0.05 & res_hery_4$log2FoldChange < -2)),
                          rownames(subset(res_hery_5, res_hery_5$padj < 0.05 & res_hery_5$log2FoldChange < -2))))
hery_zygotic <- unique(c(rownames(subset(res_hery_1, res_hery_1$padj < 0.05 & res_hery_1$log2FoldChange >2)),
                         rownames(subset(res_hery_2, res_hery_2$padj < 0.05 & res_hery_2$log2FoldChange >2)),
                         rownames(subset(res_hery_3, res_hery_3$padj < 0.05 & res_hery_3$log2FoldChange >2)),
                         rownames(subset(res_hery_4, res_hery_4$padj < 0.05 & res_hery_4$log2FoldChange >2)),
                         rownames(subset(res_hery_5, res_hery_5$padj < 0.05 & res_hery_5$log2FoldChange >2))))

contrast_htub_1 = c("Stages_htub", "stage2", "stage1")
contrast_htub_2 = c("Stages_htub", "stage3", "stage1")
contrast_htub_3 = c("Stages_htub", "stage4", "stage1")
contrast_htub_4 = c("Stages_htub", "stage3", "stage2")
contrast_htub_5 = c("Stages_htub", "stage4", "stage3")
res_htub_1 = results(ddssva_htub, contrast_htub_1, alpha = 0.05)
res_htub_2 = results(ddssva_htub, contrast_htub_2, alpha = 0.05)
res_htub_3 = results(ddssva_htub, contrast_htub_3, alpha = 0.05)
res_htub_4 = results(ddssva_htub, contrast_htub_4, alpha = 0.05)
res_htub_5 = results(ddssva_htub, contrast_htub_5, alpha = 0.05)
htub_maternal <- unique(c(rownames(subset(res_htub_1, res_htub_1$padj < 0.05 & res_htub_1$log2FoldChange < -2)),
                          rownames(subset(res_htub_2, res_htub_2$padj < 0.05 & res_htub_2$log2FoldChange < -2)),
                          rownames(subset(res_htub_3, res_htub_3$padj < 0.05 & res_htub_3$log2FoldChange < -2)),
                          rownames(subset(res_htub_4, res_htub_4$padj < 0.05 & res_htub_4$log2FoldChange < -2)),
                          rownames(subset(res_htub_5, res_htub_5$padj < 0.05 & res_htub_5$log2FoldChange < -2))))
htub_zygotic <- unique(c(rownames(subset(res_htub_1, res_htub_1$padj < 0.05 & res_htub_1$log2FoldChange >2)),
                         rownames(subset(res_htub_2, res_htub_2$padj < 0.05 & res_htub_2$log2FoldChange >2)),
                         rownames(subset(res_htub_3, res_htub_3$padj < 0.05 & res_htub_3$log2FoldChange >2)),
                         rownames(subset(res_htub_4, res_htub_4$padj < 0.05 & res_htub_4$log2FoldChange >2)),
                         rownames(subset(res_htub_5, res_htub_5$padj < 0.05 & res_htub_5$log2FoldChange >2))))

contrast_lvar_1 = c("Stages_lvar", "stage2", "stage1")
contrast_lvar_2 = c("Stages_lvar", "stage3", "stage1")
contrast_lvar_3 = c("Stages_lvar", "stage4", "stage1")
contrast_lvar_4 = c("Stages_lvar", "stage3", "stage2")
contrast_lvar_5 = c("Stages_lvar", "stage4", "stage3")
res_lvar_1 = results(ddssva_lvar, contrast_lvar_1, alpha = 0.05)
res_lvar_2 = results(ddssva_lvar, contrast_lvar_2, alpha = 0.05)
res_lvar_3 = results(ddssva_lvar, contrast_lvar_3, alpha = 0.05)
res_lvar_4 = results(ddssva_lvar, contrast_lvar_4, alpha = 0.05)
res_lvar_5 = results(ddssva_lvar, contrast_lvar_5, alpha = 0.05)
lvar_maternal <- unique(c(rownames(subset(res_lvar_1, res_lvar_1$padj < 0.05 & res_lvar_1$log2FoldChange < -2)),
                          rownames(subset(res_lvar_2, res_lvar_2$padj < 0.05 & res_lvar_2$log2FoldChange < -2)),
                          rownames(subset(res_lvar_3, res_lvar_3$padj < 0.05 & res_lvar_3$log2FoldChange < -2)),
                          rownames(subset(res_lvar_4, res_lvar_4$padj < 0.05 & res_lvar_4$log2FoldChange < -2)),
                          rownames(subset(res_lvar_5, res_lvar_5$padj < 0.05 & res_lvar_5$log2FoldChange < -2))))
lvar_zygotic <- unique(c(rownames(subset(res_lvar_1, res_lvar_1$padj < 0.05 & res_lvar_1$log2FoldChange >2)),
                         rownames(subset(res_lvar_2, res_lvar_2$padj < 0.05 & res_lvar_2$log2FoldChange >2)),
                         rownames(subset(res_lvar_3, res_lvar_3$padj < 0.05 & res_lvar_3$log2FoldChange >2)),
                         rownames(subset(res_lvar_4, res_lvar_4$padj < 0.05 & res_lvar_4$log2FoldChange >2)),
                         rownames(subset(res_lvar_5, res_lvar_5$padj < 0.05 & res_lvar_5$log2FoldChange >2))))

contrast_mcap_1 = c("Stages_mcap", "stage2", "stage1")
contrast_mcap_2 = c("Stages_mcap", "stage3", "stage1")
contrast_mcap_3 = c("Stages_mcap", "stage4", "stage1")
contrast_mcap_4 = c("Stages_mcap", "stage3", "stage2")
contrast_mcap_5 = c("Stages_mcap", "stage4", "stage3")
res_mcap_1 = results(ddssva_mcap, contrast_mcap_1, alpha = 0.05)
res_mcap_2 = results(ddssva_mcap, contrast_mcap_2, alpha = 0.05)
res_mcap_3 = results(ddssva_mcap, contrast_mcap_3, alpha = 0.05)
res_mcap_4 = results(ddssva_mcap, contrast_mcap_4, alpha = 0.05)
res_mcap_5 = results(ddssva_mcap, contrast_mcap_5, alpha = 0.05)
mcap_maternal <- unique(c(rownames(subset(res_mcap_1, res_mcap_1$padj < 0.05 & res_mcap_1$log2FoldChange < -2)),
                          rownames(subset(res_mcap_2, res_mcap_2$padj < 0.05 & res_mcap_2$log2FoldChange < -2)),
                          rownames(subset(res_mcap_3, res_mcap_3$padj < 0.05 & res_mcap_3$log2FoldChange < -2)),
rownames(subset(res_mcap_4, res_mcap_4$padj < 0.05 & res_mcap_4$log2FoldChange < -2)),
rownames(subset(res_mcap_5, res_mcap_5$padj < 0.05 & res_mcap_5$log2FoldChange < -2))))
mcap_zygotic <- unique(c(rownames(subset(res_mcap_1, res_mcap_1$padj < 0.05 & res_mcap_1$log2FoldChange >2)),
                         rownames(subset(res_mcap_2, res_mcap_2$padj < 0.05 & res_mcap_2$log2FoldChange >2)),
                         rownames(subset(res_mcap_3, res_mcap_3$padj < 0.05 & res_mcap_3$log2FoldChange >2)),
rownames(subset(res_mcap_4, res_mcap_4$padj < 0.05 & res_mcap_4$log2FoldChange >2)),
rownames(subset(res_mcap_5, res_mcap_5$padj < 0.05 & res_mcap_5$log2FoldChange >2))))

contrast_mfra_1 = c("Stages_mfra", "stage2", "stage1")
contrast_mfra_2 = c("Stages_mfra", "stage3", "stage1")
contrast_mfra_3 = c("Stages_mfra", "stage3", "stage2")
res_mfra_1 = results(ddssva_mfra, contrast_mfra_1, alpha = 0.05)
res_mfra_2 = results(ddssva_mfra, contrast_mfra_2, alpha = 0.05)
res_mfra_3 = results(ddssva_mfra, contrast_mfra_3, alpha = 0.05)
mfra_maternal <- unique(c(rownames(subset(res_mfra_1, res_mfra_1$padj < 0.05 & res_mfra_1$log2FoldChange < -2)),
                          rownames(subset(res_mfra_2, res_mfra_2$padj < 0.05 & res_mfra_2$log2FoldChange < -2)),
                          rownames(subset(res_mfra_3, res_mfra_3$padj < 0.05 & res_mfra_3$log2FoldChange < -2))))
mfra_zygotic <- unique(c(rownames(subset(res_mfra_1, res_mfra_1$padj < 0.05 & res_mfra_1$log2FoldChange >2)),
                         rownames(subset(res_mfra_2, res_mfra_2$padj < 0.05 & res_mfra_2$log2FoldChange >2)),
                         rownames(subset(res_mfra_3, res_mfra_3$padj < 0.05 & res_mfra_3$log2FoldChange >2))))

contrast_mlei = c("Stages_mlei", "stage2", "stage1")
res_mlei = results(ddssva_mlei, contrast_mlei, alpha = 0.05)
mlei_maternal <- unique(c(as.character(rownames(subset(res_mlei, padj < 0.05 & log2FoldChange < -2)))))
mlei_zygotic <- unique(c(as.character(rownames(subset(res_mlei, padj < 0.05 & log2FoldChange >2)))))

contrast_pdum_1 = c("Stages_pdum", "stage2", "stage1")
res_pdum_1 = results(ddssva_pdum, contrast_pdum_1, alpha = 0.05)
res_pdum_maternal_1 <- subset(res_pdum_1, padj < 0.05 & log2FoldChange < -2)
pdum_maternal <- unique(c(rownames(res_pdum_maternal_1)))
pdum_zygotic <- unique(c(rownames(subset(res_pdum_1, padj < 0.05 & log2FoldChange >2))))

contrast_pliv_1 = c("Stages_pliv", "stage2", "stage1")
contrast_pliv_2 = c("Stages_pliv", "stage3", "stage1")
contrast_pliv_3 = c("Stages_pliv", "stage3", "stage2")
res_pliv_1 = results(ddssva_pliv, contrast_pliv_1, alpha = 0.05)
res_pliv_2 = results(ddssva_pliv, contrast_pliv_2, alpha = 0.05)
res_pliv_3 = results(ddssva_pliv, contrast_pliv_3, alpha = 0.05)
pliv_maternal <- unique(c(rownames(subset(res_pliv_1, res_pliv_1$padj < 0.05 & res_pliv_1$log2FoldChange < -2)),
                          rownames(subset(res_pliv_2, res_pliv_2$padj < 0.05 & res_pliv_2$log2FoldChange < -2)),
                          rownames(subset(res_pliv_3, res_pliv_3$padj < 0.05 & res_pliv_3$log2FoldChange < -2))))
pliv_zygotic <- unique(c(rownames(subset(res_pliv_1, res_pliv_1$padj < 0.05 & res_pliv_1$log2FoldChange >2)),
                         rownames(subset(res_pliv_2, res_pliv_2$padj < 0.05 & res_pliv_2$log2FoldChange >2)),
                         rownames(subset(res_pliv_3, res_pliv_3$padj < 0.05 & res_pliv_3$log2FoldChange >2))))

contrast_pmin_1 = c("Stages_pmin", "stage2", "stage1")
contrast_pmin_2 = c("Stages_pmin", "stage3", "stage1")
contrast_pmin_3 = c("Stages_pmin", "stage3", "stage2")
res_pmin_1 = results(ddssva_pmin, contrast_pmin_1, alpha = 0.05)
res_pmin_2 = results(ddssva_pmin, contrast_pmin_2, alpha = 0.05)
res_pmin_3 = results(ddssva_pmin, contrast_pmin_3, alpha = 0.05)
pmin_maternal <- unique(c(rownames(subset(res_pmin_1, res_pmin_1$padj < 0.05 & res_pmin_1$log2FoldChange < -2)),
                          rownames(subset(res_pmin_2, res_pmin_2$padj < 0.05 & res_pmin_2$log2FoldChange < -2)),
                          rownames(subset(res_pmin_3, res_pmin_3$padj < 0.05 & res_pmin_3$log2FoldChange < -2))))
pmin_zygotic <- unique(c(rownames(subset(res_pmin_1, res_pmin_1$padj < 0.05 & res_pmin_1$log2FoldChange >2)),
                         rownames(subset(res_pmin_2, res_pmin_2$padj < 0.05 & res_pmin_2$log2FoldChange >2)),
                         rownames(subset(res_pmin_3, res_pmin_3$padj < 0.05 & res_pmin_3$log2FoldChange >2))))

contrast_sscr_1 = c("Stages_sscr", "stage2", "stage1")
contrast_sscr_2 = c("Stages_sscr", "stage3", "stage1")
contrast_sscr_3 = c("Stages_sscr", "stage3", "stage2")
res_sscr_1 = results(ddssva_sscr, contrast_sscr_1, alpha = 0.05)
res_sscr_2 = results(ddssva_sscr, contrast_sscr_2, alpha = 0.05)
res_sscr_3 = results(ddssva_sscr, contrast_sscr_3, alpha = 0.05)
sscr_maternal <- unique(c(rownames(subset(res_sscr_1, res_sscr_1$padj < 0.05 & res_sscr_1$log2FoldChange < -2)),
                          rownames(subset(res_sscr_2, res_sscr_2$padj < 0.05 & res_sscr_2$log2FoldChange < -2)),
                          rownames(subset(res_sscr_3, res_sscr_3$padj < 0.05 & res_sscr_3$log2FoldChange < -2))))
sscr_zygotic <- unique(c(rownames(subset(res_sscr_1, res_sscr_1$padj < 0.05 & res_sscr_1$log2FoldChange >2)),
                         rownames(subset(res_sscr_2, res_sscr_2$padj < 0.05 & res_sscr_2$log2FoldChange >2)),
                         rownames(subset(res_sscr_3, res_sscr_3$padj < 0.05 & res_sscr_3$log2FoldChange >2))))
sscr_maternal <- str_remove_all(sscr_maternal, "\\.[0-9]+$")
sscr_zygotic <- str_remove_all(sscr_zygotic, "\\.[0-9]+$")


contrast_ttra_1 = c("Stages_ttra", "stage2", "stage1")
contrast_ttra_2 = c("Stages_ttra", "stage3", "stage1")
contrast_ttra_3 = c("Stages_ttra", "stage3", "stage2")
res_ttra_1 = results(ddssva_ttra, contrast_ttra_1, alpha = 0.05)
res_ttra_2 = results(ddssva_ttra, contrast_ttra_2, alpha = 0.05)
res_ttra_3 = results(ddssva_ttra, contrast_ttra_3, alpha = 0.05)
ttra_maternal <- unique(c(rownames(subset(res_ttra_1, res_ttra_1$padj < 0.05 & res_ttra_1$log2FoldChange < -2)),
                          rownames(subset(res_ttra_2, res_ttra_2$padj < 0.05 & res_ttra_2$log2FoldChange < -2)),
                          rownames(subset(res_ttra_3, res_ttra_3$padj < 0.05 & res_ttra_3$log2FoldChange < -2))))
ttra_zygotic <- unique(c(rownames(subset(res_ttra_1, res_ttra_1$padj < 0.05 & res_ttra_1$log2FoldChange >2)),
                         rownames(subset(res_ttra_2, res_ttra_2$padj < 0.05 & res_ttra_2$log2FoldChange >2)),
                         rownames(subset(res_ttra_3, res_ttra_3$padj < 0.05 & res_ttra_3$log2FoldChange >2))))

##### Results checking ######
#Annotation table from FLY-FISH (http://fly-fish.ccbr.utoronto.ca/exports/)
annot <- read.csv("~/Documents/Gene_expr_evol/Intermediate_files/annotation.csv")
annot$gene <- str_remove_all(annot$gene, " .*")
annot$gene <- str_remove_all(annot$gene, "\\?")

annot$asd <- tibble(annot) %>% 
  mutate(conv = mapIds(org.Dm.eg.db, annot$gene, 'SYMBOL', 'FLYBASECG')) %>%
  mutate(gene = coalesce(conv, gene)) %>%
  dplyr::select(-conv)
    
library(org.Dm.eg.db)
keytypes(org.Dm.eg.db)

dmel_maternal_x <- mapIds(org.Dm.eg.db, dmel_maternal, 'SYMBOL', 'FLYBASE')
dmel_zygotic_x <- mapIds(org.Dm.eg.db, dmel_zygotic, 'SYMBOL', 'FLYBASE')

#How many differentialy expressed genes are present in the database?
summary(dmel_maternal_x %in% annot$gene)
summary(annot$gene %in% dmel_maternal_x)

#How many in the database do have maternal expression and are not found in DGE?
test <- subset(annot, !(annot$gene %in% dmel_maternal_x))
test <- subset(test, test$term == "Maternal" & test$stage.tissue == "stage 1-3")
notfound <- length(unique(test$gene))

test <- subset(annot, annot$gene %in% dmel_maternal_x)
test <- subset(test, test$term == "Maternal" & test$stage.tissue == "stage 1-3")
found <- length(unique(test$gene))

total <- subset(annot, annot$term == "Maternal" & test$stage.tissue == "stage 1-3")
total <- subset(total, str_detect(total$term, "Degrad"))
total <- length(unique(total$gene))

#Subset for genes present in the database
test <- subset(annot, annot$gene %in% dmel_maternal_x)
#How many genes have Maternal/stage 1-3 annotatio?
withexp <- length(unique(subset(test, test$term == "Maternal" | test$stage.tissue == "stage 1-3")[, 1]))
withoutexp <- length(test[!(test$Probe %in% unique(subset(test, test$term == "Maternal" | test$stage.tissue == "stage 1-3")[, 1])), ])
#What about the missing ones?
test[!(test$Probe %in% unique(subset(test, test$term == "Maternal" | test$stage.tissue == "stage 1-3")[, 1])), ]


##### SAVING #####
#save.image(file='/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files/DGE_script_enviorment.RData')
#load('~/Documents/Gene_expr_evol/Intermediate_files/DGE_script_enviorment.RData')

mat_ids <- readRDS("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files/maternal_IDs.RDS")
saveRDS(deg_ids, "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files/downregulated_IDs_NFE-T.RDS")

#Loop export
for(i in 1:length(variables)){
  temp <- get(variables[i])
  
  saveRDS(temp, paste("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files", paste0(variables[i], "_NFE-T.RDS"), sep = "/"))
  
}
