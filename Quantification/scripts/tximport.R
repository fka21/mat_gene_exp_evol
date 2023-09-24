
# Script for importing quantification files and defining maternal genes based on TPM

##### LOAD PACKAGES #####
Libraries <- c("DESeq2", "pheatmap", "vsn", "hexbin", "cowplot", "apeglm", "tximport", 
               "viridis", "PoiClaClu", "genefilter", "vidger", "reshape2", "ggplot2", "stringr", 
               "gplots", "sva", "tidyverse")
lapply(Libraries, library, character.only = TRUE)

##### FUNCTION #####
source("~/Documents/Gene_expr_evol/Scripts/functions.R")

##### TXIMPORT #####
#Initialize vector for all maternal gene IDs
mat_ids_final <- c()

#Manually importing every quantification file for each species and formatting developmental stagings for later

files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/A_stephensi",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/A_stephensi", files_ext, "quant.sf")
names(files) <- c("stage3_1", "stage4_1", "stage1_1", "stage2_1", "stage2_2", "stage1_2", "stage4_2", "stage3_2", "stage2_3")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/A_stephensi/tx2gene.tsv", header = F)
txi_aste <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/A_suum",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/A_suum", files_ext, "quant.sf")
names(files) <- c("stage1_1", "stage1_2", "stage1_3", "stage2_1", "stage2_2", "stage2_3", "stage2_4", "stage3_1", "stage3_2", "stage3_3",
                  "stage4_1", "stage4_2")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/A_suum/tx2gene.tsv", header = F)
txi_asum <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_germanica",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_germanica", files_ext, "quant.sf")
names(files) <- c("stage1_1", "stage1_2", "stage2_1", "stage2_2", "stage3_1", "stage3_2", "stage4_1", "stage4_2")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_germanica/tx2gene.tsv", header = F)
txi_bger <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_lancelatum",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_lancelatum", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage2"), 1:5, sep = "_"),
                  paste(rep("stage3"), 1:2, sep = "_"),
                  paste(rep("stage1"), 1:2, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_lancelatum/tx2gene.tsv", header = F)
txi_blan <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_mori",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_mori", files_ext, "quant.sf")
names(files) <- c("stage2_1", "stage2_2", "stage3_1", "stage3_2", "stage4_1", "stage4_2", "stage5_1", "stage5_2", "stage1_1", "stage1_2", "stage1_3")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_mori/tx2gene.tsv", header = F)
txi_bmor <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_taurus",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_taurus", files_ext, "quant.sf")
names(files) <- c("stage1_1", "stage2_1", "stage2_2", "stage1_2", "stage2_3", "stage1_3", "stage1_4")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/B_taurus/tx2gene.tsv", header = F)
txi_btau <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_angaria_2017",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_angaria_2017", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:4, sep = "_"),
                  paste(rep("stage2"), 1:4, sep = "_"),
                  paste(rep("stage3"), 1:4, sep = "_"),
                  paste(rep("stage4"), 1:4, sep = "_"),
                  paste(rep("stage5"), 1:4, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_angaria_2017/tx2gene.tsv", header = F)
txi_cana <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_elegans_2017",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_elegans_2017", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:4, sep = "_"),
                  paste(rep("stage2"), 1:4, sep = "_"),
                  paste(rep("stage3"), 1:4, sep = "_"),
                  paste(rep("stage4"), 1:4, sep = "_"),
                  paste(rep("stage5"), 1:4, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_elegans_2017/tx2gene.tsv", header = F)
txi_cele <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_gigas",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_gigas", files_ext, "quant.sf")
names(files) <- c("stage2_1", "stage2_2", "stage1_1", "stage1_2", "stage2_3", "stage1_3", "stage1_4", "stage2_4")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_gigas/tx2gene.tsv", header = F)
txi_cgig <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_hircus",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_hircus", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:3, sep = "_"), paste(rep("stage2"), 1:3, sep = "_"), paste(rep("stage3"), 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_hircus/tx2gene.tsv", header = F)
txi_chir <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_intestinalis",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_intestinalis", files_ext, "quant.sf")
names(files) <- c("stage2_1", "stage3_1", "stage1_1", "stage2_2", "stage3_2", "stage1_2")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/C_intestinalis/tx2gene.tsv", header = F)
txi_cint <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_ananassae_2018",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_ananassae_2018", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:3, sep = "_"), paste(rep("stage2"), 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_ananassae_2018/tx2gene.tsv", header = F)
txi_dana <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_erecta_2018",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_erecta_2018", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:3, sep = "_"), paste(rep("stage2"), 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_erecta_2018/tx2gene.tsv", header = F)
txi_dere <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_melanogaster_2018",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_melanogaster_2018", files_ext, "quant.sf")
names(files) <- c("stage1_1", "stage2_1", "stage1_2", "stage3_1", "stage1_3", "stage3_2", "stage2_2", "stage2_3", "stage2_4", "stage3_3", "stage1_4", "stage1_5", "stage3_4", "stage3_5")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_melanogaster_2018/tx2gene.tsv", header = F)
txi_dmel <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))

files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_mojavensis_2018",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_mojavensis_2018", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:6, sep = "_"), paste(rep("stage2"), 1:5, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_mojavensis_2018/tx2gene.tsv", header = F)
txi_dmoj <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_persimilis_2018",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_persimilis_2018", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:6, sep = "_"), paste(rep("stage2"), 1:5, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_persimilis_2018/tx2gene.tsv", header = F)
txi_dper <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_rerio",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_rerio", files_ext, "quant.sf")
names(files) <- c(paste("stage2", 1:3, sep = "_"), paste("stage5", 1:3, sep = "_"), paste("stage4", 1, sep = "_"), paste("stage1", 1, sep = "_"),
                  paste("stage4", 2:3, sep = "_"), paste("stage3", 1:3, sep = "_"), paste("stage1", 2:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_rerio/tx2gene.tsv", header = F)
txi_drer <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_simulans_2018",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_simulans_2018", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:3, sep = "_"), paste("stage2", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_simulans_2018/tx2gene.tsv", header = F)
txi_dsim <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_virilis_2018",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_virilis_2018", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:7, sep = "_"), paste("stage2", 1:6, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_virilis_2018/tx2gene.tsv", header = F)
txi_dvir <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_willistoni_2018",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_willistoni_2018", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:6, sep = "_"), paste("stage2", 1:6, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_willistoni_2018/tx2gene.tsv", header = F)
txi_dwil <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_yakuba",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_yakuba", files_ext, "quant.sf")
names(files) <- c(paste("stage2", 1:6, sep = "_"), paste("stage1", 1:2, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/D_yakuba/tx2gene.tsv", header = F)
txi_dyak <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))

files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/G_gallus",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/G_gallus", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:3, sep = "_"), paste("stage2", 1:3, sep = "_"), paste("stage3", 1:3, sep = "_"), 
                  paste("stage4", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/G_gallus/tx2gene.tsv", header = F)
txi_ggal <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_dujardini",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_dujardini", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:4, sep = "_"), paste("stage2", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_dujardini/tx2gene.tsv", header = F)
txi_hduj <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_erythrogramma",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_erythrogramma", files_ext, "quant.sf")
names(files) <- c(paste("stage3", 1:3, sep = "_"), paste("stage4", 1:3, sep = "_"), paste("stage2", 1:3, sep = "_"),
                  paste("stage1", 1:3, sep = "_"), paste("stage6", 1:3, sep = "_"), paste("stage5", 1:3, sep = "_"),
                  paste("stage7", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_erythrogramma/tx2gene.tsv", header = F)
txi_hery <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_sapiens",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_sapiens", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:2, sep = "_"), paste("stage2", 1:3, sep = "_"), paste("stage3", 1:4, sep = "_"),
                  paste("stage4", 1:11, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_sapiens/tx2gene.tsv", header = F)
txi_hsap <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_tuberculata",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_tuberculata", files_ext, "quant.sf")
names(files) <- c(paste("stage3", 1:3, sep = "_"), paste("stage4", 1:3, sep = "_"), paste("stage2", 1:3, sep = "_"),
                  paste("stage1", 1:3, sep = "_"), paste("stage6", 1:3, sep = "_"), paste("stage5", 1:3, sep = "_"),
                  paste("stage7", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/H_tuberculata/tx2gene.tsv", header = F)
txi_htub <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/L_variegatus",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/L_variegatus", files_ext, "quant.sf")
names(files) <- c(paste("stage3", 1:3, sep = "_"), paste("stage4", 1:3, sep = "_"), paste("stage2", 1:3, sep = "_"),
                  paste("stage1", 1:3, sep = "_"), paste("stage6", 1:3, sep = "_"), paste("stage5", 1:3, sep = "_"),
                  paste("stage7", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/L_variegatus/tx2gene.tsv", header = F)
txi_lvar <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_capitata",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_capitata", files_ext, "quant.sf")
names(files) <- c(paste("stage2", 1, sep = "_"), paste("stage5", 1:2, sep = "_"), paste("stage4", 1:2, sep = "_"),
                  paste("stage6", 1:2, sep = "_"), paste("stage2", 2, sep = "_"), paste("stage6", 3:6, sep = "_"),
                  paste("stage7", 1:4, sep = "_"), paste("stage1", 1, sep = "_"), paste("stage7", 5:6, sep = "_"),
                  paste("stage3", 1:3, sep = "_"), paste("stage1", 2:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_capitata/tx2gene.tsv", header = F)
txi_mcap <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_franciscanus",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_franciscanus", files_ext, "quant.sf")
names(files) <- c("stage4_1", "stage4_2", "stage1_1", "stage1_2", "stage2_1", "stage1_3", "stage2_2", "stage2_3", "stage3_1", "stage3_2", "stage4_3", "stage3_3")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_franciscanus/tx2gene.tsv", header = F)
txi_mfra <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_leidyi",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_leidyi", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:4, sep = "_"), paste("stage2", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_leidyi/tx2gene.tsv", header = F)
txi_mlei <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_mulatta",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_mulatta", files_ext, "quant.sf")
names(files) <- c("stage2_1", "stage4_1", paste("stage1", 1:3, sep = "_"), paste("stage2", 2:4, sep = "_"), paste("stage3", 1:5, sep = "_"), paste("stage4", 2:4, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_mulatta/tx2gene.tsv", header = F)
txi_mmul <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_musculus",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_musculus", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:3, sep = "_"), paste("stage2", 1:3, sep = "_"), paste("stage3", 1:3, sep = "_"),
                  paste("stage4", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/M_musculus/tx2gene.tsv", header = F)
txi_mmus <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))

files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/N_vectensis",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/N_vectensis", files_ext, "quant.sf")
names(files) <- c("stage1_1", "stage2_1", "stage3_1", "stage1_2", "stage2_2", "stage3_2")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/N_vectensis/tx2gene.tsv", header = F)
txi_nvec <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_caudatus",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_caudatus", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:2, sep = "_"), paste("stage2", 1:2, sep = "_"), paste("stage3", 1:2, sep = "_"),
                  paste("stage4", 1:2, sep = "_"), paste("stage5", 1:2, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_caudatus/tx2gene.tsv", header = F)
txi_pcau <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_dumerili",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_dumerili", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:7, sep = "_"), paste("stage2", 1:7, sep = "_"), paste("stage3", 1:7, sep = "_"),
                  paste("stage4", 1:7, sep = "_"), paste("stage5", 1:29, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_dumerili/tx2gene.tsv", header = F)
txi_pdum <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_lividus",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_lividus", files_ext, "quant.sf")
names(files) <- c(paste("stage8", 1:3, sep = "_"), paste("stage7", 1:3, sep = "_"), paste("stage6", 1:3, sep = "_"),
                  paste("stage5", 1:3, sep = "_"), paste("stage4", 1:3, sep = "_"), paste("stage3", 1:3, sep = "_"),
                  paste("stage2", 1:3, sep = "_"), paste("stage1", 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_lividus/tx2gene.tsv", header = F)
txi_pliv <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))

files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_miniata",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_miniata", files_ext, "quant.sf")
names(files) <- c(paste("stage2", 1:2, sep = "_"), paste("stage3", 1, sep = "_"), paste("stage1", 1:3, sep = "_"),
                  paste("stage4", 1:2, sep = "_"), paste("stage3", 2:3, sep = "_"), paste("stage4", 3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/P_miniata/tx2gene.tsv", header = F)
txi_pmin <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_carpocapsae_2017",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_carpocapsae_2017", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:4, sep = "_"),
                  paste(rep("stage2"), 1:4, sep = "_"),
                  paste(rep("stage3"), 1:4, sep = "_"),
                  paste(rep("stage4"), 1:4, sep = "_"),
                  paste(rep("stage5"), 1:4, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_carpocapsae_2017/tx2gene.tsv", header = F)
txi_scar <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_feltiae_2017",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_feltiae_2017", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:4, sep = "_"),
                  paste(rep("stage2"), 1:4, sep = "_"),
                  paste(rep("stage3"), 1:4, sep = "_"),
                  paste(rep("stage4"), 1:4, sep = "_"),
                  paste(rep("stage5"), 1:4, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_feltiae_2017/tx2gene.tsv", header = F)
txi_sfel <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_scrofa",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_scrofa", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage2"), 1:6, sep = "_"),
                  paste(rep("stage3"), 1:13, sep = "_"),
                  paste(rep("stage4"), 1:23, sep = "_"),
                  paste(rep("stage1"), 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/S_scrofa/tx2gene.tsv", header = F)
txi_sscr <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))



files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/T_transversa",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/T_transversa", files_ext, "quant.sf")
names(files) <- names(files) <- c("stage1_1", "stage2_1", "stage3_1", "stage4_1", "stage5_1", "stage6_1", "stage8_1", "stage10_1", "stage11_1",
                                  "stage12_1", "stage13_1", "stage14_1", "stage1_2", "stage2_2", "stage3_2", "stage4_2", "stage6_2", "stage9_1",
                                  "stage10_2", "stage11_2", "stage14_2", "stage7_1", "stage9_2", "stage5_2", "stage7_2", "stage8_2", "stage12_2", "stage13_2")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/T_transversa/tx2gene.tsv", header = F)
txi_ttra <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/T_castenum",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/T_castenum", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1"), 1:2, sep = "_"),
                  paste(rep("stage4"), 1:2, sep = "_"),
                  paste(rep("stage5"), 1:2, sep = "_"),
                  paste(rep("stage2"), 1:12, sep = "_"),
                  paste(rep("stage3"), 1:12, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/T_castenum/tx2gene.tsv", header = F)
txi_tcas <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


files_ext <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/X_tropicalis",  pattern = "*salmon")
files <- file.path("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/X_tropicalis", files_ext, "quant.sf")
names(files) <- c("stage1_1", "stage2_1", "stage3_1", "stage4_1", "stage5_1", "stage6_1",
                  "stage1_3", "stage2_2", "stage3_2", "stage4_2", "stage5_2", 
                  "stage1_2", "stage2_3", "stage3_3", "stage4_3", "stage5_3", "stage6_3")
all(file.exists(files))
tx2gene <- read.table("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/X_tropicalis/tx2gene.tsv", header = F)
txi_xtro <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

#Export abundance tables
export.txtable(ls(pattern = "txi_...."))
#Extracting IDs which have TPM > 2 in oocyte stages
mat_ids_final <- c(mat_ids_final, expr.extr(ls(pattern = "txi_....")))
#Working with big files, to save up space
rm(list = ls(pattern = "txi_...."))


##### SALMON EXPORT #####
saveRDS(mat_ids_final, "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files/maternal_IDs.RDS")

