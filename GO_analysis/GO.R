###### LIBRARIES #####
library(DESeq2)
library(clusterProfiler)
library(tidyverse)
library(ggplotify)
library(viridis)

####### READ IN DATA #######
load("~/Documents/Gene_expr_evol/Intermediate_files/DGE_script_enviorment.RData") #DGE analysis output
rm(list = setdiff(ls(), ls(pattern = "res_*"))) #removing parts to save memory

#Read in GO annotations
variables <- list.files("~/Documents/Gene_expr_evol/GO/Annotations/", pattern = "*go.tsv",
                        full.names = T)
#Read in maternal and degraded IDs for subsetting later on
mat_ids <- readRDS("~/Documents/Gene_expr_evol/Intermediate_files/maternal_IDs.RDS")
deg_ids <- readRDS("~/Documents/Gene_expr_evol/Intermediate_files/downregulated_IDs_NFE-T.RDS")

#Crosreference table for later loops
orgdbs <- data.frame(code = c("hsap", "drer", "xtro", "ggal", "mmul", "mmus", "cele", "dmel", "sscr", "btau"),
                    orgdb = c("org.Hs.eg.db", "org.Dr.eg.db", "org.Xl.eg.db", "org.Gg.eg.db", "org.Mmu.eg.db",
                              "org.Mm.eg.db", "org.Ce.eg.db", "org.Dm.eg.db", "org.Ss.eg.db", "org.Bt.eg.db"))
id_crosref <- read.table("~/Documents/Gene_expr_evol/GO/Annotations/uniprot_go_crosref.tsv",
                         header = F, sep = "\t")

#Initializing empty lists
dge_list <- vector(mode = "list", length = length(variables))
mat_list <- vector(mode = "list", length = length(variables))
df <- data.frame()

##### GO ENRICHMENTS #######

for(i in 1:length(variables)){
  temp_name <- str_remove_all(basename(variables[i]), "_go.tsv$")
  
  #Reading in data and formatting it
  temp <- read_tsv(variables[i]) %>% filter(ARGOT_PPV >= 0.5)
  temp$goid <- paste0("GO:", temp$goid)
  temp <- temp %>% dplyr::select(qpid, goid, desc)
  
  #Adjusting IDs to be able to overlap protein IDs and gene IDs
  if(temp_name == "bger"){
    crosref <- read.csv("~/Documents/Gene_expr_evol/GO/bger_prot2geneid.csv", header = T)
    crosref$Locus.tag[!(crosref$Locus == "")] <- crosref$Locus[!(crosref$Locus == "")]
    temp$qpid <- crosref$Locus.tag[match(temp$qpid, crosref$Protein.product)]
    rm(crosref)
  } else if(temp_name == "pmin"){
    crosref <- read.csv("~/Documents/Gene_expr_evol/GO/pmin_prot2geneid.csv", header = T)
    temp$qpid <- crosref$Locus[match(temp$qpid, crosref$Protein.product)]
    rm(crosref)
  } else if (sum(str_detect(temp$qpid, "t[0-9]+$")) > 0 & !(temp_name == "sfel")){
    temp$qpid <- str_remove_all(temp$qpid, "t[0-9]+$")
  } else if(!(temp_name == "sfel")){
    temp$qpid <- str_remove_all(temp$qpid, "\\.[0-9]+$")
  }
  
  if(temp_name == "mcap" | temp_name == "cana"){
    temp$qpid <- str_remove_all(temp$qpid, "\\.$")
  } else if(temp_name == "asum"){
    temp$qpid <- str_remove_all(temp$qpid, "_$")
  }
  

  #Assigning table for enricher use
  assign(paste0(temp_name, "_annot"), temp)
  
  #Creating gene list for each species first for degraded
  temp_var <- ls(pattern = paste0("res_", (temp_name)))
  temp_glst <- data.frame()
  for(k in 1:length(temp_var)){
    int_df <- as.data.frame(get(temp_var[k]))
    temp_glst <- rbind(temp_glst, int_df)
  }
  
  temp_glst$id <- rownames(temp_glst)
  #Adjusting ID
  if(!(temp_name == "sfel")){
    temp_glst$id <- str_remove_all(temp_glst$id, "\\.[0-9]+$") 
  }
  #Building table for log2FoldChange plots
  logfold <- temp_glst; logfold$sp <- temp_name
  
  #Keeping only degraded ones
  temp_glst <- subset(temp_glst, temp_glst$log2FoldChange <= -2 & temp_glst$padj <= 0.05)
  temp_glst <- tibble(temp_glst) %>% group_by(id) %>% slice(which.max(abs(log2FoldChange))) %>%
    dplyr::select(log2FoldChange, id) %>% arrange(log2FoldChange)
  
  
  int_df <- abs(temp_glst$log2FoldChange); names(int_df) <- temp_glst$id
  
  assign(paste0(temp_name, "_dge_lst"), int_df)
  
  #Creating list for maternal genes
  pth <- paste0("~/Documents/Gene_expr_evol/Intermediate_files/",
                paste0("txi_", paste0(temp_name, ".RDS")))
  txi <- readRDS(pth)
  txi <- as.data.frame(t(txi$abundance)) #Extracting and transforming abundance table
  txi_metadata <- rownames(txi) #Get metadata for counts
  txi$replicate <- as.factor(str_extract(txi_metadata, "[0-9]+$")) #Add replicate metadata to tbl
  txi$stage <- as.factor(str_remove_all(txi_metadata, "_[0-9]+")) #Add stage metadata to tbl
  
  txi <- tibble(txi) %>%
    filter(txi$stage == "stage1") %>%
    dplyr::select(-stage, -replicate)
  
  txi <- as.data.frame(t(txi[, -1]))
  txi <- subset(txi, rowMeans(txi) >= 2) #using maternal genes which are not in degraded category
  txi <- txi[order(rowMeans(txi), decreasing = T), ]
  if(!(temp_name == "sfel")){
    rownames(txi) <- str_remove_all(rownames(txi), "\\.[0-9]+$") 
  }
  
  int_df <- rowMeans(txi); names(int_df) <- rownames(txi)
  
  #Add TPM values for logfold table and format that one
  logfold <- tibble(logfold)  %>%
    distinct( .keep_all=T) %>%
    mutate(category = case_when(id %in% deg_ids ~ "degraded",
                                id %in% mat_ids & !(id %in% deg_ids) ~ "maternal")) %>%
    dplyr::select(log2FoldChange, sp, category, padj, id)
  logfold$tpm <- int_df[match(logfold$id, names(int_df))]
  df <- rbind(df, logfold)
  
  assign(paste0(temp_name, "_mat_lst"), int_df)
 
  #Enrichment
  #Enrichment where information is available
  genelist_dge <- get(paste0(temp_name, "_dge_lst"))
  genelist_mat <- get(paste0(temp_name, "_mat_lst"))
  if(temp_name %in% c("hsap", "drer", "ggal", "mmus", "cele", "dmel", "btau")){
    library(orgdbs$orgdb[orgdbs$code %in% temp_name], character.only=TRUE)
    
    #Retrieve entrezIDs and build genelist for degraded
    entrez <- select(get(orgdbs$orgdb[orgdbs$code %in% temp_name]), keys=names(genelist_dge), columns="ENTREZID", keytype="ENSEMBL")
    names(genelist_dge) <- entrez$ENTREZID[match(names(genelist_dge), entrez$ENSEMBL)]
    
    #Enrichment and filtering for qValue
    temp_enrich_dge <- enrichGO(gene = names(genelist_dge), OrgDb = get(orgdbs$orgdb[orgdbs$code %in% temp_name]), 
                         ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
    
    #Retrieve entrezIDs and build genelist for maternal
    entrez <- select(get(orgdbs$orgdb[orgdbs$code %in% temp_name]), keys=names(genelist_mat), columns="ENTREZID", keytype="ENSEMBL")
    names(genelist_mat) <- entrez$ENTREZID[match(names(genelist_mat), entrez$ENSEMBL)]
    
    #Enrichment and filtering for qValue
    temp_enrich_mat <- enrichGO(gene = names(genelist_mat), OrgDb = get(orgdbs$orgdb[orgdbs$code %in% temp_name]), 
                                ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
    
    unloadNamespace(orgdbs$orgdb[orgdbs$code %in% temp_name])
    rm(entrez)
    #Enrichment for uniprot annotated transcriptomes
  } else if(temp_name %in% c("hduj", "hery", "htub", "lvar", "pdum", "pliv", "ttra")){
    crosref <- read.table(paste0("~/Documents/Gene_expr_evol/GO/Annotations/", paste0(temp_name, "_crosref.tsv")), header = F, sep = " ")
    temp_crosref <- id_crosref[id_crosref$V1 %in% crosref$V2, ]
    temp_crosref$V1 <- crosref$V1[match(temp_crosref$V1, crosref$V2)]
    
    temp_enrich_dge <- enricher(names(genelist_dge),
                           TERM2GENE = temp_crosref[,c(2,1)], TERM2NAME = temp_crosref[, c(2,3)])
    temp_enrich_mat <- enricher(names(genelist_mat),
                            TERM2GENE = temp_crosref[ ,c(2,1)], TERM2NAME = temp_crosref[ ,c(2,3)])
    rm(crosref)
  } else {
    temp_enrich_dge <-enricher(names(genelist_dge),
                               TERM2GENE = temp[,c(2,1)], TERM2NAME = temp[, c(2,3)])
    
    temp_enrich_mat <- enricher(names(genelist_mat),
                                TERM2GENE = temp[ ,c(2,1)], TERM2NAME = temp[ ,c(2,3)])
  }
  
  #Saving enrichments
  if(dim(temp_enrich_dge@result[temp_enrich_dge@result$p.adjust <= 0.05, ])[1] > 0){
    assign(paste0(temp_name, "_dge_enrich"), temp_enrich_dge)
    dge_list[[i]] <- temp_enrich_dge
    names(dge_list)[i] <- temp_name
  }

  if(dim(temp_enrich_mat@result[temp_enrich_mat@result$p.adjust <= 0.05, ])[1] > 0){
    assign(paste0(temp_name, "_mat_enrich"), temp_enrich_mat)
    mat_list[[i]] <- temp_enrich_mat
    names(mat_list)[i] <- temp_name
  }
   
}

###### GO PLOT ######
#Bring together species specific enrichment results
dge_list_final <- merge_result(dge_list)
mat_list_final <- merge_result(mat_list)

#Adjusting factor levels to species tree later in plotting
dge_list_final@compareClusterResult$Cluster <- factor(dge_list_final@compareClusterResult$Cluster, levels = c("dere", "dyak", "dmel", "dsim", "dana", "dper", "dwil", "dmoj", "dvir",
                         "aste", "bmor", "tcas", "bger", "hduj", "cana", "cele", "scar", "sfel",
                         "asum", "pcau", "cgig", "ttra", "pdum", "btau", "chir", "sscr", "hsap",
                         "mmul", "mmus", "ggal", "xtro", "drer", "cint", "blan", "hery", "htub",
                         "lvar", "mfra", "pliv", "pmin", "mcap", "nvec", "mlei"))
mat_list_final@compareClusterResult$Cluster <- factor(mat_list_final@compareClusterResult$Cluster, levels = c("dere", "dyak", "dmel", "dsim", "dana", "dper", "dwil", "dmoj", "dvir",
                                                                  "aste", "bmor", "tcas", "bger", "hduj", "cana", "cele", "scar", "sfel",
                                                                  "asum", "pcau", "cgig", "ttra", "pdum", "btau", "chir", "sscr", "hsap",
                                                                  "mmul", "mmus", "ggal", "xtro", "drer", "cint", "blan", "hery", "htub",
                                                                  "lvar", "mfra", "pliv", "pmin", "mcap", "nvec", "mlei"))


#Subsetting for significant enrichments
data_plot <- mat_list_final@compareClusterResult[mat_list_final@compareClusterResult$p.adjust <= 0.05, ]
data_plot <- dge_list_final@compareClusterResult[dge_list_final@compareClusterResult$p.adjust <= 0.05, ]

#Prepare data for plotting
data_plot <- data_plot %>% tibble() %>% mutate(ratio1 = as.numeric(sapply(strsplit(GeneRatio,"/"), `[`, 1))) %>%
  mutate(ratio2 = as.numeric(sapply(strsplit(GeneRatio,"/"), `[`, 2))) %>%
  mutate(GeneRatio = ratio1/ratio2) %>%
  dplyr::select(-ONTOLOGY) %>%
  dplyr::group_by(Description) %>%
  mutate(freq = n()) %>% 
  ungroup() %>%
  dplyr::group_by(Cluster) %>%
  mutate(Size = n())
  

data_plot$Description <- factor(data_plot$Description, 
                                levels = unique(data_plot$Description[order(data_plot$freq)]), 
                                ordered = T)

#Subsetting for top enriched terms for better visualizations
data_plot <- data_plot %>% group_by(freq) %>% top_n(n = 5, wt = Description)

##### PLOTTING GO ENRICHMENT ######

#Enrichment plots
svg("~/Desktop/Figures/go_mat.svg", width = 14, height = 9.09)
print(ggplot(data_plot, aes(x = Cluster, y = Description, color = p.adjust, size = GeneRatio)) +
  geom_point()+
  scale_size_continuous(range = c(2, 6)) +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  xlab("") +
  ylab("") +
  viridis::scale_color_viridis() +
  theme_bw() +
  theme(text = element_text(size = 15),
        plot.margin=unit(c(0.3,0.1,0.1,0),"cm")))
dev.off()
ggsave(file="~/Desktop/Figures/go_mat.svg", width = 14, height = 9.09)


###### LOG2FOLDCHANGE PLOT ######
#Adjusting factor levels
df$sp <- factor(df$sp, levels = c("dere", "dyak", "dmel", "dsim", "dana", "dper", "dwil", "dmoj", "dvir",
                                  "aste", "bmor", "tcas", "bger", "hduj", "cana", "cele", "scar", "sfel",
                                  "asum", "pcau", "cgig", "ttra", "pdum", "btau", "chir", "sscr", "hsap",
                                  "mmul", "mmus", "ggal", "xtro", "drer", "cint", "blan", "hery", "htub",
                                  "lvar", "mfra", "pliv", "pmin", "mcap", "nvec", "mlei"))
df_asd <- df[!(is.na(df$category)), ] 
colnames(df_asd) <- c("log2FoldChange", "sp", "Category", "padj", "id", "TPM")

ggplot(df_asd, aes(x = log2FoldChange, y = sp, color = log(TPM), shape = Category))+
  geom_jitter(alpha = 0.6) +
  theme_bw() +
  scale_color_viridis(direction = -1) +
  geom_vline(xintercept = 2, linetype="dotted", 
             color = "black", size=1, alpha = 0.5) +
  geom_vline(xintercept =  -2, linetype="dotted", 
             color = "black", size=1, alpha = 0.5) +
  geom_vline(xintercept = 0, size = 1, alpha = 0.5) +
  ylab("") +
  theme(text = element_text(size = 15))
ggsave(file="~/Desktop/log2FoldChange.svg")

##### SAVE ENVIORMENT #####

#save.image("~/Documents/Gene_expr_evol/Intermediate_files/GO.RData")
#load("~/Documents/Gene_expr_evol/Intermediate_files/GO.RData")
