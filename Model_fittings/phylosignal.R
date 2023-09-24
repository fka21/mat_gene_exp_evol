#Script for estimating phylogenetic signal present in expression and fold change data
#As input it requires the dated species tree, the expression and fold change data, and the categorization of orthogroups
#As output it will output a table with phylosignal metrics and plots

##### LIBRARYIES ######
library(motmot)
library(castor)
library(MuMIn)
library(tidyverse)
library(geiger)
library(phytools)

##### FUNCTION ######
source("~/Documents/Gene_expr_evol/Scripts/functions.R")

##### READ TREE #####
#For reproducibility
set.seed(7)

tree <- read.tree('/home/ferenkagan/Documents/Gene_expr_evol/Models/congr_sp_tree.dated.tre')
tree$tip.label[tree$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"

OG_df <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Models/batch_OG_DF_norm-final.tsv")
OG_cat <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Models/OG_categories.tsv")
fc <- t(read.table("~/Documents/Gene_expr_evol/fc.tsv", sep = "\t"))

##### MODEL SELECTION #####
variables <- ls(pattern = "*_OG_df")

#Collapsing replicates into single values and calculating standard errors for each species measurement
OG_df_se <- t(SEcollapseReplicate(OG_df))
OG_df <- t(collapseReplicate(OG_df))

phylosignal_res <- data.frame(OG = as.character(),
                              K = as.double(),
                              pval = as.double())

##### TPM PHYLOSIGNAL #####

for(k in 1:dim(OG_df)[2]){

  #Extract orthogroup with its standard errors
  temp <- OG_df[, k, drop = F]
  temp_se <- OG_df_se[, k, drop = F]

  index <- OG_cat %>% tibble %>% filter(ID == colnames(temp) & Category == "M") %>% unique()
    
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
    
  sortedData <- sortTraitData(phy = tree, y = temp, log.trait = F, pass.ultrametric = T)
    
  phylosignal_res[k, 1] <- colnames(temp)
  phylosignal <- phylosig(tree = sortedData$phy,
                          x = sortedData$trait, 
                          se = temp_se_int,
                          test = T,
                          nsim = 1000)
    
  phylosignal_res[k, 2] <- phylosignal$K
  phylosignal_res[k, 3] <- phylosignal$P
  
  phylosignal_res <- na.omit(phylosignal_res)
    
}



##### FC PHYLOSIGNAL #####

phylosignal_res_fc <- data.frame(OG = as.character(),
                              K = as.double(),
                              pval = as.double())

for(k in 1:dim(fc)[2]){
  
  #Extract orthogroup with its standard errors
  temp <- fc[, k, drop = F]
  
  index <- OG_cat %>% tibble %>% filter(ID == colnames(temp) & Category == "M") %>% unique()
  
  #Subsetting for OGs found in maternal gene category
  temp <- as.matrix(subset(temp, rownames(temp) %in% index$species))
  temp <- na.omit(temp)

  #Pointless to use small trees as models will get misspecified
  #Using a tree size of 15 as cutoff
  if(length(temp[, 1]) < 15){next}
  
  sortedData <- sortTraitData(phy = tree, y = temp, log.trait = F, pass.ultrametric = T)
  
  phylosignal_res_fc[k, 1] <- colnames(temp)
  phylosignal <- phylosig(tree = sortedData$phy,
                          x = sortedData$trait, 
                          test = T,
                          nsim = 1000)
  
  phylosignal_res_fc[k, 2] <- phylosignal$K
  phylosignal_res_fc[k, 3] <- phylosignal$P
  
}

save.image("~/Documents/Gene_expr_evol/Intermediate_files/phylosignal.RData")
write.table(phylosignal_res, "~/Desktop/Publication_plots/Intermediate_data/phylosignal_res.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(phylosignal_res_fc, "~/Desktop/Publication_plots/Intermediate_data/phylosignal_res_fc.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


#Plot significant phylogenetic signal

pie_df <- data.frame(
  signif = c("Phylogenetic signal present", "Phylogenetic signal absent"),
  nr = c(dim(tibble(phylosignal_res) %>% filter(pval <= 0.05))[1], dim(tibble(phylosignal_res) %>% filter(pval > 0.05))[1])
)

pie_df <- pie_df %>% 
  arrange(desc(signif)) %>%
  mutate(prop = nr / sum(pie_df$nr) ) %>%
  mutate(ypos = c(2000, 5500) )

svg("~/Desktop/Publication_plots/Blombergs_signif_prop.svg")
ggplot(pie_df, aes(x = "", y = nr, fill = signif)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title= element_blank()) +
  geom_text(aes(x = c(1.1, 1.1), y = ypos, label = scales::percent(prop)), color = "white", size=10) +
  scale_fill_brewer(palette="Set2")
dev.off()

svg("~/Desktop/Publication_plots/Blombergs_K_statistic.svg")
tibble(phylosignal_res) %>% filter(pval < 0.05) %>% 
  ggplot(aes(x = K)) + 
  geom_histogram(fill = "#FC8D62") +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1.5) +
  theme_bw() +
  xlab("K statistic") +
  ylab("Orthogroup count") +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
dev.off()



pie_df <- data.frame(
  signif = c("Phylogenetic signal present", "Phylogenetic signal absent"),
  nr = c(dim(tibble(phylosignal_res_fc) %>% filter(pval <= 0.05))[1], dim(tibble(phylosignal_res_fc) %>% filter(pval > 0.05))[1])
)

pie_df <- pie_df %>% 
  arrange(desc(signif)) %>%
  mutate(prop = nr / sum(pie_df$nr) ) %>%
  mutate(ypos = c(2000, 5500) )

svg("~/Desktop/Publication_plots/Blombergs_signif_prop_fc.svg")
ggplot(pie_df, aes(x = "", y = nr, fill = signif)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title= element_blank()) +
  geom_text(aes(x = c(1.1, 1.1), y = ypos, label = scales::percent(prop)), color = "white", size=10) +
  scale_fill_brewer(palette="Set2")
dev.off()

svg("~/Desktop/Publication_plots/Blombergs_K_statistic_fc.svg")
tibble(phylosignal_res_fc) %>% filter(pval < 0.05) %>% 
  ggplot(aes(x = K)) + 
  geom_histogram(fill = "#FC8D62") +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1.5) +
  theme_bw() +
  xlab("K statistic") +
  ylab("Orthogroup count") +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
dev.off()


