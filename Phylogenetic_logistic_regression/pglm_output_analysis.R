#Script used to analyse the output from the logistic model fittings
#It will use the fitted models as input, will filter them according to fitting criteria and perform GO analysis
#Apart from that it requires the orthology mapping file and the GO crossreference files
#As output it will output GO enrichment plot and a summarising figure

#Read in libraries
library(tidyverse)
library(RColorBrewer)
library(ggplotify)
library(clusterProfiler)

#Set working directory
setwd("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/pglm/")

#Read in annotation data
#Annotation for orthogroups using a python scripts from:
#https://github.com/davidemms/OrthoFinder/issues/451
annot <- read.table("~/Documents/Gene_expr_evol/Models/N0_blasted.tsv", header = T, sep = "\t", quote = "")
annot <- tibble(annot) %>% 
  filter(Best_hit != "NO_HIT") %>%
  mutate(across(E_val:len, as.double)) %>% 
  mutate(perc = iden/len) %>% 
  group_by(OG) %>% 
  arrange(E_val, desc(perc)) %>% 
  top_n(1) %>% 
  arrange(OG) 

#Get output tables
files <- list.files("./output_fc/", pattern = "pglm", full.names = T)

#Initialize empty data frame
pglm <- data.frame()

#Loop along output tables and read them in
for(i in seq_along(files)){
  
  #Adding information of what was tested extracted from table names
  category <- str_split(files[i], "_")[[1]][4]
  repr <- str_remove(str_split(files[i], "_")[[1]][5], ".tsv")
  
  temp_files <- read_tsv(files[i], col_names = T) %>%
    mutate(Category = category,
           Mode = repr) 

  pglm <- rbind(temp_files, pglm)
    
}

#Filtering for significant fits
pglm <- pglm %>%
  filter(LRT <= 0.05 & 
         wAIC >= 0.5 & 
         binOne > 3 &
         !(is.na(Model)))

pglm$annot <- annot$Best_hit[match(pglm$OG, annot$OG)]

#How many orthogroups are left after filtering?
pglm %>%
  group_by(Category, Mode) %>%
  summarise(N = n())

##### GO ANALYSIS #####

#Read in data from accs.R
t2g <- read.table("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/GO/TERM2GENE.tsv",
                  sep = "\t", header = T)
t2n <- read.table("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/GO/TERM2NAME.tsv",
                  sep = "\t", header = T)

#Initialize empty vector to use for merge_results downstream
#Enrichment results are saved in this list
list_out_expr <- vector(mode = "list", length = length(unique(pglm$Mode)))

#Loop over each model and find enriched terms
for(i in 1:length(unique(pglm$Mode))){
  
  temp_mode <- unique(pglm$Mode)[i]
  temp_df <- pglm %>%
    filter(Mode == temp_mode & Category == "expr") %>%
    mutate(annot = str_remove_all(str_extract(annot, "sp\\|.*\\|"), "sp\\||\\|"))
  
  enrich_res <- enricher(temp_df$annot, 
                         pvalueCutoff = 0.05, 
                         pAdjustMethod = "fdr", 
                         qvalueCutoff = 0.05, 
                         TERM2GENE = t2g, 
                         TERM2NAME = t2n)
  
  enrich_res@result <- subset(enrich_res@result, enrich_res@result$qvalue <= 0.05)
  
  
  list_out_expr[[i]] <- enrich_res
  names(list_out_expr)[i] <- paste0(temp_mode, "_expr")
}


list_out_fc <- vector(mode = "list", length = length(unique(pglm$Mode)))

#Loop over each model and find enriched terms
for(i in 1:length(unique(pglm$Mode))){
  
  temp_mode <- unique(pglm$Mode)[i]
  temp_df <- pglm %>%
    filter(Mode == temp_mode & Category == "fc") %>%
    mutate(annot = str_remove_all(str_extract(annot, "sp\\|.*\\|"), "sp\\||\\|"))
  
  enrich_res <- enricher(temp_df$annot, 
                         pvalueCutoff = 0.05, 
                         pAdjustMethod = "fdr", 
                         qvalueCutoff = 0.05, 
                         TERM2GENE = t2g, 
                         TERM2NAME = t2n)
  
  enrich_res@result <- subset(enrich_res@result, enrich_res@result$qvalue <= 0.05)
  
  
  list_out_fc[[i]] <- enrich_res
  names(list_out_fc)[i] <- paste0(temp_mode, "_fc")
}

list_out <- c(list_out_expr, list_out_fc)

#Merge enriched results for comparisons with a dotplot
enrich_merge <- merge_result(list_out)

#Plot
dotplot(enrich_merge, showCategory = 10) +
  theme(axis.text.y = element_text(size = 12))

##### VISUALIZATIONS #####

p1 <- pglm %>%
  group_by(Mode, Model) %>%
  mutate(Model = factor(Model, levels = unique(Model))) %>%
  mutate(Model = str_remove(Model, "fit"),
         Category = case_when(Category == "fc" ~ "Fold change",
                              Category == "expr" ~ "Expression"),
         Mode = case_when(Mode == "vivipar" ~ "Hemotrophic\r\nviviparous",
                          Mode == "ovipar" ~ "Oviparous",
                          Mode == "ovulipar" ~ "Ovuliparous")) %>%
  ggplot(aes(x = Mode, fill = Model)) +
  geom_bar(color = "black", alpha = 0.7, width = 0.7) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        #legend.position = c(0.2, 0.85),
        legend.background = element_rect(size = 0.5, colour = NA),
        strip.background = element_rect(fill="white")) +
  xlab(NULL) +
  ylab("Count of orthogroups")



p2 <- pglm %>%
  group_by(Mode, Model) %>%
  mutate(Model = factor(Model, levels = unique(Model))) %>%
  mutate(Model = str_remove(Model, "fit"),
         Category = case_when(Category == "fc" ~ "Fold change",
                              Category == "expr" ~ "Expression"),
         Mode = case_when(Mode == "vivipar" ~ "Hemotrophic\r\nviviparous",
                          Mode == "ovipar" ~ "Oviparous",
                          Mode == "ovulipar" ~ "Ovuliparous")) %>%
  mutate(coef_slope_percent = exp(Coeff_slope)) %>%
  ggplot(aes(x = Mode, y = coef_slope_percent, color = Model)) +
  geom_jitter(size = 2, alpha = 0.7, position = position_jitter(width = 0.3, height = 0.4)) +
  facet_wrap(~Category, nrow = 2) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        legend.position = 'none',
        legend.background = element_rect(size = 0.5, colour = "black"),
        strip.background = element_rect(fill="white", color = NA)) +
  xlab(NULL) +
  ylab("Odds ratio") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1)

svg("~/Desktop/Publication_plots/Pglm_model_summary.svg",
    height = 6, width = 10)
p1 + p2 + plot_annotation(tag_levels = 'A') 
dev.off()

