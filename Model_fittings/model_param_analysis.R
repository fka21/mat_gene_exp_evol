#Script used for analysing the output of the model fitting phase
#As input it requires the table outputs from the model fitting, the orthogroup annotations, the phylosignal analysis, and for the enrichment the generated GO crossreference tables 
#As output it will generate the plots found in the manuscript

##### LOAD LIBRARIES #####
library(tidyverse)
library(broom)
library(RColorBrewer)
library(ggplot2)
library(rstatix)
library(patchwork)

#### READ IN DATA #####

#Read in and prepare parameter estimates for each orthogroup
param <- read.table("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/Model_fittings/multiregime_sandbox/output/expr_param_output.tsv",
                    sep = "\t", header = T)
param_fc <- read.table("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/Model_fittings/multiregime_sandbox_fc/output/fc_param_output.tsv",
                       sep = "\t", header = T)
param$Category <- "Expression"
param_fc$Category <- "Fold change"

param_df <- rbind(param, param_fc)
param_df <- param_df[!(is.na(param_df$Model)), ]

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

param_df$annot <- annot$Best_hit[match(param_df$OG, annot$OG)]

#Read phylosignal
phylosig <- read.table("~/Desktop/Publication_plots/Intermediate_data/phylosignal_res.tsv", sep = "\t", header = T)
colnames(phylosig) <- c("OG_phylosig", "Blomberg_K", "Blomberg_pval")

#Combine with phylosig
param_df <- tibble(param_df) %>% 
  full_join(., tibble(phylosig), by = c("OG" = "OG_phylosig"), keep = T)

#Read in and prepare fitting statistics for each orthogroup
res_df_exprs <- read.table("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/Model_fittings/multiregime_sandbox/exprs_cont_fits.tsv", 
                           sep = "\t", header = T) %>%
  mutate(Category = "Expression")
res_df_fc <- read.table("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/Model_fittings/multiregime_sandbox_fc/exprs_fc_fits.tsv", 
                           sep = "\t", header = T) %>%
  mutate(Category = "Fold change")

res_df <- rbind(res_df_exprs, res_df_fc)
res_df$Model <- str_remove(res_df$Model, "^fit")

#What is the distribution of the AICc weights for each model before filtering?
svg("~/Desktop/Publication_plots/aicc_multimodel.svg", height = 6, width = 10)
as_tibble(res_df) %>% 
  filter(!(is.na(OG))) %>% 
  ggplot(aes(x = aic.weight, fill = Category)) +
  geom_histogram(alpha = 0.8, color = "black") +
  facet_wrap(~Model, scales = "free_y", ncol = 2) +
  theme_bw() +
  xlab("AICc weights") +
  ylab("Count of orthogroups") +
  scale_fill_brewer(palette="Set1") +
    theme(text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.position = c(0.08, 0.11),
          legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
          strip.background = element_rect(fill="white"),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.5, 'cm')) +
  geom_vline(aes(xintercept = 0.5), linewidth = 1, linetype = "dashed")
dev.off()    

##### EXPLORATORY PLOTS #####

param_df <- filter(param_df, !(is.na(Model)))

#First: how many of the fitted models actually show phylogeenetic signal according to the permutation test
#Discarding orthogroups where there is was no Blomberg's K performed due to calculation errors
p1 <- param_df %>% 
  filter(Category == "Expression") %>%
  mutate(`Phylogenetic signal` = case_when(Blomberg_pval <= 0.05 ~ "Significant", Blomberg_pval > 0.05 ~ "Not significant")) %>% 
  group_by(Model, `Phylogenetic signal`) %>% 
  dplyr::summarise(Proportion = n()) %>% 
  na.omit() %>%
  ggplot(aes(y = reorder(Model, Proportion), x = Proportion, fill = `Phylogenetic signal`)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  theme_bw() +
  ylab(NULL) +
  xlab("Count of orthogroups") +
  scale_fill_brewer(palette="Set2") +
  theme(legend.position = c(0.75, 0.17),
        text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

#Checking the distribution of K statistic for winning models
library(ggridges)

p2 <- param_df %>% 
  filter(Category == "Expression") %>%
  filter(Blomberg_pval <= 0.05) %>%
  ggplot(aes(x = Blomberg_K, y = Model, fill = Model)) +
  geom_density_ridges(alpha = 0.65, stat = 'binline', scale = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1.2) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  ylab(NULL) +
  xlab("Blomberg's K") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.background = element_rect(size = 0.5, colour = "black")) +
  xlim(c(0, 1.6))

svg("~/Desktop/Publication_plots/Winning_model_phylosignal_composite.svg", height = 6, width = 10)
p1 + p2 + plot_annotation(tag_levels = "A")
dev.off()

#What is the support for the winning models?
svg("~/Desktop/Publication_plots/aicc_multimodel.svg")
param_df %>% 
  ggplot(aes(x = Blomberg_K, y = Model, fill = Model)) +
  geom_density_ridges(alpha = 0.65, stat = 'binline', scale = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1.2) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  ylab(NULL) +
  xlab("Blomberg's K") +
  theme(legend.position = 'none',
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.background = element_rect(size = 0.5, colour = "black")) +
  xlim(c(0, 1.6))
dev.off()


##### GO ANALYSIS #####
library(clusterProfiler)

#Read in data from accs.R
t2g <- read.table("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/GO/TERM2GENE.tsv",
                  sep = "\t", header = T)
t2n <- read.table("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/GO/TERM2NAME.tsv",
                  sep = "\t", header = T)

#Initialize empty vector to use for merge_results downstream
#Enrichment results are saved in this list
list_out <- vector(mode = "list", length = length(unique(param_df$Model)))

#Loop over each model and find enriched terms
  for(i in 1:length(unique(param_df$Model))){
    
    temp_model <- unique(param_df$Model)[i]
    temp_df <- param_df
    
    temp_df$annot <- str_remove_all(str_extract(temp_df$annot, "sp\\|.*\\|"), "sp\\||\\|")
    temp_df <- subset(temp_df, temp_df$Model %in% temp_model)
    
    enrich_res <- enricher(temp_df$annot, 
                           pvalueCutoff = 0.05, 
                           pAdjustMethod = "fdr", 
                           qvalueCutoff = 0.05, 
                           TERM2GENE = t2g, 
                           TERM2NAME = t2n)
    
    enrich_res@result <- subset(enrich_res@result, enrich_res@result$qvalue <= 0.05)
    
    
    list_out[[i]] <- enrich_res
    names(list_out)[i] <- temp_model
  }

#Merge enriched results for comparisons with a dotplot
enrich_merge <- merge_result(list_out)

#Plot
svg("~/Desktop/Publication_plots/GO_signif_models.svg")
dotplot(enrich_merge, showCategory = 8) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  xlab(NULL) +
  labs(size="Gene ratio", color = "Adjusted p-value") 
dev.off()

#Which terms are unique in each model?

enrich_merge@compareClusterResult %>% 
  select(1, 3) %>% 
  as_tibble() %>% 
  group_by(Description) %>% 
  filter(n() == 1) %>% 
  View()

##### INSPECT BASICS #####

#How many shared orthogroups between expression and fold change datasets
library(ggvenn)

venn <- list(Expression = param_df$OG[param_df$Category == "Expression"],
             "Fold change" = param_df$OG[param_df$Category == "Fold change"])

#What models are best fitting for orthogroups where the expression and fold change data have the same 
#best fitting model?
overlap <- venn$Expression[venn$Expression %in% venn$`Fold change`]
param_df %>%
  filter(OG %in% overlap) %>%
  select(OG, Model, Category) %>%
  spread(key = OG, value = Model) %>%
  select(-1) %>%
  t() %>% as.data.frame() %>%
  rename("E" = 1, "F" = 2) %>% as.tibble() %>%
  mutate(asd = case_when(E == F ~ 1,
                         E != F ~ 0)) %>%
  filter(asd == 1) %>%
  group_by(E) %>%
  summarise(N = n()) %>%
  arrange(desc(N))
  
 
p4 <- ggvenn(
  venn, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  set_name_size = 4, text_size = 4,
  auto_scale = T
)


#Inspect how many models for each category
nrs <- param_df %>% 
  as_tibble() %>% 
  filter(!(is.na(OG))) %>%
  group_by(Model, Category) %>% 
  summarize(N = n()) %>% 
  arrange(Category, desc(N))

p1 <- nrs %>%
  filter(Category == "Expression") %>%
  ggplot(aes(y = reorder(Model, N), x = N, fill = Model)) +
  geom_bar(stat = 'identity',  color = "black", alpha = 0.8) +
  facet_wrap(~Category, nrow = 2) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 13),
        legend.position = 'none') +
  ylab("Winning model") +
  xlab("Number of orthogroups")

p2 <- nrs %>%
  filter(Category == "Fold change") %>%
  ggplot(aes(y = reorder(Model, N), x = N, fill = Model)) +
  geom_bar(stat = 'identity', color = "black", alpha = 0.8) +
  facet_wrap(~Category, nrow = 2) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 13),
        legend.position = 'none') +
  xlab("Number of orthogroups") +
  ylab(NULL)

p3 <- param_df %>% 
  filter(!(is.na(OG))) %>%
  ggplot(aes(x = Model, y = species_nr, fill = Model)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        legend.position = 'none') +
  ylab("Number of species") +
  xlab(NULL)

svg("~/Desktop/Publication_plots/winning_signif_models.svg", height = 6, width = 10)
(p1|p2) / (p3 + p4)  + plot_annotation(tag_levels = 'A') 
dev.off()



##### EXAMINE MULTIPLE THETA  ####

p1 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUM") & 
           str_detect(name, "theta") & 
           !(str_detect(name, "se"))) %>%
  na.omit() %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  ggplot(aes(y = value, x = mode, fill = mode)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Category, scale = "free_y", nrow = 2) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        legend.position = 'none',
        strip.background = element_rect(fill="white")) +
  xlab(NULL) +
  ylab(expression(Theta)) 


#Plot differences in optima
param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUM", "OUMV", "OUMA", "OUMVA") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Category) %>%
  mutate(diff1 = abs(diff(value, lag = 1)[1]),
         diff2 = abs(diff(value, lag = 1)[2]),
         diff3 = abs(diff(value, lag = 2)[1])) %>%
  select(-value) %>%
  pivot_longer(diff1:diff3) %>%
  ggplot(aes(y = value, x = mode, fill = mode)) +
  geom_violin(color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Model, nrow = 2) +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("abs(", paste(Delta, ")")))) +
  xlab(NULL) 


#Check which reproductive modes have highest optima across categories
p2 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUM") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.max(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.7) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        strip.background = element_rect(fill="white")) +
  ylab("Count of highest optima") +
  xlab(NULL) 

#Same, but for the minimum optimum
p3 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUM") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.min(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.7) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab("Count of lowest optima") +
  xlab(NULL) 

svg("~/Desktop/Publication_plots/OUM_thetas.svg",
    height = 6, width = 10)
p1 + (p2 / p3) + plot_annotation(tag_levels = 'A') 
dev.off()

#KS test for modes and values

#Are thetas normally distributed?
param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUM", "OUMV", "OUMA", "OUMVA") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(Model, Category) %>%
  shapiro_test(value) # >>> non-normal distributions >>> using non-parametric version of ANOVA: Kruskal-Wallis test

svg("~/Desktop/Publication_plots/theta_cohensD_table.svg",
    height = 7, width = 10)
param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUM", "OUMV", "OUMA", "OUMVA") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(Model, Category) %>%
  cohens_d(value ~ mode) %>%
  select(-1) %>%
  grid.table(., row = NULL)
dev.off()


##### EXAMINE MULTIPLE ALPHA  #####

#Same plots are generated as above with a different parameter
p1 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMA") & str_detect(name, "alpha")) %>%
  na.omit() %>%
  filter(value > 1.000000e-08) %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  ggplot(aes(y = log(value), x = mode, fill = mode)) +
  geom_violin(color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("log(", paste(alpha, ")")))) +
  xlab(NULL)

p2 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMA") & str_detect(name, "alpha")) %>%
  filter(value > 1.000000e-08) %>%
  mutate(value = log(value)) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.max(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("Count of highest ", alpha))) +
  xlab(NULL)

p3 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMA") & str_detect(name, "alpha")) %>%
  filter(value > 1.000000e-08) %>%
  mutate(value = log(value)) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.min(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("Count of lowest ", alpha))) +
  xlab(NULL)

p4 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMA") & 
           str_detect(name, "theta") & 
           !(str_detect(name, "se"))) %>%
  na.omit() %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  ggplot(aes(y = value, x = mode, fill = mode)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Category, scale = "free_y", nrow = 2) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        legend.position = 'none',
        strip.background = element_rect(fill="white")) +
  xlab(NULL) +
  ylab(expression(Theta)) 

#Check which reproductive modes have highest optima across categories
p5 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMA") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.max(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab("Count of highest optima") +
  xlab(NULL) 

#Same, but for the minimum optimum
p6 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMA") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.min(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab("Count of lowest optima") +
  xlab(NULL) 

svg("~/Desktop/Publication_plots/OUMA_alphas.svg",
    height = 11, width = 10)
(p4 + (p5 /p6)) / (p1 + (p2 / p3)) + plot_annotation(tag_levels = 'A') 
dev.off()

svg("~/Desktop/Publication_plots/alpha_cohensD_table.svg",
    height = 8, width = 10)
param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMA") & str_detect(name, "alpha")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(value > 1.000000e-08) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode),
         value = log(value)) %>%
  group_by(Model, Category) %>%
  cohens_d(value ~ mode) %>%
  select(-1) %>%
  grid.table(., row = NULL)
dev.off()


##### EXAMINE MULTIPLE SIGMA2 #####

#Same plots are generated as above with a different parameter
p1 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMV") & str_detect(name, "sigma")) %>%
  na.omit() %>%
  filter(value > 1.000000e-08) %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  ggplot(aes(y = log(value), x = mode, fill = mode)) +
  geom_violin(color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("log(", paste(sigma^2, ")")))) +
  xlab(NULL)

p2 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMV") & str_detect(name, "sigma")) %>%
  filter(value > 1.000000e-08) %>%
  mutate(value = log(value)) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.max(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("Count of highest ", sigma^2))) +
  xlab(NULL)

p3 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMV") & str_detect(name, "sigma")) %>%
  filter(value > 1.000000e-08) %>%
  mutate(value = log(value)) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.min(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("Count of lowest ", sigma^2))) +
  xlab(NULL)

p4 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMV") & 
           str_detect(name, "theta") & 
           !(str_detect(name, "se"))) %>%
  na.omit() %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  ggplot(aes(y = value, x = mode, fill = mode)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Category, scale = "free_y", nrow = 2) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        legend.position = 'none',
        strip.background = element_rect(fill="white")) +
  xlab(NULL) +
  ylab(expression(Theta)) 

#Check which reproductive modes have highest optima across categories
p5 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMV") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.max(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab("Count of highest optima") +
  xlab(NULL) 

#Same, but for the minimum optimum
p6 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMV") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.min(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab("Count of lowest optima") +
  xlab(NULL) 

svg("~/Desktop/Publication_plots/OUMV_sigmas.svg",
    height = 11, width = 10)
(p4 + (p5 /p6)) / (p1 + (p2 / p3)) + plot_annotation(tag_levels = 'A') 
dev.off()

svg("~/Desktop/Publication_plots/sigma_cohensD_table.svg",
    height = 8, width = 10)
param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("OUMV") & str_detect(name, "sigma")) %>%
  filter(str_detect(name, "se", negate = T)) %>%
  filter(value > 1.000000e-08) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode),
         value = log(value)) %>%
  group_by(Model, Category) %>%
  cohens_d(value ~ mode) %>%
  select(-1) %>%
  grid.table(., row = NULL)
dev.off()


##### EXAMINE BMS #####

#Same plots are generated as above with different parameter
p1 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("BMS") & str_detect(name, "sigma")) %>%
  na.omit() %>%
  filter(value > 1.000000e-08) %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  ggplot(aes(y = log(value), x = mode, fill = mode)) +
  geom_violin(color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("log(", paste(sigma^2, ")")))) +
  xlab(NULL)

p2 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("BMS") & str_detect(name, "sigma")) %>%
  filter(value > 1.000000e-08) %>%
  mutate(value = log(value)) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.max(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("Count of highest ", sigma^2))) +
  xlab(NULL)

p3 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("BMS") & str_detect(name, "sigma")) %>%
  filter(value > 1.000000e-08) %>%
  mutate(value = log(value)) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.min(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab(expression(paste("Count of lowest ", sigma^2))) +
  xlab(NULL)

p4 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("BMS") & 
           str_detect(name, "theta") & 
           !(str_detect(name, "se|anc"))) %>%
  na.omit() %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  ggplot(aes(y = log(value), x = mode, fill = mode)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Category, scale = "free_y", nrow = 2) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        legend.position = 'none',
        strip.background = element_rect(fill="white")) +
  xlab(NULL) +
  ylab(expression(log(Theta))) 

#Check which reproductive modes have highest optima across categories
p5 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("BMS") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se|anc", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.max(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab("Count of highest optima") +
  xlab(NULL) 

#Same, but for the minimum optimum
p6 <- param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("BMS") & str_detect(name, "theta")) %>%
  filter(str_detect(name, "se|anc", negate = T)) %>%
  filter(case_when(Category == "Expression" ~ (value > -5 & value < 13),
                   Category == "Fold change" ~ (value < 10 & value > -10))) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", "\r\n")) %>%
  mutate(mode = str_to_title(mode)) %>%
  group_by(OG, Model, Category) %>%
  slice(which.min(value)) %>%
  ungroup() %>%
  group_by(mode, Category) %>% 
  summarise(N = n()) %>%
  ggplot(aes(y = N, x = mode, fill = mode)) +
  geom_bar(stat = 'identity', color = "black", width = 0.6) +
  facet_wrap(~Category, nrow = 2, scale = "free_y") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.background = element_rect(fill="white")) +
  ylab("Count of lowest optima") +
  xlab(NULL) 

svg("~/Desktop/Publication_plots/BMS_sigmas.svg",
    height = 11, width = 10)
(p4 + (p5 /p6)) / (p1 + (p2 / p3)) + plot_annotation(tag_levels = 'A') 
dev.off()

svg("~/Desktop/Publication_plots/BMS_sigma_cohensD_table.svg",
    height = 8, width = 10)
param_df %>% 
  as_tibble() %>%
  pivot_longer(mean_expr:theta.anc.se) %>%
  filter(Model %in% c("BMS") & str_detect(name, "sigma")) %>%
  filter(str_detect(name, "se|anc", negate = T)) %>%
  filter(value > 1.000000e-08) %>%
  na.omit() %>%
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode),
         value = log(value)) %>%
  group_by(Model, Category) %>%
  cohens_d(value ~ mode) %>%
  select(-1) %>%
  grid.table(., row = NULL)
dev.off()


##### EXAMINE CORRELATION AND ASSOCIATIONS BETWEEN PARAMTERES AND EXPRESSION/TREE SIZE #####
library(ggcorrplot)

#Single regime model correlations with expression and tree size
p1 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c( "OU1")) %>%
  select(mean_expr, species_nr, sigma2, alpha, theta) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Sigma" = 3,
         "Alpha" = 4,
         "Theta" = 5) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(legend.position = 'none',
        text = element_text(size = 12))

p2 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c( "BM1")) %>%
  select(mean_expr, species_nr, sigma2) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Sigma" = 3) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12))


svg("~/Desktop/Publication_plots/single_rate_models_corr.svg")
p1+p2 + plot_annotation(tag_levels = 'A') 
dev.off()

#Multi-regime models correlations with expression and tree size
#OUM
  
p2 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUM") & Category == "Fold change") %>%
  select(c(mean_expr, species_nr, ovuliparity.theta, oviparity.theta, hemotrophic_viviparous.theta,
           sigma2, alpha)) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Ovuliparity\r\ntheta" = 3,
         "Oviparity\r\ntheta" = 4,
         "H.v.\r\ntheta" = 5,
         "Sigma" = 6,
         "Alpha" = 7) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12))

p1 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUM") & Category == "Expression") %>%
  select(c(mean_expr, species_nr, ovuliparity.theta, oviparity.theta, hemotrophic_viviparous.theta,
           sigma2, alpha)) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Ovuliparity\r\ntheta" = 3,
         "Oviparity\r\ntheta" = 4,
         "H.v.\r\ntheta" = 5,
         "Sigma" = 6,
         "Alpha" = 7) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12),
        legend.position = 'none',)

p3 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUM")) %>%
  pivot_longer(ovuliparity.theta:theta.anc.se) %>%
  filter(value >= -5 & value <= 13 & !(str_detect(name, "se"))) %>%  #Filtering out irrealistic thetas such as TPM > 1.000.000 or TPM 0.01 (latter was set artificially before log transform)
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode)) %>%
  ggplot(aes(x = value, y = mean_expr, fill = "firebrick")) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "firebrick") +
  facet_wrap(~Category + mode, ncol = 3, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17,),
        axis.title.y = element_text(size = 17, margin=margin(r=-20)),
        axis.title.x = element_text(margin=margin(t=10)),
        strip.background = element_rect(fill="white")) +
  xlab(expression(Theta)) +
  ylab("Mean log(TPM)")

svg("~/Desktop/Publication_plots/corr_param_OUM.svg",
    height = 10, width = 10)
(p1 + p2) / p3 + plot_annotation(tag_levels = 'A') 
dev.off()

param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUM") & Category %in% c("Fold change")) %>%
  pivot_longer(ovuliparity.theta:theta.anc.se) %>%
  filter(value >= -5 & value <= 13 & !(str_detect(name, "se"))) %>%  #Filtering out irrealistic thetas such as TPM > 1.000.000 or TPM 0.01 (latter was set artificially before log transform)
  filter(name == "ovuliparity.theta") %>%
  group_modify(~broom::tidy(lm(mean_expr ~ value, data = .x)))

#BMS

p1 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("BMS") & Category == "Expression") %>%
  select(c(mean_expr, species_nr, ovuliparity.sigma2, oviparity.sigma2, hemotrophic_viviparous.sigma2,
           ovuliparity.theta, oviparity.theta, hemotrophic_viviparous.theta)) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Ovuliparity\r\nsigma2" = 3,
         "Oviparity\r\nsigma2" = 4,
         "H.v.\r\nsigma2" = 5,
         "Ovuliparity\r\ntheta" = 6,
         "Oviparity\r\ntheta" = 7,
         "H.v.\r\ntheta" = 8) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        legend.position = 'none')

p2 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("BMS") & Category == "Fold change") %>%
  select(c(mean_expr, species_nr, ovuliparity.sigma2, oviparity.sigma2, hemotrophic_viviparous.sigma2,
           ovuliparity.theta, oviparity.theta, hemotrophic_viviparous.theta)) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Ovuliparity\r\nsigma2" = 3,
         "Oviparity\r\nsigma2" = 4,
         "H.v.\r\nsigma2" = 5,
         "Ovuliparity\r\ntheta" = 6,
         "Oviparity\r\ntheta" = 7,
         "H.v.\r\ntheta" = 8) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 11))

p3 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("BMS")) %>%
  mutate(across(ovuliparity.sigma2:hemotrophic_viviparous.sigma2, log)) %>%
  pivot_longer(ovuliparity.sigma2:theta.anc.se) %>%
  filter(!(str_detect(name, "se|anc")) & !(is.na(value))) %>% 
  group_by(Category, name) %>%
  filter(!(abs(value - median(value)) > 2*sd(value))) %>% #Removing outliers
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode),
         mode = case_when(mode == "Hemotrophic Viviparous" ~ "H.v.",
                          T ~ mode)) %>%
  ggplot(aes(x = value, y = mean_expr, fill = "firebrick")) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "firebrick") +
  facet_wrap(~Category + mode + param, ncol =  4, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.background = element_rect(fill="white")) +
  xlab(NULL) +
  ylab("Mean log(TPM)")

svg("~/Desktop/Publication_plots/corr_param_BMS.svg",
    height = 13, width = 10)
(p1 + p2) / p3 + plot_annotation(tag_levels = 'A') 
dev.off()

#OUMV
p1 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUMV") & Category == "Expression") %>%
  select(c(mean_expr, species_nr, ovuliparity.sigma2, oviparity.sigma2, hemotrophic_viviparous.sigma2,
           ovuliparity.theta, oviparity.theta, hemotrophic_viviparous.theta)) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Ovuliparity\r\nsigma2" = 3,
         "Oviparity\r\nsigma2" = 4,
         "H.v.\r\nsigma2" = 5,
         "Ovuliparity\r\ntheta" = 6,
         "Oviparity\r\ntheta" = 7,
         "H.v.\r\ntheta" = 8) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        legend.position = 'none')

p2 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUMV") & Category == "Fold change") %>%
  select(c(mean_expr, species_nr, ovuliparity.sigma2, oviparity.sigma2, hemotrophic_viviparous.sigma2,
           ovuliparity.theta, oviparity.theta, hemotrophic_viviparous.theta)) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Ovuliparity\r\nsigma2" = 3,
         "Oviparity\r\nsigma2" = 4,
         "H.v.\r\nsigma2" = 5,
         "Ovuliparity\r\ntheta" = 6,
         "Oviparity\r\ntheta" = 7,
         "H.v.\r\ntheta" = 8) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 11))

p3 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUMV")) %>%
  mutate(across(ovuliparity.sigma2:hemotrophic_viviparous.sigma2, log)) %>%
  pivot_longer(ovuliparity.sigma2:theta.anc.se) %>%
  filter(!(str_detect(name, "se|anc")) & !(is.na(value))) %>% 
  group_by(Category, name) %>%
  filter(!(abs(value - median(value)) > 2*sd(value))) %>% #Removing outliers
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode),
         mode = case_when(mode == "Hemotrophic Viviparous" ~ "H.v.",
                          T ~ mode)) %>%
  filter(case_when(param == "theta" ~ value >= -5 & value <= 13,
                   T ~ value == value)) %>%
  ggplot(aes(x = value, y = mean_expr, fill = "firebrick")) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "firebrick") +
  facet_wrap(~Category + mode + param, ncol =  4, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.background = element_rect(fill="white")) +
  xlab(NULL) +
  ylab("Mean log(TPM)")

svg("~/Desktop/Publication_plots/corr_param_OUMV.svg",
    height = 13, width = 10)
(p1 + p2) / p3 + plot_annotation(tag_levels = 'A') 
dev.off()


#OUMA

p1 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUMA") & Category == "Expression") %>%
  select(c(mean_expr, species_nr, ovuliparity.alpha, oviparity.alpha, hemotrophic_viviparous.alpha,
           ovuliparity.theta, oviparity.theta, hemotrophic_viviparous.theta)) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Ovuliparity\r\nalpha" = 3,
         "Oviparity\r\nalpha" = 4,
         "H.v.\r\nalpha" = 5,
         "Ovuliparity\r\ntheta" = 6,
         "Oviparity\r\ntheta" = 7,
         "H.v.\r\ntheta" = 8) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        legend.position = 'none')

p2 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUMA") & Category == "Fold change") %>%
  select(c(mean_expr, species_nr, ovuliparity.alpha, oviparity.alpha, hemotrophic_viviparous.alpha,
           ovuliparity.theta, oviparity.theta, hemotrophic_viviparous.theta)) %>%
  rename("Mean\r\nexpression" = 1,
         "Species\r\nnumber"= 2,
         "Ovuliparity\r\nalpha" = 3,
         "Oviparity\r\nalpha" = 4,
         "H.v.\r\nalpha" = 5,
         "Ovuliparity\r\ntheta" = 6,
         "Oviparity\r\ntheta" = 7,
         "H.v.\r\ntheta" = 8) %>%
  cor_mat(., method = "spearman") %>%
  ggcorrplot(., hc.order = TRUE, outline.col = "white",
             type = "upper", lab = T, insig = "blank") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 11))

p3 <- param_df %>% 
  as_tibble() %>%
  filter(Model %in% c("OUMA")) %>%
  mutate(across(ovuliparity.alpha:hemotrophic_viviparous.alpha, log)) %>%
  pivot_longer(ovuliparity.alpha:theta.anc.se) %>%
  filter(!(str_detect(name, "se|anc")) & !(is.na(value))) %>% 
  group_by(Category, name) %>%
  filter(!(abs(value - median(value)) > 2*sd(value))) %>% #Removing outliers
  separate(name, "\\.", into = c("mode", "param")) %>%
  mutate(mode = str_replace_all(mode, "_", " ")) %>%
  mutate(mode = str_to_title(mode),
         mode = case_when(mode == "Hemotrophic Viviparous" ~ "H.v.",
                          T ~ mode)) %>%
  filter(case_when(param == "theta" ~ value >= -5 & value <= 13,
                   T ~ value == value)) %>%
  ggplot(aes(x = value, y = mean_expr, fill = "firebrick")) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "firebrick") +
  facet_wrap(~Category + mode + param, ncol =  4, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.background = element_rect(fill="white")) +
  xlab(NULL) +
  ylab("Mean log(TPM)")

svg("~/Desktop/Publication_plots/corr_param_OUMA.svg",
    height = 13, width = 10)
(p1 + p2) / p3 + plot_annotation(tag_levels = 'A') 
dev.off()
