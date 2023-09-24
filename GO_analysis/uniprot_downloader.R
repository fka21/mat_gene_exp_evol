#Script used to download GO annotation for all blasted protein accesion from UniPxyrot database
#Requires the N0_blasted.tsv as input and uses the protein accesion present to query the UniProt database
#Outputs two tables required for clusterprofiler enricher() function

#Load required libraries
library(tidyverse)
library(UniprotR)

setwd("~/Documents/Gene_expr_evol/Gene_exprs_modelling_2023/GO/")

#Read in annotation data as it will be used downstream
annot <- read.table("~/Documents/Gene_expr_evol/Models/N0_blasted.tsv", header = T, sep = "\t", quote = "")
annot <- tibble(annot) %>% 
  filter(Best_hit != "NO_HIT") %>%
  mutate(across(E_val:len, as.double)) %>% 
  mutate(perc = iden/len) %>% 
  group_by(OG) %>% 
  arrange(E_val, desc(perc)) %>% 
  top_n(1) %>% 
  arrange(OG)

#Extract uniprot IDs to be identifiable by UniprotR
accs <- unique(str_remove_all(str_extract(annot$Best_hit, "sp\\|.*\\|"), "sp\\||\\|"))

#Loop over IDs to retrieve GO information
go_tbl <- lapply(accs, function(x){ GetProteinGOInfo(x)})

#Combine outputs
go_tbl <- do.call(rbind, go_tbl)
go_tbl$IDs <- rownames(go_tbl)

#Create a custom TERM2GENE table for enricher from clusterprofiler
t2g <- go_tbl %>% 
  separate(Gene.Ontology.IDs, 
           sep = ";", 
           remove = F,
           into = as.character(c(1:max(na.omit(str_count(go_tbl$Gene.Ontology.IDs, "GO")))))) %>%
  select(-Gene.Ontology.IDs) %>% 
  pivot_longer(1:213) %>% 
  na.omit() %>% 
  select(value, IDs) %>%
  mutate(value = str_remove_all(value, " ")) %>%
  rename("Term" = 1, "Gene" = 2)

#Create a custom TERM2NAME table for enricher from clusterprofiler
t2n <- go_tbl %>% 
  select(Gene.Ontology..GO.) %>%
  separate(Gene.Ontology..GO., 
           sep = ";", 
           remove = T,
           into = as.character(c(1:max(na.omit(str_count(go_tbl$Gene.Ontology.IDs, "GO")))))) %>%
  pivot_longer(1:dim(.)[2]) %>% 
  na.omit() %>%
  mutate(value = str_remove_all(value, "^ ")) %>% 
  select(value) %>%
  separate(value,
           sep = " \\[",
           into = c("Name", "Term")) %>%
  mutate(Term = str_remove_all(Term, "\\]")) %>%
  select(Term, Name) 

#Export tables
write.table(t2g, "./TERM2GENE.tsv", 
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

write.table(t2n, "./TERM2NAME.tsv", 
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

save.image("./accs.RData")
