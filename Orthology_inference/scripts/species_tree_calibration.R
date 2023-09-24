library(geiger)
library(phytools)
library(ape)
library(tidyverse)
library(taxize)
library(ggpubr)
library(rotl)

#Using timetree.org to acquire time calibrated tree
tree <- read.tree("~/Downloads/TOL_for_timetree.nwk")

plot(ladderize(tree))
axisPhylo()

#Need to rename a species as gbresolve() recognizes it under a different name
tree$tip.label[tree$tip.label == "Mesocentrotus_franciscanus"] <- "Strongylocentrotus_franciscanus"
tree$tip.label[tree$tip.label == "Nereis_denhamensis"] <- "Platynereis_dumerilii"
tree$tip.label[tree$tip.label == "Thulinius_stephaniae"] <- "Hypsibius_dujardini"

#Fetch ultrameric tree generated with IQTree2 under fixed topology of ROTL generated tree
species_tree <-  read.tree("~/Documents/Gene_exp_evol/Gene_expr_evolution/Orthofinder/SpTree_from_data_fixed_topology.partitions.treefile")
plot(ladderize(species_tree))
axisPhylo()

#Need to reroot this tree as IQTree2 could not handle outgroups
species_tree <- root(species_tree, "Mnemiopsis_leidyi", resolve.root = T)

#Need to rename a species as gbresolve() recognizes it under a different name
species_tree$tip.label[species_tree$tip.label == "Mesocentrotus_franciscanus"] <- "Strongylocentrotus_franciscanus"

#Following the tutorial with own data
trl <- gbresolve(tree, within="Metazoa", rank=c("genus", "family"))$tax 
ex <- subset(tree, tax=trl, rank="genus")
print(cbind(original=tree$tip.label[1:6], exemplar=ex$tip.label[1:6]))

plot.phylo(ladderize(ex, right=FALSE), type="phylogram", cex=1, label.offset=3) 
axisPhylo(cex.axis=0.75) 

pw <- ladderize(species_tree, right=FALSE) 
tmp <- gbresolve(pw, within="Metazoa", rank=c("genus", "class")) 
tax <- tmp$tax 
head(tax[,c("family", "order")]) 

plot.phylo(pw, type="fan", show.tip=FALSE, edge.width=0.2, edge.color="gray", no.margin=TRUE)

res1 <- congruify.phylo(reference=ex, target=pw, taxonomy=tax, scale=NA)

#Need to set PATH to PATHd8
curr_env <- Sys.getenv("PATH")
new_env <- paste(curr_env, "/home/ferenkagan/PATHd8", sep = ":")
Sys.setenv(PATH = new_env)

res <- congruify.phylo(reference=ex, target=pw, taxonomy=tax, scale="PATHd8")

cal <- res$calibrations 

#Get # of sites from IQTree2 log file
write.treePL(pw, cal, base="congr_sp_tree", nsites=58733, opts=list(smooth=0.1, nthreads=2, opt=1, optad=1, thorough=TRUE)) 

#After this just run:
#treePL congr_sp_tree.infile

#Comparing to NCBI based fixed topology (has polytomies)
tree1 <- read.tree("~/Metazoa.dated.tre")
tree2 <- read.tree("~/congr_sp_tree.dated.tre")

par(mfrow = c(1,2))
plot(ladderize(tree1), cex = 0.8)
axisPhylo(cex = 0.8)
plot(ladderize(tree2), cex = 0.8)
axisPhylo(cex = 0.8)
par(mfrow = c(1,1))

#Need to rename to match species names of expression table
tree2$tip.label[tree2$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"
write.tree(tree2, "~/congr_sp_tree.dated.tre")
