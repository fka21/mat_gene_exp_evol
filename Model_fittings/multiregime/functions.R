#Function used in different scripts

pgls.SEy<-function(model,data,corClass=corBrownian,tree=tree,
                   se=NULL,method=c("REML","ML"),interval=c(0,1000),...){
  corfunc<-corClass
  ## preliminaries
  data<-data[tree$tip.label,]
  if(is.null(se)) se<-setNames(rep(0,Ntip(tree),
                                   tree$tip.label))
  ## likelihood function
  lk<-function(sig2e,data,tree,model,ve){
    tree$edge.length<-tree$edge.length*sig2e
    ii<-sapply(1:Ntip(tree),function(x,e) which(e==x),
               e=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+
      ve[tree$tip.label]
    v<-diag(vcv(tree))
    vf<-varFixed(~v)
    COR<-corfunc(0.1,tree,...)
    fit<-gls(model,data=data,correlation=COR,weights=vf)
    -logLik(fit)
  }
  ## estimate sig2[e]
  fit<-optimize(lk,interval=interval,
                data=data,tree=tree,model=model,ve=se^2)
  tree$edge.length<-tree$edge.length*fit$minimum
  ii<-sapply(1:Ntip(tree),function(x,e) which(e==x),
             e=tree$edge[,2])
  tree$edge.length[ii]<-tree$edge.length[ii]+
    se[tree$tip.label]^2
  v<-diag(vcv(tree))
  vf<-varFixed(~v)
  ## fit & return model
  gls(model,data,correlation=corfunc(0.1,tree),weights=vf,
      method=method)
}

#Calculate std error
std_mean <- function(x) sd(x)/sqrt(length(x))

#A function which collapses technical replicates into their means
#Given a data.frame where replicates are indicated in the names by "_[0-9]+" it calculates the means across replaicates
#and returns a single mean value for that grouping
collapseReplicate <- function(df){
  temp <- as.data.frame(t(df))
  
  temp$ID <- str_remove_all(rownames(temp), "_[0-9]+")
  
  temp2 <- tibble(temp) %>% dplyr::group_by(ID) %>%
    dplyr::summarize(across(1:dim(.)[2]-1, \(x) mean(x, na.rm = TRUE)))
  
  temp <- as.data.frame(temp2)[, -1]
  rownames(temp) <- temp2$ID; temp <- t(temp)
  return(temp)
}

#A function which collapses technical replicates into their standard errors
#Given a data.frame where replicates are indicated in the names by "_[0-9]+" it calculates the standard errors across
#replaicates and returns a single standard error value for that grouping
SEcollapseReplicate <- function(df){
  temp <- as.data.frame(t(df))
  
  temp$ID <- rownames(temp)
  temp$ID <- str_remove_all(temp$ID, "_[0-9]+")
  
  temp_2 <- tibble(temp) %>% group_by(ID) %>% dplyr::summarize(across(1:dim(.)[2]-1, std_mean))
  temp <- as.data.frame(temp_2)[,-1, drop = F]
  rownames(temp) <- temp_2$ID; temp <- t(temp)
  return(temp)
}



#Function to estimate MZT timing
#Inputs are variance stabilized objects and a integer
#Adjust k to cut sample distance dendogram clusters (k >= 2)

mzt.dist <- function(vst_obj, k){
  
  #Calculate distance of samples
  sampleDists <- as.matrix(dist(t(assay(vst_obj)), method = 'euclidean'))
  
  #Setup rownames and plot
  rownames(sampleDists) <- as.character(rownames(colData(vst_obj)))
  colnames(sampleDists) <- NULL
  (p <- pheatmap(sampleDists, col=colors, fontsize = 20))
  
  #Cut dendogram at given clusters
  temp <- sort(cutree(p$tree_row, k))
  
  #To return to samples beloning to the cluster with stage1
  temp_grp <- temp[grepl(names(temp), pattern = "oocyte")]
  temp <- temp[temp == temp_grp[1]]
  temp <- unique(str_remove_all(names(temp), "_[0-9]+$"))
  print(temp)
}

#Loop through all tximport objects, deteremine expressed genes with cut-off approach and return the expressed gene names
#As input it requires a vector or tximport object names (not the object itself)
#It removes tximport object from enviorment (memory heavy due to technical replicates)

expr.extr <- function(variables){
  mat_ids <- c() #Initialize empty result vector
  
  for(i in 1:length(variables)){
    temp <- get(variables[i])$abundance
    colnames(temp) <- str_remove_all(colnames(temp), "_[0-9]+$")
    
    #Using cutoff of TPM >= 2
    temp_expr <- rownames(temp[rowMeans(as.matrix(temp)[, grepl("stage1", colnames(temp)), drop = F]) >= 2, ])
    
    #Subsetting TPM table for expressed genes
    temp <- subset(temp, rownames(temp) %in% temp_expr)
    #To match OrthoFinder output need to format gene IDs in some datasets
    if(!(variables == "txi_sfel")){
      temp <- str_remove_all(rownames(temp), "\\.[0-9]+$")
    } else {temp <- rownames(temp)}
    
    mat_ids <- c(mat_ids, temp)
  }
  
  return(mat_ids)
}

#Loop through all tximport objects and save RDS of the abundance table
#As input it requires a vector or tximport object names (not the object itself)

export.txtable <- function(variables){
  for(i in 1:length(variables)){
    temp <- get(variables[i])
    
    saveRDS(temp, paste("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files", paste0(variables[i], ".RDS"), sep = "/"))
    assign(variables[i], temp$abundance)
    
  }
}
