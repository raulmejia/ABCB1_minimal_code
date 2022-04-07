# This program receives a folder's path and retrieves KEGG data base
# in its more updated version in a tabular way and a R-object in both KEGG ids and
# genesymbols
#########################################
##### Installing required libraries #####
#########################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library(BiocManager)
}
if (!require("gage")) {
  BiocManager::install("gage", ask =FALSE)
  library(gage)
}
if (!require("KEGGREST")) {
  BiocManager::install("KEGGREST", ask =FALSE)
  library("KEGGREST")
}
if (!require("stringr")) {
  install.packages("stringr", ask =FALSE)
  library("stringr")
}
if (!require("utils")) {
  install.packages("utils", ask =FALSE)
  library("utils")
}
##########################################
#### Defining user  variables #######
##########################################

args <- commandArgs(trailingOnly = TRUE)
results_path <- args[1]
# results_path<-c("~/Documents/4Andre") # You can set a different path to save your results
# choose.dir(default = "", caption = "Select folder")
results_path<-normalizePath(results_path)
################################################
####### Defining functions that will be used ###
################################################

list_of_chr_to_df<-function(X){
  # function that converts a list of character vectors to a data frame filled up with NA in blank spaces
  # X is a list of 310 lists, each "sublist" is a a list of 1 element, that element is a character vector
  mymax<-max(unlist(lapply(X,length)))
  Num_pathways<-length(lapply(X,length))
  mylist<-list()
  for(k in 1:Num_pathways){
    N_faltantes<-mymax-length(X[[k]])
    mylist[[k]]<-c(X[[k]],rep("NA",N_faltantes))
  }
  names(mylist)<-names(X)
  dd<-data.frame(matrix(nrow=Num_pathways,ncol=mymax))
  for(i in 1:length(mylist)){
    dd[i,] <- mylist[[i]]
  }
  rownames(dd)<-names(mylist)
  return(dd)
}

#######################################
###### Data manipulation ##############
#######################################
dir.create(results_path, recursive = TRUE)

kegg_gsets <-kegg.gsets(species = "hsa", id.type = "kegg") # Downloading the most recent keggdb
# WS #  str(kegg_gsets)
# WS #names(kegg_gsets)
# WS #  kegg_gsets$kg.sets
# WS # kegg_gsets$kg.sets[kegg_gsets$sig.idx]

save( kegg_gsets,file=paste0(results_path,"/","uptdateKeggPathways.RData" ))

KEGG_pathways_in_df <-list_of_chr_to_df(kegg_gsets$kg.sets) # shaping the pathway's data into a data.frame (The genes are in KEGG ids)

Kegghsa <- KEGGREST::keggList("hsa") # retrieving  KEGG ids and genesymbols of the human genes
# WS # Kegghsa[1:5]

GeneSymbol_kid <- str_extract(Kegghsa, "[:alnum:]+") # Extracting only the genesymbols
names(GeneSymbol_kid) <- gsub("hsa:","",names(Kegghsa))

kegg_sets_kid_gs <- kegg_gsets$kg.sets
for( w in 1:length(kegg_gsets$kg.sets) ){
  kegg_sets_kid_gs[[ names(kegg_gsets$kg.sets)[w] ]] <- GeneSymbol_kid[which(  names(GeneSymbol_kid) %in%  kegg_gsets$kg.sets[[ names(kegg_gsets$kg.sets)[w]   ]])]
} # Changing the KEGG pathways from KEGG ids to genesymbols
# WS # kegg_sets_kid_gs

kegg_sets_kid_gs <- lapply(kegg_sets_kid_gs,unique) # deleting duplicated gene symbols from the KEGG pathway's definition

KEGG_pathways_in_df_genesymbols <-list_of_chr_to_df( kegg_sets_kid_gs) # shaping the pathway's data into a data.frame (The genes are in genesymbols)
#KEGG_pathways_in_df_genesymbols[1:5,1:5]
#KEGG_pathways_in_df[1:5,1:5]

write.table(KEGG_pathways_in_df_genesymbols, file=paste0(results_path,"/","KEGG_pathways_in_df_genesymbol.tsv") ,sep="\t",col.names = FALSE,quote = FALSE)
write.table(KEGG_pathways_in_df, file=paste0(results_path,"/","KEGG_pathways_in_df_in_keggids.tsv") ,sep="\t",col.names = FALSE,quote = FALSE)

save(KEGG_pathways_in_df_genesymbols,file=paste0(results_path,"/","KEGG_pathways_in_df_in_genesymbols.RData"))
save(KEGG_pathways_in_df,file=paste0(results_path,"/","KEGG_pathways_in_df.RData"))

# WS # save.image(file="myworkspace.RData")
