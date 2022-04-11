################################################################################
##  This program receives a data frame of your preference with pathways as rownames
##  a tabular definition of such pathways and a Data frame with DEG 
##  and retrieves a data frame of Pathway, Gene, Drugs
##  Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
##
# Input description: 
#    Df_your_Dereg_Pw 
#    PtwDf_defintions 
#    DfDEG

# to do change: colnames(Bigdf)[1:2] <- c("Pathway", "Target") by colnames(Bigdf)[1:2] <- c("Pathway", "Target-Gene")
# keep all the columns from: Bigdf <- ldply(list_GnDg, data.frame)

################################################################################
##   Installing and loading the required libraries              ################
################################################################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("utils")) {
  install.packages("utils", ask =FALSE)
  library("utils")
}
if (!require("rDGIdb")) {
  BiocManager::install("rDGIdb", ask =FALSE)
  library(rDGIdb)
}
if (!require("plyr")) {
  BiocManager::install("plyr", ask =FALSE)
  library(plyr)
}
################################################
# Data given by the user
###########################################
args <- commandArgs(trailingOnly = TRUE)
Path_to_Ptwydb <- args[1] # The path to the information about the samples
# Path_to_Ptwydb <-c("../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv"); Path_to_your_Matrix <-c("../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix.txt_top_50.tsv"); Path_your_DEG <- c("../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv"); Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/BigDfPTD/TCGA/") ; Label_for_Results <-"TCGA_Basal_under_percentile_25_ALLpw_DEG_log2only"
# Path_to_Ptwydb <-"/home/rmejia/Documents/ABCB1_minimal/Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv"

Path_to_your_Matrix <- args[2] # The path to your matrix
# Path_to_your_Matrix <-c("TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt");
# Path_to_your_Matrix <- "/home/rmejia/Documents/ABCB1_minimal/Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt"

Path_your_DEG <- args[3]
# Path_your_DEG <- c("../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv");
# Path_your_DEG <-"/home/rmejia/Documents/ABCB1_minimal/Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv"

Path_of_Code <- args[4] # The path to your code
# Path_of_Code<-c("./")
# Path_of_Code <- "/home/rmejia/Documents/ABCB1_minimal/ABCB1_minimal_code"

Path_of_Results <-args[5] # # where do you want to save your results?
# Path_of_Results<-c("../Results/BigDfPTD/TCGA/")
# Path_of_Results <- "/home/rmejia/Documents/ABCB1_minimal/Results/BigDfPTD/TCGA/"

Label_for_Results<-args[6] # Label for your results
# Label_for_Results <- "Df_pw_gn_TCGA"
# Label_for_Results <- "TCGA_Basal_under_percentile_25_top20_DEG_log2only"

## NOTE: for windows users #  If you use windows and you have troubles with your paths, please try out this tool: # Path_of_x <-choose.dir()

##########################################################
## Functions that will be used 
##########################################################
AnottateDrugsForthisCharacter <- function(mygenes){
  myquery <- queryDGIdb(mygenes)
  return(detailedResults(myquery))
}

##########################################################
### Reading the data #####################################
##########################################################
Df_your_Dereg_Pw <- read.table( Path_to_your_Matrix, sep= "\t", header=TRUE, row.names = 1, quote = "" , stringsAsFactors = FALSE)
PtwDf_defintions <- read.table( Path_to_Ptwydb , sep= "\t", header=TRUE, quote = "" , row.names = 1 , stringsAsFactors = FALSE)
DfDEG <- read.table( Path_your_DEG , sep= "\t", header=TRUE, quote = "" , row.names = 1, stringsAsFactors = FALSE)

# Deleting the indicator row (of "Normal") in DEG data frame
NORMALpos <- which(rownames(DfDEG) %in% "NORMAL")
DfDEG <- DfDEG[ -NORMALpos,]

##########################################################
### Getting the genes inside your Deregulated pathways ####
##########################################################
Coincident_positions <- which( rownames( Df_your_Dereg_Pw ) %in% rownames(PtwDf_defintions) )
Df_yourPw_with_genes <- PtwDf_defintions[Coincident_positions, ]

##########################################################
###### Drugs 
##########################################################
list_GnDg <-apply(Df_yourPw_with_genes ,1 ,AnottateDrugsForthisCharacter)
Bigdf <- ldply(list_GnDg, data.frame) # Making a big data.frame
Bigdf <- Bigdf[,c(1,2,4,5)]
colnames(Bigdf)[1:2] <- c("Pathway", "Target")

##########################################################
###### Adding log2FC and DGE padj column 
##########################################################
Bigdf$log2FC <- rep( NA, dim( Bigdf )[1] )
Bigdf$DEGpvalue <- rep( NA, dim( Bigdf )[1] )
for( k in rownames(DfDEG ) ){
  positions <- which( as.character(Bigdf[,"Target"])  %in% k )
  Bigdf$log2FC[positions] <- DfDEG[k,"log2FoldChange"]
  Bigdf$DEGpvalue[positions] <- DfDEG[k,"padj"]
}

##########################################################
###### Saving the results
##########################################################
dir.create( Path_of_Results, recursive = TRUE)
write.table(Bigdf ,
            file=paste0(Path_of_Results,"DF_PwTgtDrugInteracLFCPval_includingNoDEG",Label_for_Results,".tsv" ),
            sep="\t", quote=FALSE , row.names= FALSE, col.names= TRUE  )

BigdfonlyDEG <- Bigdf[!is.na(Bigdf[,"DEGpvalue"]) ,]
write.table(BigdfonlyDEG ,
            file=paste0(Path_of_Results,"DF_PwTgtDrugInteracLFCPval_OnlyDEG",Label_for_Results,".tsv" ),
            sep="\t", quote=FALSE , row.names= FALSE, col.names= TRUE  )
