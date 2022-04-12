################################################################################
##  This program receives a data frame of your patways definitions in GMT format
##  PDS table, DEG table  
##  and retrieves a data frame of Pathway, Gene, PDS, logFC, adjPval
##  Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
##
# Input description: 
#    Df_your_Dereg_Pw 
#    PtwDf_defintions 
#    DfDEG

# to do
#   -h options with an example
#   Describe better the inputs, be consistent using gene symbols for example
#   change: colnames(Bigdf)[1:2] <- c("Pathway", "Target") by colnames(Bigdf)[1:2] <- c("Pathway", "Target-Gene")
# Use this style for pasing the arguments https://github.com/raulmejia/General_statistical_tools/blob/main/Correlations/corr_between_only_numerical_variables/Correlation_matrix_n_Heatmap_out_of_a_numeric_matrix.R
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

##########################################
### Functions defined by the user
##########################################
eliminating_nonDEG_from_PwDefinitionDf <- function(input_PwDF_defitions,  myDfDEG){
  # This program convert in NAs all the non DEGs according to all the row names from a given DEG table (myDfDEG) 
  for( k in 1: dim(input_PwDF_defitions )[1] ){
    my1rowdf <- input_PwDF_defitions[k,] 
    DEGs_in_that_pw_positions <- which(  my1rowdf %in% rownames(myDfDEG)) # identifying the positions with DEGs
    DEGs_in_that_pw_contents_bk <-   myvec[ , DEGs_in_that_pw_positions ] # now the DEGs names
    DEGs_in_that_pw_contents_vec <-   as.character(myvec[ , DEGs_in_that_pw_positions ]) 
    my1rowdf_filledNAs <- my1rowdf
    my1rowdf_filledNAs[k,] <- rep(NA,dim(my1rowdf)[2])
    my1rowdf_filledNAs[k,DEGs_in_that_pw_positions ] <- DEGs_in_that_pw_contents_vec
    input_PwDF_defitions[k,]  <- my1rowdf_filledNAs[k,]
  }
  return(input_PwDF_defitions) 
}

################################################
# Data given by the user
###########################################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-d", "--pwdefinitions", type="character", 
                    help="path to the GMT file with with pathway's definitions")
parser$add_argument("-m", "--matrixpds", type="character", 
                    help="file with your PDS (pathifier) matrix, PDS then zscore agains controls and median afterwards")
parser$add_argument("-t", "--dgetable", type="character", 
                    help="file with your differential expression table")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-o", "--outputfile", type="character", 
                    help="output file where you want to store your correlation plots")

args <- parser$parse_args( )

##########################################################
### Reading the data #####################################
##########################################################
PtwDf_defintions <- read.table( file=args$pwdefinitions , sep= "\t", header=FALSE, quote = "" , row.names = 1 , stringsAsFactors = FALSE)
# Path_to_Ptwydb <-"/home/rmejia/Documents/ABCB1_minimal/Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv"
# PtwDf_defintions <- read.table( file=Path_to_Ptwydb , sep= "\t", header=FALSE, quote = "" , row.names = 1 , stringsAsFactors = FALSE)

Df_your_Dereg_Pw <- read.table( file = args$matrixpds , sep= "\t", header=TRUE, row.names = 1, quote = "" , stringsAsFactors = FALSE)
# Path_to_your_Matrix <- "/home/rmejia/Documents/ABCB1_minimal/Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt"
# Df_your_Dereg_Pw <- read.table( file = Path_to_your_Matrix , sep= "\t", header=TRUE, row.names = 1, quote = "" , stringsAsFactors = FALSE)

DfDEG <- read.table( file=args$dgetable , sep= "\t", header=TRUE, quote = "" , row.names = 1, stringsAsFactors = FALSE)
# Path_your_DEG <-"/home/rmejia/Documents/ABCB1_minimal/Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv"
# DfDEG <- read.table( Path_your_DEG , sep= "\t", header=TRUE, quote = "" , row.names = 1, stringsAsFactors = FALSE)

# Deleting the indicator row (of "Normal") in DEG data frame
NORMALpos <- which(rownames(DfDEG) %in% "NORMAL")
DfDEG <- DfDEG[ -NORMALpos , ]

##########################################################
### Getting the genes inside your Deregulated pathways ###
##########################################################
Coincident_positions <- which( rownames( Df_your_Dereg_Pw ) %in% rownames( PtwDf_defintions ) ) # identifying the Dereg pathways (acording to your dataframe) inside your pw definitions
Pw_definitions_with_deregpw  <- PtwDf_defintions[ Coincident_positions , ] # keeping those pathways

#######################################
### Doing sparse your pw definitions keeping just the pw with at least one DEG 
#######################################
Sparse_PwDef_df <- eliminating_nonDEG_from_PwDefinitionDf( Pw_definitions_with_deregpw  , DfDEG)
     # Sparse_PwDef_df[1,1:70] # Just to take a look

a <- Sparse_PwDef_df # Doing the previos Df compact (erasing the rows with all NAs)
Sparse_PwDef_df_ALLNArows_deleted <- a[rowSums(is.na( a ) ) != ncol( a ),  ]

# Converting the compact Df to a list of characters format
compact_df <- t( Sparse_PwDef_df_ALLNArows_deleted ) 
transposed_Df_yourPw_with_DEGs <- as.data.frame( compact_df )
yourPwDef_as_lists <- as.list(  transposed_Df_yourPw_with_DEGs )

# eliminating the NA entries
test_list <- yourPwDef_as_lists
result_list <- lapply(test_list, na.omit)
# converting the list into a vector
vector_char_DEGs_in_DegPws <- unlist( result_list )

Deg_Pathways_with_at_least_one_DEG <- rep(names(result_list), times= as.numeric(lapply(result_list, length)))

#### Pathways from your definitions with at least a Differential Expressed gene
DF_Pws_DGEs <- cbind( Deg_Pathways_with_at_least_one_DEG,  vector_char_DEGs_in_DegPws)
colnames( DF_Pws_DGEs) <- c("Pathways","Target")
rownames( DF_Pws_DGEs) <- NULL

#######################
### Inserting the PDS
#######################
str(DF_Pws_DGEs)
DF_Pws_DGEs[,1]
DF_Pws_DGEs[,2]

#Positions <- which( rownames(Df_your_Dereg_Pw) %in% DF_Pws_DGEs[ ,"Pathways"] )
#Df_your_Dereg_Pw[Positions,] 

Positions_in_DF <- which(  DF_Pws_DGEs[ ,"Pathways"] %in%   rownames(Df_your_Dereg_Pw) )
DF_Pws_DGEs_PDS <- DF_Pws_DGEs[ Positions_in_DF, ]
DF_Pws_DGEs_PDS <- as.data.frame(DF_Pws_DGEs_PDS)

DF_Pws_DGEs_PDS$PDS <- rep(NA,dim(DF_Pws_DGEs_PDS)[1])

for( k in rownames( Df_your_Dereg_Pw ) ){
  positions <- which( as.character( DF_Pws_DGEs_PDS[,"Pathways"])  %in% k )
  DF_Pws_DGEs_PDS$PDS[positions] <- Df_your_Dereg_Pw[k,2]
}




#######################
### Inserting LogFc & Padj
#######################
DF_Pws_DGEs_PDS_LogFC_Padj <- DF_Pws_DGEs_PDS
DF_Pws_DGEs_PDS_LogFC_Padj$log2FC     <- rep( NA, dim( DF_Pws_DGEs_PDS_LogFC_Padj )[1] )
DF_Pws_DGEs_PDS_LogFC_Padj$DEGpvalue  <- rep( NA, dim( DF_Pws_DGEs_PDS_LogFC_Padj )[1] )

# !!! Error why the gene is not in DfDEG ???
k="PDHA1"
for( k in rownames(DfDEG ) ){
  positions <- which( as.character( DF_Pws_DGEs_PDS_LogFC_Padj[,"Target"])  %in% k )
  DF_Pws_DGEs_PDS_LogFC_Padj$log2FC[positions] <- DfDEG[k,"log2FoldChange"] # The DEG is not here !! (in the DfDEG[700:810,]) why ?? it was supposed to be present by construction!
  DF_Pws_DGEs_PDS_LogFC_Padj$DEGpvalue[positions] <- DfDEG[k,"padj"]
}

Bigdf <- DF_Pws_DGEs_PDS_LogFC_Padj



# DfDEG 

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
