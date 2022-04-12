################################################################################
##  This program receives a data frame with gene symbols column and 
##   annotate their pharmacological targets according to rDGIdb
##  Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
##
#
# Note: pwdefinitions your PtwDf_defintions and your PtwDf_defintions tables should be in the same Pathway database (KEGG , GO , etc )
#
# Input description: 
#    -d dataframe (with a column with genesybols, no rownames ):
#         Pathways        Target  MedianPDSzs     log2FC  DEGpvalue
#         hsa00230 Purine metabolism      AMPD1   1.48432605718692        -2.02502990184976       4.54904680923331e-07
#         hsa00230 Purine metabolism      AK5     1.48432605718692        -1.10711372264512       0.00158257005550845
#         hsa00230 Purine metabolism      ADCY5   1.48432605718692        -1.78027935014617       8.21124977167617e-11
#
#    -n columnname ( name of your genesymbol column) 
#
#   -l  label (the character string without spaces that you want to label the name of your results)
#   
#   -o  outputfolder ( the path where you want to store your results)
#
# Example: 
#
#  Rscript /path/to/Drugtarget_annotation_over_a_Dataframe_with_genesymbols.R \
#  -d /path/to/your/df/with_gs_column.tsv \
#  -n Name_of_your_column_with_genesymbols \
#  -l label_for_your_results \
#  -o /path/to/save/your/results 
#
# to do: improve -h option

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
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("DescTools")) {
  install.packages("DescTools", ask =FALSE)
  library("DescTools")
}

##########################################################
## Functions that will be used 
##########################################################
AnottateDrugsForthisCharacter <- function(mygenes){
  myquery <- queryDGIdb(mygenes)
  return(detailedResults(myquery))
}

###########################################
# Data given by the user          #########
###########################################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-d", "--dataframe", type="character", 
                    help="Your dataframe (with a column with genesybols)")
parser$add_argument("-n", "--columnname", type="character", 
                    help="-n Name of your genesymbol column")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="output file where you want to store your correlation plots")

args <- parser$parse_args( )

##########################################################
### Reading the data #####################################
##########################################################
yourdf <- read.table( file=args$dataframe , sep= "\t", header=TRUE, quote = "" , row.names = NULL , stringsAsFactors = FALSE )
# Path_to_df <-"/DF_Pw_Target_PDS_LogFC_adjpvalTCGA_Basal_under_percentile_25_top20_DEG_log2only.tsv" 
# yourdf <- read.table( file=Path_to_df , sep= "\t", header=TRUE, quote = "" , row.names = NULL , stringsAsFactors = FALSE )

columnwithgenesymbols <-args$columnname
# columnwithgenesymbols <- "Target"

Path_of_Results <- normalizePath( args$outputfolder )
# Path_of_Results <- "Results/BigDfPTD/TCGA"
dir.create( Path_of_Results, recursive = TRUE)

Label_for_Results <- args$label
# Label_for_Results  <- "Drug-target_annotated"

################################
### Annotating the drugs
################################
DF_annotated_gene_symbols <- AnottateDrugsForthisCharacter(  yourdf[,columnwithgenesymbols ] )

Comment_1 <- paste( length(unique(DF_annotated_gene_symbols$Gene)), "gene symbols had an annotation in rDGIdb", "out of the",
       length(unique(yourdf[,columnwithgenesymbols ])),"gene symbols in your data frame")
print(Comment_1)

# Deleting the orphan GeneSymbols (the ones without an Drug-target annotation)
Positions_annotated <- which( yourdf[,columnwithgenesymbols ]  %in% DF_annotated_gene_symbols$Gene )
df_not_orphans <- yourdf[ Positions_annotated,  ]
rownames(df_not_orphans) = NULL

##############################
### Looking for the clave
##############################
La_clave <- table( DF_annotated_gene_symbols$Gene )

df_not_orphans$clave <- rep( NA, dim(df_not_orphans)[1])

for( k in df_not_orphans[ , "Target"]  ){
  positions <- which( df_not_orphans[ , "Target"] %in%  k )
  df_not_orphans$clave[positions] <- La_clave[k]
}

######################################################
### Replicating the vectors according to "La clave"
######################################################
replicating <- function(myvector ){
  rep(myvector,as.numeric(myvector["clave"]))
}

b <- apply( df_not_orphans , 1, replicating ) 

vector2matrix <- function( myvec ){ 
  Reshaped_to_a_matrix <- t(matrix( myvec , nrow= ncol(df_not_orphans) )  )
  colnames( Reshaped_to_a_matrix ) <- names( myvec ) [1:ncol(df_not_orphans)]
  return( Reshaped_to_a_matrix)
  }

list_of_matrices <- lapply( b , vector2matrix)

df_not_orphans_space_4_annotations <- do.call( rbind, list_of_matrices )

######################################################
### Filling the bigdf now with space for annotations 
######################################################
Bigdf <- df_not_orphans_space_4_annotations
Bigdf <- as.data.frame( Bigdf)

Bigdf$Drug <- rep( NA, dim( Bigdf )[1] )
Bigdf$InteractionType <- rep( NA, dim( Bigdf )[1] )
Bigdf$Source <- rep( NA, dim( Bigdf )[1] )
Bigdf$PMIDs <- rep( NA, dim( Bigdf )[1] )

for( k in unique( DF_annotated_gene_symbols$Gene )  ){
  positionsBigDf <-       which(  as.character( Bigdf[ , columnwithgenesymbols ] )  %in% k )
  positionsDF_annot_gs <- which(  as.character( DF_annotated_gene_symbols$Gene  )   %in% k )
  Bigdf[ positionsBigDf, c("Drug","InteractionType","Source","PMIDs") ] <- DF_annotated_gene_symbols[ positionsDF_annot_gs , c("Drug","InteractionType","Source","PMIDs") ]
}

if( any(grepl( "clave", colnames(Bigdf)) ) == TRUE ) {
  Bigdf <- Bigdf[ , -grep( "clave", colnames(Bigdf) ) ]  
}

##########################################################
###### Saving the results
##########################################################
write.table(Bigdf ,
            file=paste0(Path_of_Results,"/","DF_Pw_Target_PDS_LogFC_a.pval_Drug_InteractionType_Source_PMIDs",Label_for_Results,".tsv" ),
            sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE  )
write.table( DF_annotated_gene_symbols ,
            file=paste0(Path_of_Results,"/","DF_annotated_rDGIdbgene_symbols",Label_for_Results,".tsv" ),
            sep="\t", quote=FALSE , row.names= FALSE, col.names= TRUE  )
