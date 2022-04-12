######################################################################################
# This program receives a GMT file with the pathways (pws) definitions, a tsv with deregulated Pws and their scores, a DEG table
# the programs does the intersection and retrieves a data frame with the Pw, target(gene), PDS (or what ever single dereg value), LogFC and pval (from the DEG) 
##  Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
######################################################################################
#
# Note: pwdefinitions your PtwDf_defintions and your PtwDf_defintions tables should be in the same Pathway database (KEGG , GO , etc )
#
# Input description: 
#    -d pwdefinitions (a tsv file with a GMT format ):
#        hsa00010 Glycolysis / Gluconeogenesis   LDHAL6A ALDH3B2 PGAM2 ...
#        hsa00020 Citrate cycle (TCA cycle)      PDHA2   IDH3G   MDH1 ...
#        sa00030 Pentose phosphate pathway      TKTL2   IDNK    PRPS1L1 ...
#
#    -m matrixpds ( tsv, rownames,  matrix with deregulation scores like the median zcores PDS from pathifier, the scores should be under the second column ): 
#        Controls        TCGA_Basal_under_percentile_25_stbl_10
#        hsa03440 Homologous recombination       -0.49621966119327       1.56078456712318 ...
#        hsa05169 Epstein-Barr virus infection   -0.484698560481039      1.51871924292226 ...
#        hsa03013 RNA transport  -0.488858085179727      1.49153097860408 ...
#
#    -t dgetable (tsv file, dataframe with rownames, column names log2FoldChange and padj) 
#        log2FoldChange  padj
#        A2ML1   1.36385896974289        9.19178949837641e-09
#        AACSP1  1.92210465512522        0.0219454396948858
#        AADAC   -2.23896221600004       2.30100976593795e-06
#
#   -l  label (the character string without spaces that you want to label the name of your results)
#   
#   -o  outputfolder ( the path where you want to store your results)
#
#   Example: 
#     Rscript $path/DfPwDefinitions_PDS_DGEtable_2_Df_of_PwGenesPDS.R \
#     -d $path/KEGG_pathways_in_df_genesymbol.tsv \
#     -m /path/to/SubtypevsConstrol_median_PDSz_ordered_matrix_Top20.txt \
#     -t /path/to/table_of_DEGs_lfc2_of0_6_padjof0_05.tsv \
#     -l some_label \
#     -o /path/to/save/your/resultfile

# to do
#   improve -h options with an example

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

##########################################
### Functions defined by the user     ####
##########################################
eliminating_nonDEG_from_Pwdf_GMTformat <- function( input_Pwdf ,  myDfDEG ){
  # This program convert in NAs all the non-DEGs according to all the row names from a given DEG table 
  mydfresult <- input_Pwdf
  for( k in 1: dim( input_Pwdf )[1] ){
    my1rowdf <- input_Pwdf[k,] # extracting the a pathway/row
    DEGs_in_that_pw_positions <- which(  my1rowdf %in% rownames(myDfDEG)) # identifying the positions with DEGs
    DEGs_in_that_pw_contents_vec <-   as.character(my1rowdf[ , DEGs_in_that_pw_positions ]) # extracting the DEGs names
    
    my1rowdf_filledNAs <- my1rowdf # filling the previous in an object to deliver
    my1rowdf_filledNAs[1,] <- rep(NA,dim(my1rowdf)[2])
    my1rowdf_filledNAs[1,DEGs_in_that_pw_positions ] <- DEGs_in_that_pw_contents_vec
    mydfresult[k,]  <- my1rowdf_filledNAs[1,]
  }
  return(  mydfresult ) 
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
parser$add_argument("-d", "--pwdefinitions", type="character", 
                    help="path to the GMT file with with pathway's definitions")
parser$add_argument("-m", "--matrixpds", type="character", 
                    help="file with your PDS (pathifier) matrix, PDS then zscore against controls and median afterwards")
parser$add_argument("-t", "--dgetable", type="character", 
                    help="file with your differential expression table")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="output file where you want to store your correlation plots")

args <- parser$parse_args( )

##########################################################
### Reading the data #####################################
##########################################################
PtwDf_defintions <- read.table( file=args$pwdefinitions , sep= "\t", header=FALSE, quote = "" , row.names = 1 , stringsAsFactors = FALSE )
# Path_to_Ptwydb <-"/home/rmejia/Documents/ABCB1_minimal/Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv" ; PtwDf_defintions <- read.table( file=Path_to_Ptwydb , sep= "\t", header=FALSE, quote = "" , row.names = 1 , stringsAsFactors = FALSE)

Df_your_Dereg_Pw <- read.table( file = args$matrixpds , sep= "\t", header=TRUE, row.names = 1, quote = "" , stringsAsFactors = FALSE )
# Path_to_your_Matrix <- "/home/rmejia/Documents/ABCB1_minimal/Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt" ; Df_your_Dereg_Pw <- read.table( file = Path_to_your_Matrix , sep= "\t", header=TRUE, row.names = 1, quote = "" , stringsAsFactors = FALSE)

DfDEG <- read.table( file=args$dgetable , sep= "\t", header=TRUE, quote = "" , row.names = 1, stringsAsFactors = FALSE )
# Path_your_DEG <-"/home/rmejia/Documents/ABCB1_minimal/Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv" ; DfDEG <- read.table( Path_your_DEG , sep= "\t", header=TRUE, quote = "" , row.names = 1, stringsAsFactors = FALSE)

# Deleting the indicator row (of "Normal") in DEG data frame
NORMALpos <- which( rownames(DfDEG) %in% "NORMAL" )
DfDEG <- DfDEG[ -NORMALpos , ]

Path_of_Results <- normalizePath( args$outputfolder )
dir.create( Path_of_Results, recursive = TRUE)

Label_for_Results <- args$label

########################################################################################################
### Let's identify the elements of your deregulated pathways according with your Pw definition file  ###
########################################################################################################
Coincident_positions <- which( rownames( PtwDf_defintions )  %in% rownames(Df_your_Dereg_Pw )  ) # identifying the Dereg pathways (acording to your dataframe) inside your pw definitions
Pw_definitions_with_deregpw  <- PtwDf_defintions[ Coincident_positions , ] # Dataframe containing the elements of your Dereg pathways

Report1 <- paste( dim( Df_your_Dereg_Pw )[1] , " Pws were annotated according to your Pw definitions from a total of " , dim( Df_your_Dereg_Pw )[1] , " Deregulated pathways"  ) ; print(Report1)

#################################################################################
### Doing deleting all the non-DEG members from the Dereg pathways 
### Deleting the pathways with no DEG (i.e. all NA)
#################################################################################
Sparse_PwDef_df <- eliminating_nonDEG_from_Pwdf_GMTformat( Pw_definitions_with_deregpw  , DfDEG ) # Eliminate non DEGs

a <- Sparse_PwDef_df # Doing the previous Df compact ( erasing the rows with all NAs )
Sparse_PwDef_df_ALLNArows_deleted <- a[ rowSums( is.na( a ) ) != ncol( a ) ,  ]
Report2 <- paste( "only ", dim(Sparse_PwDef_df_ALLNArows_deleted)[1] , "Pws had at least one DEG. From a total of " , dim(a)[1] , "Deregulated pathways"   ) ; print(Report2)

# Converting the DF to a list of characters 
compact_df <- t( Sparse_PwDef_df_ALLNArows_deleted ) # transposing as a preparation for the reshaping
transposed_Df_yourPw_with_DEGs <- as.data.frame( compact_df )
yourPwDef_as_lists <- as.list(  transposed_Df_yourPw_with_DEGs ) # list of characters

compact_list <- lapply( yourPwDef_as_lists , na.omit) # eliminating the NA elements from the lists entries

############################################################################################
######## Once having the compact list let's Build the Dataframe with all the relevant info
#############################################################################################
vector_char_DEGs_in_DegPws <- unlist( compact_list ) # converting the list into a single character vector with names
Deg_Pathways_with_at_least_one_DEG <- rep( names( compact_list ), times= as.numeric( lapply( compact_list , length ) ) ) # Preparing the names of the previous vector (Pathways) to build the BigDf

DF_Pws_DGEs <- cbind( Deg_Pathways_with_at_least_one_DEG,  vector_char_DEGs_in_DegPws) # Let's do the df
colnames( DF_Pws_DGEs) <- c("Pathways","Target")
rownames( DF_Pws_DGEs) <- NULL # I erased it to avoid confusion bc the rows had names like PwA1, PwA2, PwA3, ...

#######################
### Inserting the PDS
#######################
DF_Pws_DGEs_PDS <- as.data.frame(DF_Pws_DGEs)
DF_Pws_DGEs_PDS$MedianPDSzs <- rep(NA,dim(DF_Pws_DGEs_PDS)[1])

for( k in DF_Pws_DGEs_PDS[ , "Pathways"]  ){
  positions <- which( as.character( DF_Pws_DGEs_PDS[ , "Pathways"] )  %in% k )
  DF_Pws_DGEs_PDS$MedianPDSzs[positions] <- Df_your_Dereg_Pw[k,2]
}

############################
### Inserting LogFc & Padj
#############################
DF_Pws_DGEs_PDS_LogFC_Padj <- DF_Pws_DGEs_PDS
DF_Pws_DGEs_PDS_LogFC_Padj$log2FC     <- rep( NA, dim( DF_Pws_DGEs_PDS_LogFC_Padj )[1] )
DF_Pws_DGEs_PDS_LogFC_Padj$DEGpvalue  <- rep( NA, dim( DF_Pws_DGEs_PDS_LogFC_Padj )[1] )

for( k in rownames(DfDEG ) ){
  positions <- which( as.character( DF_Pws_DGEs_PDS_LogFC_Padj[,"Target"])  %in% k )
  DF_Pws_DGEs_PDS_LogFC_Padj$log2FC[positions] <- DfDEG[k,"log2FoldChange"] # The DEG is not here !! (in the DfDEG[700:810,]) why ?? it was supposed to be present by construction!
  DF_Pws_DGEs_PDS_LogFC_Padj$DEGpvalue[positions] <- DfDEG[k,"padj"]
}
Bigdf <- DF_Pws_DGEs_PDS_LogFC_Padj

##########################################################
###### Saving the results
##########################################################
write.table(Bigdf ,
            file=paste0(Path_of_Results,"/","DF_Pw_Target_PDS_LogFC_adjpval",Label_for_Results,".tsv" ),
            sep="\t", quote=FALSE , row.names= FALSE, col.names= TRUE  )
