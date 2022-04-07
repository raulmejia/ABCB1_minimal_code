################################################################################
# This script retrieves the database Basal subtype's Breast Cancer samples
# giving as a result the expression matrix and the Lehman subtype of the samples 
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
################################################################################
# Installing and loading the required libraries
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("breastCancerMAINZ")) {
  BiocManager::install("breastCancerMAINZ", ask =FALSE)
  library("breastCancerMAINZ")
}
if (!require("STROMA4")) {
  BiocManager::install("STROMA4", ask =FALSE)
  library("STROMA4")
}
if (!require("Biobase")) {
  BiocManager::install("Biobase", ask =FALSE)
  library("Biobase")
}
if (!require("genefilter")) {
  BiocManager::install("genefilter", ask =FALSE)
  library("genefilter")
}
if (!require("robustbase")) {
  BiocManager::install("robustbase", ask =FALSE)
  library("robustbase")
}

###########################################
# DATA INPUT given by the user
###########################################
args <- commandArgs(trailingOnly = TRUE)
Path_to_save_your_ExpMatrix <-args[1] # The path to save your expression matrix
Path_to_save_your_Lehmann_subtypes <-args[2] # The path to save your matrix of Lehmann subtypes
Path_of_Code <-args[3] # The path to your code, scripts and so on ...
Label_for_Results <-args[4] # Label for your results

## Running in the command line (you can uncomment selecting the following 5 lines and then Crtl+Shift+C
# Rscript $path/to/script/getting_MAINZ_DB_V1.0.R \
# /Path/to/results/extracted_basals_no_controls_matrix/Splited/Submatrices_ohne_controls/ \
# /Path/to/Results-from-STROMA-algorithm/Lehmann-STROMA4/MAINZ/ \
# /Folder/with/the/code \
# Label

# Note: If you use windows and you have troubles with your paths, please try out this tool: # Path_of_x <-choose.dir()

##########################################
## Creating the folders 
##########################################
dir.create(Path_to_save_your_ExpMatrix, recursive = TRUE)
Path_to_save_your_ExpMatrix <- normalizePath( Path_to_save_your_ExpMatrix)
dir.create(Path_to_save_your_Lehmann_subtypes , recursive = TRUE)
Path_to_save_your_Lehmann_subtypes <- normalizePath( Path_to_save_your_Lehmann_subtypes)

###########################################
# Loading the data
###########################################
data(mainz, package='breastCancerMAINZ')

EmMainz<- cbind( fData(mainz)[,c("probe","Gene.symbol")] , Biobase::exprs( mainz ))
EmMainz[1:5,1:5]
rownames( EmMainz ) <- make.names( EmMainz[, 2], unique=TRUE) ## putting the second column (GenSym) as rownames

EmMainz$Gene.symbol <- as.factor(EmMainz$Gene.symbol)
EmMainz_splitted <- split( EmMainz[,-c(1,2)], EmMainz$Gene.symbol)
make_me_a_matrix <- function(someDF){
  someDF <- as.matrix(someDF)
}

lapply(EmMainz_splitted, dim)
EmMainz_splitted_matrices <- lapply(EmMainz_splitted, make_me_a_matrix)
lapply(EmMainz_splitted_matrices, dim)
EmMainz_splitted_matrices_collapsed <- lapply(EmMainz_splitted_matrices, robustbase::colMedians)
EmMainz_splitted_matrices_collapsed[[3]][1:5]
EmMainz_collapsed <- data.frame(matrix(unlist(EmMainz_splitted_matrices_collapsed), nrow= length(EmMainz_splitted_matrices_collapsed), byrow=T),stringsAsFactors=FALSE)
rownames(EmMainz_collapsed) <- names( EmMainz_splitted_matrices_collapsed)
colnames(EmMainz_collapsed) <- names(EmMainz_splitted_matrices_collapsed[[1]])
EmMainz_collapsed[1:5,1:5]
###
EmMainz <- EmMainz[, -c(1,2)]
EmMainz[1:5,1:5]

write.table( EmMainz ,           file = paste0( Path_to_save_your_ExpMatrix ,"/", "Subexpression_matrix_Basal_from_MAINZ_HCGS_.tsv" ), sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE )
write.table( EmMainz_collapsed , file = paste0( Path_to_save_your_ExpMatrix , "/", "Subexpression_matrix_Basal_from_MAINZ_HCGS_collapsed.tsv" ) , sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE )
write.table( mainz ,             file = paste0( Path_to_save_your_ExpMatrix , "/", "Subexpression_matrix_Basal_from_MAINZ_atIDs_.tsv" ), sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE )



###########################################
# Processing
###########################################
all.properties <- assign.properties(ESet=mainz, geneID.column="Gene.symbol",
                                    genelists=c("Stroma4", "TNBCType"), n=10, mc.cores=1)
allp <- all.properties
# Creating the folders to save the data
dir.create(Path_to_save_your_ExpMatrix ,recursive = TRUE) ; dir.create(Path_to_save_your_Lehmann_subtypes , recursive = TRUE)

###################################
##  converting the results from STROMA4 ("high", "intermediate", "low") to  numerical (-1, 0 , 1)
##################################
numericalallp <- pData(allp)[, c("D.stroma.property", "MSL.property","M.property","B.stroma.property","LAR.property","IM.property","T.stroma.property","BL1.property","E.stroma.property","BL2.property")]
char2num <- function(row){
  row <-gsub("low" ,-1,row)
  row <-gsub("intermediate" , 0 ,row)
  row <-gsub("high" , 1 ,row)
  return(as.numeric(row))
}
numericalallp <- apply(numericalallp,1,char2num)
head(numericalallp)
numericalallp <- t(numericalallp)
colnames(numericalallp) <- colnames( pData(allp)[, c("D.stroma.property", "MSL.property","M.property","B.stroma.property","LAR.property","IM.property","T.stroma.property","BL1.property","E.stroma.property","BL2.property")] )

write.table( numericalallp , file = paste0(  Path_to_save_your_Lehmann_subtypes ,"Lehmann's_Subt_and_properties_Numerical_", Label_for_Results,".tsv" ), sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE )



