# This program Split an expression matrix (EM) in submatrices according to a given data frame that
# contains the column names of the EM in the same order that the labels. 

#################################################
### Loading and installing required libraries
#################################################
args <- commandArgs(trailingOnly = TRUE)

#if (!require("splitstackshape")) {
#  install.packages("splitstackshape", ask =FALSE)
#  library("splitstackshape")
#}

###########################################
### Data given by the user
###########################################
# Path_to_your_Matrix<-c("../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt")
Path_to_your_Matrix <-args[1]
# Path_to_Labels<-c("../Data/Labels_Controls_and_Normal_separated_TCGA.txt")
Path_to_Labels <-args[2]
# Path_of_Results<-c("../Results/Splited/Matrices/")
Path_of_Results <-args[3]
# Title_of_Results <- c("TCGA")
Title_of_Results <-args[4]
############################################
### Reading the data
############################################
MyEMatrix <- read.table(Path_to_your_Matrix,header=TRUE, sep= "\t", quote = "")
Labels<-read.table(Path_to_Labels)

############################################
### Quality control 
############################################
if (length(colnames(MyEMatrix)) / ( sum(colnames(MyEMatrix) == rownames(Labels)) ) == 1 ){
  print("Your expression matrix and Labels have the same names in the same order")
} else{
  print("Your expression matrix and Labels DOESN'T have the same names in the same order")
} # Checking consistency between the give data

############################################
### Splitting the original data frame in subdataframes
############################################
### Preparing the data for the splitting function
EMLabeled <- cbind(as.character(Labels$Labels), as.data.frame(t(MyEMatrix) , StringAsFactors = FALSE) ) # tuvimos que convertirlo a Df's porque sólo corta data frames
colnames(EMLabeled)[1] <- "Labels" 
list_of_submatrices <- split( EMLabeled , EMLabeled$Labels )

### Deleting the first column "Labels"
for( k in 1:length(list_of_submatrices)){
  DeleteThisPositions <- which(colnames(list_of_submatrices[[k]]) %in% c("Labels"))
  list_of_submatrices[[k]] <- list_of_submatrices[[k]][, -DeleteThisPositions]
}

### Transposing matrices to their original orientation
list_of_submatrices <- lapply( list_of_submatrices, t)

###########################################
######### Save the results
###########################################
dir.create( Path_of_Results, recursive = TRUE)
saveRDS(list_of_submatrices , file = paste0(Path_of_Results,"List_of_sub_expression_matrices_of_",Title_of_Results,".RDS" ))
for( k in names(list_of_submatrices)  ){
  write.table(list_of_submatrices[[k]], col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE,
              file = paste0(Path_of_Results,"Subexpression_matrix_",k,"_from_",Title_of_Results,"_.tsv" ))
}
