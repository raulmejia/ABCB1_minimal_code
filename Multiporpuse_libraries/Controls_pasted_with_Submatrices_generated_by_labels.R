# This program Split a expression matrix (EM) in submatrices that include
# the controls submatrix with the submatrix of certain label (in that order)
# according to a given data frame that contains the column names of the EM in the same order that the labels. 

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
  print("Ok: Your expression matrix and Labels have the same names in the same order")
} else{
  print("ERROR: Your expression matrix and Labels DOESN'T have the same names in the same order")
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
#########  Pasting the submatrix of "Controls" to the other submatrices
###########################################
Position_of_controls <- which(names(list_of_submatrices) %in% "Control")
if (  (names(list_of_submatrices)[Position_of_controls] == "Control") == TRUE ){
  print("Ok: You have a Label of -Control- that will be used for furhter processing")
} else{
  print("ERROR: You DON'T have a Label of -Control- (Exact match)")
} # Checking consistency between the give data

### Pasting the submatrices with the submatrix of Control
list_of_submatrices_pasted_with_Control <- list()
for( k in setdiff(names(list_of_submatrices),"Control")  ){
  list_of_submatrices_pasted_with_Control[[k]] <- cbind(list_of_submatrices[["Control"]],list_of_submatrices[[k]])
}

###########################################
######### Save the results
###########################################
dir.create( Path_of_Results, recursive = TRUE)
saveRDS(list_of_submatrices_pasted_with_Control , file = paste0(Path_of_Results,"List_of_Control_with_subexpression_matrix_of_",Title_of_Results,".RDS" ))
for( k in names(list_of_submatrices_pasted_with_Control)  ){
  write.table( list_of_submatrices_pasted_with_Control[[k]], col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE,
              file = paste0(Path_of_Results,"Control_with_subexpression_matrix_of_",k,"_from_",Title_of_Results,"_.tsv" ))
}
