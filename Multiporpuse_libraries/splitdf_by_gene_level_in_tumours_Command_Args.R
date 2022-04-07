# This program Split a expression matrix (EM) in four submatrices according to
# the expression level of a particular gene (above the mean, under the mean)
# over quantile 75 , below quantile 25 and the paste each of such submatrices after a submatrix of controls
# It requires a indicator row to depict which columns are healthies (1 by convention) and disease (0)

#################################################
### Loading and installing required libraries
#################################################

#################################################
### Defining the required functions
#################################################

Split_by_thisgene_level<- function( mymatrix, gene){
  # This funtion receives only the matrix of tumours and split it by that gene value 
  mymatrix_25th_top_high <- mymatrix[ , mymatrix[gene,] > quantile( as.matrix(mymatrix[gene,]) , .75) ]
  mymatrix_high <- mymatrix[, mymatrix[gene,] > quantile( as.matrix(mymatrix[gene,]) , .50) ]
  # mymatrix_high <- mymatrix[, mymatrix[gene,] > mean(as.matrix(mymatrix[gene,]))]
  #mymatrix_low <- mymatrix[,mymatrix[gene,] < mean(as.matrix( mymatrix[gene,]))]
  mymatrix_low <- mymatrix[,mymatrix[gene,] < quantile( as.matrix(mymatrix[gene,]) , .50) ]
  mymatrix_25th_top_low <- mymatrix[ , mymatrix[gene,] < quantile(as.matrix(mymatrix[gene,]) , .25) ]
  result<- list(mymatrix_25th_top_high, mymatrix_high, mymatrix_low, mymatrix_25th_top_low)
  names(result) <- c("25th_top_high","high","low","25th_top_low")
  return(result)
}

merge_with_healthem <- function(list_splited, healthem){
  #This funcion paste the previously splitted tumours with the "Healthies"
  list_splited[["25th_top_high"]] <- cbind(healthem, list_splited[["25th_top_high"]]  )
  list_splited[["high"]] <- cbind(healthem, list_splited[["high"]] )
  list_splited[["25th_top_low"]] <- cbind(healthem, list_splited[["25th_top_low"]] )
  list_splited[["low"]] <- cbind(healthem, list_splited[["low"]] )
  return(list_splited)
}

splitdf_by_gene_level_in_tumours <- function(df, healthy_positions, tumour_positions,mygene){
  #This function split a data frame in health and disease, and create sub data frames of the
  # disease for example the 4th df has eliminated all the patients (columns) if 
  # their level of the choosed gene is below than the mean, below than 25 th percentil
  # then that submatrix will be pasted to the Healthy df
  # the function returns a list whose elements are each data frame  according that condition
  df_Health <- df[,healthy_positions]
  df_Tumour <- df[,tumour_positions]
  list_of_splited_df_Tumour_by_thisgene_levels <- Split_by_thisgene_level(df_Tumour,mygene)
  list_of_splited_df_by_thisgene_levels <- merge_with_healthem(list_of_splited_df_Tumour_by_thisgene_levels, df_Health)
}

###########################################
### Data given by the user
###########################################
args <- commandArgs(trailingOnly = TRUE)
# Path_to_your_Matrix<-c("../Results/Splited/SubMatrices_with_controls/Control_with_subexpression_matrix_of_Basal_from_TCGA_.tsv")
Path_to_your_Matrix <-args[1]
# Mygene <-c("ABCB1")
Mygene <-args[2]
# Path_of_Results <- c("../Results/Matrices_splited_by_gene/ABCB1/")
Path_of_Results <- args[3]
# Title_of_Results <- c("TCGA_Basal")
Title_of_Results <-args[4]
############################################
### Reading the data
############################################
MyEMatrix <- read.table(Path_to_your_Matrix,header=TRUE, sep= "\t", quote = "")
 
dir.create(Path_of_Results, recursive = TRUE)
############################################
### Reading the data
############################################¨
NumberofControls <- sum(MyEMatrix[1,] == 1)
NumberofCases <- sum(MyEMatrix[1,] == 0)
list_TCGA_splited_by_ABCB1_levels_in_tumours <- splitdf_by_gene_level_in_tumours(MyEMatrix , 1:NumberofControls, 113:(NumberofControls+NumberofCases), "ABCB1")

############################################
### Saving the data
############################################
for(  k in 1: length( list_TCGA_splited_by_ABCB1_levels_in_tumours )){
  write.table( list_TCGA_splited_by_ABCB1_levels_in_tumours[[k]]  ,
              file = paste0( Path_of_Results ,Title_of_Results,"_splited_by_the_expression_of_the_gene_", Mygene,"_",names(list_TCGA_splited_by_ABCB1_levels_in_tumours)[k],".tsv" ),
              sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE  )    
}  