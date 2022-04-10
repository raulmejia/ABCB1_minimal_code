########################################################################
## This script get a dataframe (describing phenotypes) and return the 
## data frame with a "Triple negative column" added
# Writed by Raúl Alejandro Mejía Pedroza                              ##
########################################################################
#
# To do:
#     Implement DESeq2, test the DEG results against gold standard
#
# Your Input matrix:
#     Tab separated
#     Indicator Row: "NORMAL 1 1 1 0 0 ..."
#     Shouldn't be log2 transformed (the program does it)
#     Example: 
#         SampleA SampleB SampleC...  
#         NORMAL  1       1       0 ....
#         A1BG    154.473 109.3484        15.7841 19.0722 41.5445 ...
#         A1BG-AS1        171.5365        131.197 33.7322 33.5206 69.545 ...

# Your input Labels data Frame:
# Note: The reference category is called "Control"
#     Tab separated:
#     "Labels"
#     "SampleA" "Control"
#     "SampleB" "Control"
#     "SampleC" "SomeCondition"

# The folder were your code is located /Folder_Code/Selecting_genes/DEG/

# Folder to save the results /Folder_to_Results/DEG/TCGA/log2only/

# Some label to your results

# LogFC

# P.adj

# Number of cores to use

# Example of run:
#  Rscript $Folder_Code/DESeq2_with_log2_transformed_data_only.R \
#  /Folder_to_Results/Your_input_matrix.tsv \
#  /Folder_with_Data/Labels_Controls_and_Normal_separated_TCGA.txt \
#  /Folder_Code/Selecting_genes/DEG/ \
#  /Folder_to_Results/DEG/TCGA/log2only/ \
#  Somelabel_to_your_results_only_log2transformed_lfc2_of0_6_padjof0_05 \
#  0.6 \
#  0.05 \
#  7 # Cores available

########################################################################
#######   Data selected by te user #####################################
########################################################################
args <- commandArgs(trailingOnly = TRUE)
exp_mat_path <-args[1]
# exp_mat_path <-c("../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv")

Labels_path <-args[2]
# Labels_path <-c("../Data/Labels_Controls_and_Normal_separated_TCGA.txt")

Path_of_Code<-args[3]
# Path_of_Code<-c("./")

results_path <-args[4]
# results_path <-c("../Results/DEG/TCGA/log2only/")

Label_for_results <-args[5]
# Label_for_results <- c("_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed")

mylfctreshold <-args[6]
# mylfctreshold <- 0.6

myp.adj <-args[7]
# myp.adj <- 0.05

mycoresavaiable <-args[8]
# mycoresavaiable <-7

################################################
######  Coercing to numeric ####################
################################################
mylfctreshold <- as.numeric(mylfctreshold)
myp.adj <- as.numeric(myp.adj)
mycoresavaiable <- as.numeric(mycoresavaiable)
################################################
######  Libraries needed    ####################
################################################
if (!require("BiocManager")) {
  biocLite("BiocManager")
  library(BiocManager)}
if (!require("DESeq2")) {
  BiocManager::install("DESeq2")
  library(DESeq2)}
if (!require("BiocParallel")) {
  BiocManager::install("BiocParallel")
  library(BiocParallel)}
if (!require("DEFormats")) {
  BiocManager::install("DEFormats")
  library(DEFormats)}

dir.create(results_path,recursive = TRUE)
########################################################################
#######   Loading the data          ####################################
########################################################################
Exp_Mat <- read.table(exp_mat_path ,header=TRUE , sep= "\t")
Exp_Mat <- as.matrix(Exp_Mat)
MyLabels <-read.table(Labels_path)

#######################################################################
############   log2 transforming your matrix              #############
#######################################################################
Exp_Mat_bk <- Exp_Mat 
NORMAL <- Exp_Mat_bk[ grep("NORMAL",rownames(Exp_Mat_bk)),] # coping your indicator row
Exp_Mat_bk <- Exp_Mat_bk[- grep("NORMAL",rownames(Exp_Mat_bk)),] # deleting your indicator row
Exp_Mat_bk <- log2(Exp_Mat_bk +1) # log2 transforming
Exp_Mat_bk <- rbind(NORMAL , Exp_Mat_bk) # pasting your indicator row again

########################################################################
#######   Extracting the proper labels ################################
########################################################################
positions <- which( rownames(MyLabels) %in% colnames(Exp_Mat_bk) ) # Extracting just the labels related to the colnames of your input matrix
Adjusted_Labels <- data.frame(MyLabels[ positions,]) # eyo
rownames(Adjusted_Labels) <- rownames(MyLabels)[ positions] # Labels suited to your inputmatrix

#######################################################################
############   Building the DDS DESeqDataSet      #####################
#######################################################################
colnames(Adjusted_Labels)[1] <- "condition"
Adjusted_Labels$condition <- as.factor(Adjusted_Labels$condition)
Adjusted_Labels$condition <- relevel(Adjusted_Labels$condition, ref="Control")
My_dds <- DESeqDataSetFromMatrix(countData = round( Exp_Mat_bk ), colData = Adjusted_Labels, design = ~ condition) 

if( all(colnames( Exp_Mat_bk) == rownames(Adjusted_Labels)) == TRUE){
  print("Colnames of your inputmatrix and their adjusted levels match")
}
if( all(colnames( Exp_Mat_bk) == rownames(Adjusted_Labels)) != TRUE){
  print("Colnames of your inputmatrix and their adjusted levels doesn't match check again your inputs")
  break()
}

#######################################################################
############   VST of your matrix               #####################
#######################################################################
saveRDS(My_dds,file=paste0(results_path,"dds_log2_only_matrix_from",Label_for_results,"_.RDS"))
# dists <- dist(t(assay(My_dds))) # let's see the grouping of the samples
# plot(hclust(dists))

#######################################################################
############   Now the DGE function               #####################
#######################################################################
time1<-proc.time()

  register(MulticoreParam(mycoresavaiable))
  My_DESeq <- DESeq( My_dds,parallel = TRUE)
  print(" My_DESeq <- DESeq( My_dds,parallel = TRUE) finished")
  register(MulticoreParam(mycoresavaiable))
  results_DESeq <-results( My_DESeq,parallel = TRUE)
  print("results_DESeq <-results( My_DESeq,parallel = TRUE) finished")
time2 = proc.time()-time1


# Saving the results
save(My_DESeq,file=paste0(results_path,"DESeq_object_of_",Label_for_results,".RData"))
save(results_DESeq,file=paste0(results_path,"DESeq_results",Label_for_results,".RData"))

register(MulticoreParam(mycoresavaiable))
lfc1_results_DESeq <-results(My_DESeq,parallel = TRUE, lfcThreshold= mylfctreshold)

time3 = proc.time()-time1

# Saving the results and timming
save(lfc1_results_DESeq,file=paste0(results_path,"lfc1_results_DESeq",Label_for_results,".RData"))
write.table(c(time2,time3),file=paste0(results_path,Label_for_results,"Timing_runfunctions.txt"))

  pminus10_3<-grepl("TRUE",lfc1_results_DESeq$padj < myp.adj)
  padj10_3_lfc1_results_DESeq <-lfc1_results_DESeq[pminus10_3,]

#saving the results
save(padj10_3_lfc1_results_DESeq,file=paste0(results_path,"padj10_3_lfc1_results_DESeq",Label_for_results,".RData"))
write.table(padj10_3_lfc1_results_DESeq[,c("log2FoldChange","padj")], 
            file=paste0(results_path,"padj10_3_lfc1_results_DESeq",Label_for_results,".tsv"),
            sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE
            )

