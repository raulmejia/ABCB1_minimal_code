########################################################################
## This script get a dataframe (describing phenotypes) and return the 
## data frame with a "Triple negative column" added
# Writed by Raúl Alejandro Mejía Pedroza                              ##
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
# Nota la categoria de referencia se llama "Control" 
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
NORMAL <- Exp_Mat_bk[ grep("NORMAL",rownames(Exp_Mat_bk)),]
Exp_Mat_bk <- Exp_Mat_bk[- grep("NORMAL",rownames(Exp_Mat_bk)),]
Exp_Mat_bk <- log2(Exp_Mat_bk +1)
Exp_Mat_bk <- rbind(NORMAL , Exp_Mat_bk)

########################################################################
#######   Extracting the propper labels ################################
########################################################################
positions <- which( rownames(MyLabels) %in% colnames(Exp_Mat_bk) )
Adjusted_Labels <- data.frame(MyLabels[ positions,])
rownames(Adjusted_Labels) <- rownames(MyLabels)[ positions]

#######################################################################
############   Building the DDS DESeqDataSet      #####################
#######################################################################
colnames(Adjusted_Labels)[1] <- "condition"
Adjusted_Labels$condition <- relevel(Adjusted_Labels$condition, ref="Control")
My_dds <- DESeqDataSetFromMatrix(countData = round( Exp_Mat_bk ), colData = Adjusted_Labels, design = ~ condition) 
colnames( Exp_Mat_bk) == rownames(Adjusted_Labels)
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

