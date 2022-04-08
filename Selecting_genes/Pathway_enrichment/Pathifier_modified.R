################################################################################
# 2. Filtering and preprocessing of data for Pathifier
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
## MICHAL's  #_# METABRIC flavor
# DATA INPUT
# captura argumentos de la linea de comandos
args <- commandArgs(trailingOnly = TRUE)

###########################################
Path_to_your_Matrix<-args[1] # The path to your matrix
# Path_to_your_Matrix<-c("TCGA_Control_vs_Basal_indicator_10and10.tsv")
gene_sets_path<-args[2] # The path to your pathway's definition file
# gene_sets_path<-c("../Results/KEGG_pathways_in_df_genesymbol_demo.tsv")
Path_of_Code<-args[3] # The path to your code
# Path_of_Code<-c("./")
Path_of_Results<-args[4] # # where do you want to save your results?
# Path_of_Results<-c("../Results/Pathifier/")
Tumour_subtype<-args[5] # Label for your results
# Tumour_subtype<-"Basal"
Stabilizing <- args[6] # Parameter to stabilization of the adjusted curve
# Stabilizing <- 4 # 
Filter_value <- args[7] # Filter low value genes 
# Filter_value=3.75
#######################################################################
###  Some data given by the user should be numeric,         ###########
###  hence we coerce such variables to that mode            ###########
#######################################################################
Stabilizing <- as.numeric( Stabilizing )
Filter_value <- as.numeric( Filter_value )

###############################################################################
### Installing and/or loading required packages
###############################################################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("pathifier")) {
  BiocManager::install("pathifier", ask =FALSE)
  library("pathifier")
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("som")) {
  install.packages("som", dependencies = TRUE)
  library(som)
}

source(paste0(Path_of_Code,"plot_raw_Matrix_png.R"))
source(paste0(Path_of_Code,"plot_raw_Matrix_png_NO_clustering.R"))

####################################################
####   Reading the data expression matrix  #########
####     And preparing it for Pathifier    #########
####################################################
dir.create(Path_of_Results)
M.matrix<-read.table(Path_to_your_Matrix,header=TRUE,row.names = 1)

matrix <- M.matrix
NORMALS <- matrix[1,]

# Calculate MIN_STD
N.matrix <- matrix[-1,as.logical(NORMALS)]
rsd <- apply(N.matrix, 1, sd)
min_std <- quantile(rsd, 0.25)

# Filter low value genes. At least 10% of samples with values over 4
min_exp <- Filter_value #_# Min exp calculated to RDB
over4 <- apply(matrix, 1, function(x) x > min_exp) # Higher than 4
G.over4 <- apply(over4, 2, mean)
G.over4 <- names(G.over4)[G.over4 > 0.1]
matrix <- matrix[G.over4,]

# Filter low variance. Use the genes with more variance (5000 TOP genes)
V <- apply(matrix, 1, var)
V <- sort(V, decreasing = T)
matrix <- matrix[names(V[1:5000]),]

# Set every value bellow 4 to 4  #_# In my case i substitute for 2.9
matrix[matrix < Filter_value] <- Filter_value
genes.m <- rownames(matrix) # Checking genes

# Prepare DATA
matrix <- rbind(NORMALS,matrix)

# Write Table, samples/normals and genes files
write(min_std, file = paste(Path_of_Results,c("M_min_std.txt"),sep = ""))
write.table(matrix, file = paste(Path_of_Results,c("M_BRCA_Pathifier.txt"),sep = ""), quote = F, sep = "\t", 
            col.names = NA)
write.table(x = NORMALS, file = paste(Path_of_Results,c("M_NORMALS.txt"),sep=""), quote = F, sep = "\t", col.names = F)
write(genes.m, file = paste(Path_of_Results,c("genes_M.txt"),sep=""), sep = "\t")

########################################################################
####  Now Pathifier !!                       ###########################
########################################################################
my_min_std<-min_std
my_min_exp<-Filter_value

# Load expression data for PATHIFIER
exp.matrix <- matrix
### My pathifier 

# Loading Genesets annotation 
# NOTE: MUST contain SAME number of elements in all rows!!

gene_sets <- as.matrix(read.delim(file = gene_sets_path,header = F, sep = "\t", as.is = T))

#  Generate a list that contains genes in genesets
gs <- list()
for (i in 1:nrow(gene_sets)){
  a <- as.vector(gene_sets[i,2:ncol(gene_sets)])
  a <- unique(na.omit(a))
  a <- matrix(a, ncol = 1)
  gs[[length(gs)+1]] <- a
  rm(a,i)
}

# Generate a list that contains the names of the genesets used
pathwaynames <- as.list(gene_sets[,1])

# Generate a list that contains the previos two lists: genesets and their names
PATHWAYS <- list()
PATHWAYS$gs <- gs
PATHWAYS$pathwaynames <- pathwaynames

# Extract information from binary phenotypes. 1 = Normal, 0 = Tumor
normals <- as.vector(as.logical(exp.matrix[1,]))
data <- as.matrix(exp.matrix[-1, ])
dimnames(data) <- NULL
allgenes <- as.vector(rownames(exp.matrix))

# Generate a list that contains previous data: gene expression, normal status,
# and name of genes
DATASET <- list(); DATASET$allgenes <- allgenes[-1]
DATASET$normals <- normals
DATASET$data <- data

# Make some timing ;)
ptm <- proc.time()

# Run Pathifier
PDS<-quantify_pathways_deregulation(DATASET$data, DATASET$allgenes,
                                    PATHWAYS$gs,
                                    PATHWAYS$pathwaynames,
                                    DATASET$normals, attempts = Stabilizing,
                                    logfile="logfile.txt", min_std = my_min_std, min_exp =my_min_exp)

#Stop the clock
PDS_time<-proc.time() - ptm

# Remove unnecesary data
#rm(gene_sets, data, exp.matrix, allgenes, DATASET, PATHWAYS)

# Saving data 

write(PDS_time, file = paste(Path_of_Results,Tumour_subtype,c("PDS_time.txt"),sep = ""))
save.image(file = paste(Path_of_Results,Tumour_subtype,c("PDS.RData"),sep=""))

###############################################################################
####################   Plotting  Raw PDS ######################################
###############################################################################
color_labels_vector<-as.character(matrix[1,])
color_labels_vector<-gsub("1","green",color_labels_vector)
color_labels_vector<-gsub("0","red",color_labels_vector)

###############################################################################
## Load Pathifier results and turn into a matrix 
###############################################################################

#load(RData_Directory)
PDSmatrix <- mapply(FUN = c, PDS$scores)
PDSmatrix <- t(PDSmatrix)
colnames(PDSmatrix)<- colnames(matrix)
###############################################################################
## Generating and assigning labels for pathways used in the analysis
###############################################################################

###############################################################################
## Creating Custom Palette
###############################################################################

# creates a own color palette passing from blue, green yellow to dark red
my_palette <- colorRampPalette(c("blue", "cyan", "chartreuse1", "yellow",
                                 "red", "firebrick4"))(n = 1000)

# (optional) defines the color breaks manually for a "skewed" color transition
# col_breaks = c(seq(-1,0,length=100),
#                seq(0,0.8,length=100),
#                seq(0.8,1,length=100))

###############################################################################
## Clustering Methods
###############################################################################

# If you want to change the default clustering method (complete linkage method
# with Euclidean distance measure), this can be done as follows: For non-square
# matrix, we can define the distance and cluster based on our matrix data by

row.distance = dist(PDSmatrix, method = "euclidean")
row.cluster = hclust(row.distance, method = "complete")

col.distance = dist(t(PDSmatrix), method = "euclidean")
col.cluster = hclust(col.distance, method = "complete")

# Arguments for the dist() function are: euclidean (default), maximum, canberra,
# binary, minkowski, manhattan

# And arguments for hclust(): complete (default), single, average, mcquitty,
# median, centroid, ward.

# NOTE that for non-square matrices you have to define the distance and cluster
# for both row and column dendrograms separately.
# Otherwise you will get a not so pleasant Error in:
# x[rowInd, colInd] : subscript out of bounds.

###############################################################################
## Assign Colors to Columns 
###############################################################################

###############################################################################
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################

png(paste(Path_of_Results,Tumour_subtype,c("_heatmap_raw_PDS_samples_Labels.png"),sep=""), # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 5)        # font size

heatmap.2(PDSmatrix,
          main = paste0("Raw_PDS_",Tumour_subtype),          # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(10,21),     # widens margins around plot
          col = my_palette,        # use on color palette defined earlier
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize = 0.8,           # size of color key
          #Additional Options
          ## Color labeling columns (Opt. RowSideColors for rows)
          ColSideColors =color_labels_vector ,       # Grouping col-samples into two different
          
          #    breaks= col_breaks,  # enable color transition at specified limits
          dendrogram= "col",   # only draw a column dendrogram (opt. "row")
          cexRow = 0.2 + 1/log10(dim(PDSmatrix)[1]),
          cexCol= 0.2 + 1/log10(dim(PDSmatrix)[2])
          
          ## Legend for ColumnSide color labeling
          # par(lend = 1)           # square line ends for the color legend
          # legend("topright",      # location of the legend on the heatmap plot
          #        legend = c("Normals", "Tumors"), # category labels
          #        col = c("dodgerblue", "firebrick1"),  # color key
          #        lty= 1,          # line style
          #        lwd = 5, unit    # line width
)

dev.off()               # close the PNG device


################# Save the PDS_Matrix in .txt ############


colnames(PDSmatrix)<-colnames(matrix)
matrix_txt_path<-paste(Path_of_Results,Tumour_subtype,c("_rawPDSmatrix.txt"),sep="")
write.table(PDSmatrix,file=matrix_txt_path,quote=FALSE,sep="\t",row.names = TRUE,col.names = TRUE)




#########################################################
#############      PDS in z-score          ##############
#########################################################

# calculating z-scores
PDSmatrix_zscores<-som::normalize(PDSmatrix, byrow = TRUE)
# saving it
matrix_txt_path<-paste(Path_of_Results,Tumour_subtype,c("_PDSz_matrix.txt"),sep="")
write.table(PDSmatrix_zscores,file=matrix_txt_path,quote=FALSE,sep="\t",row.names = TRUE,col.names = TRUE)
#ploting
plot_raw_Matrix_png(PDSmatrix_zscores,paste0("PDSz_subtype_",Tumour_subtype),paste(Path_of_Results,Tumour_subtype,c("_heatmap_PDSz_custering.png"),sep=""),color_labels_vector)


######################################################################
#############      medians of PDS and PDSz  by Pathway  ##############
######################################################################

median_dfz<-as.matrix(data.frame(apply(PDSmatrix_zscores[,matrix[1,]==1],1,median),apply(PDSmatrix_zscores[,matrix[1,]==0],1,median)))
colnames(median_dfz)<-c("Controls",Tumour_subtype)
median_df<-as.matrix(data.frame(apply(PDSmatrix[,matrix[1,]==1],1,median),apply(PDSmatrix[,matrix[1,]==0],1,median)))
colnames(median_df)<-c("Controls",Tumour_subtype)



png(paste(Path_of_Results,Tumour_subtype,c("Boxplot_of_medians_PDS.png"),sep=""), # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 5)
boxplot(median_df,main=paste0("Boxplot of medians PDS"))
dev.off()
png(paste(Path_of_Results,Tumour_subtype,c("Boxplot_of_medians_PDSz.png"),sep=""), width = 6 * 500, height = 6 * 400,   units = "px", res = 300, pointsize = 5)
boxplot(median_dfz,main=paste0("Boxplot of medians PDSz"))
dev.off()
png(paste(Path_of_Results,Tumour_subtype,c("PDS_by_Pathways.png"),sep=""), width = 6 * 500, height = 6 * 400,   units = "px", res = 300, pointsize = 5)
boxplot(t(PDSmatrix),main=paste0("PDS by Pathway"))
dev.off()
png(paste(Path_of_Results,Tumour_subtype,c("PDSz_by_Pathways.png"),sep=""), width = 6 * 500, height = 6 * 400,   units = "px", res = 300, pointsize = 5)
boxplot(t(PDSmatrix_zscores),main=paste0("PDSz by Pathway"))
dev.off()

# Ploting the Entire matrix without clustering
# plot_raw_Matrix_png_NO_clustering(PDSmatrix,Tumour_subtype,paste(Path_of_Results,c("heatmap_Raw_PDS_samples_Labels_No_clustering_rows.png"),sep=""))
# plot_raw_Matrix_png_NO_clustering(PDSmatrix_zscores,Tumour_subtype,paste(Path_of_Results,c("heatmap_PDSz_samples_Labels_No_clustering_rows.png"),sep=""))

median_dfz_ordered <- median_dfz[order(median_dfz[,2],decreasing = TRUE),]
# plot_raw_Matrix_png_NO_clustering(median_dfz_ordered,Tumour_subtype,paste(Path_of_Results,Tumour_subtype,c("_median_PDSz_No_clustering_rows.png"),sep=""))
matrix_txt_path <- paste(Path_of_Results,Tumour_subtype,c("_median_PDSz_ordered_matrix.txt"),sep="")
write.table(median_dfz_ordered,file=matrix_txt_path,quote=FALSE,sep="\t",row.names = TRUE,col.names = TRUE)

# Plotting the top 20
# plot_raw_Matrix_png_NO_clustering(median_dfz_ordered[1:20,],Tumour_subtype,paste(Path_of_Results,Tumour_subtype,c("_TOP_20_median_PDSz_No_clustering_rows.png"),sep=""))

matrix_txt_path<-paste(Path_of_Results,Tumour_subtype,c("_median_PDSz_ordered_matrix_Top20.txt"),sep="")
write.table(median_dfz_ordered[1:20,],file=matrix_txt_path,quote=FALSE,sep="\t",row.names = TRUE,col.names = TRUE)
