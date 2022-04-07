
######## This Function creates a Heatmap of rawvalues #################################
### Just give the 3 arguemnts:
##### MyMatrix is the matrix that you what make the heatmap ###########################
##### PlotTitle is the title that you want to display in your plot
##### Pathtosave is the direction including the final name and extension for your graph
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

plot_raw_Matrix_png <- function(MyMatrix,PlotTitle,Pathtosave,mycolorside){
  # This function takes a matrix, convert it in z-scores by row; plot and save  both in several formats. png and .svg 
  

###############################################################################
## Creating Custom Palette
###############################################################################

# creates a own color palette passing from blue, green yellow to dark red
my_palette <- colorRampPalette(c("blue", "cyan", "chartreuse1", "yellow",
                                 "red", "firebrick4"))(n = 100)

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

row.distance = dist(MyMatrix, method = "euclidean")
row.cluster = hclust(row.distance, method = "complete")

col.distance = dist(t(MyMatrix), method = "euclidean")
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
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################

png(Pathtosave, # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 5)        # font size

heatmap.2(MyMatrix,
          main = PlotTitle,          # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(10,21),     # widens margins around plot
          col = my_palette,        # use on color palette defined earlier
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize = 0.8,           # size of color key
          #Additional Options
          ## Color labeling columns (Opt. RowSideColors for rows)
         ColSideColors = mycolorside ,       # Grouping col-samples into two different
          
          #    breaks= col_breaks,  # enable color transition at specified limits
          dendrogram= "col",   # only draw a column dendrogram (opt. "row")
          cexRow = 0.2 + 1/log10(dim(MyMatrix)[1]),
          cexCol= 0.2 + 1/log10(dim(MyMatrix)[2])
          
          ## Legend for ColumnSide color labeling
          # par(lend = 1)           # square line ends for the color legend
          # legend("topright",      # location of the legend on the heatmap plot
          #        legend = c("Normals", "Tumors"), # category labels
          #        col = c("dodgerblue", "firebrick1"),  # color key
          #        lty= 1,          # line style
          #        lwd = 5, unit    # line width
)

dev.off()               # close the PNG device


}