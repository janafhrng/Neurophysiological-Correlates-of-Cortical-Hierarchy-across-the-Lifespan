# library to read matlab data formats into R
# library(raveio)
library(tidyverse)
library(pheatmap)
library(dendextend)
library(scico)
library(RColorBrewer)
library(Cairo)
library(hdf5r)
install.packages('rmatio')

filepath = getwd()

##### LME MAIN EFFECT CORTICAL HIERARCHY #####
nfeat = 50

# Load lme t values
tval <-  read.csv("G:/JF/HCTSA feature gradients/CamCan/HierarchyCorrelates/DATA/lme_out_age_effect_hierarchy.csv", header = F)
tval_h <- tval$V4
tval_sort = sort(abs(tval_h),decreasing = T)[1:50]
index = order(abs(tval_h),decreasing = T)

# load features and numeric feature matrix and sort them according to this index
h5 <- H5File$new("G:/JF/HCTSA feature gradients/CamCan/big DATA files/mat_z.mat", mode = "r")
#mat_z <- read_mat(paste(filepath,"/DATA/mat_z.mat" ,sep=""),ram=FALSE)
mat_z = h5[["mat_z"]]$read()
# df <- mat_z[[1]]
mat_z_sort = mat_z[,,index]
mat_z_sort <- mat_z_sort[,,1:50]

table <- read.csv('G:/JF/HCTSA feature gradients/CamCan/HierarchyCorrelates/DATA/features_z.csv', header = F)
table_sort = table[index,]
features_sort = table_sort[1:50,]

# reshape mat into a long matrix where participants*parcels x features
dim(mat_z_sort) <- c(200 * 350 , nfeat)
mat_z_corr <- abs(cor(mat_z_sort,method = "spearman"))

# clustered heatmap
cols= str_c(features_sort[,3])
#cols= features_sort[,4]
rownames(mat_z_corr) <- cols
colnames(mat_z_corr) <- cols

annotation_col = data.frame(gsub(",.*", "", features_sort[,4]))
#rownames(annotation_col)=rownames(mat_z_corr)
#names(annotation_col) <- "Hctsa_group"

colnames <- features_sort[,4]

col_dist = rev(hcl.colors(100, "inferno"))
#col_label = list(Hctsa_group = c(correlation="#E78AC3", stationarity="#B3E2CD", FourierSpectrum="#BC80BD", information="#CAB2D6",scaling="#BEAED4",model="#66C2A5",
#                                 outliers="#8DA0CB",stepdetection="#FB9A99",entropy="#FDDAEC"))

windowsFonts(A = windowsFont("Calibri"))

Cairo::Cairo(
  45, #length
  25, #width
  file = "G:/JF/HCTSA feature gradients/CamCan/HierarchyCorrelates/PLOTS/Fig2a_cluster_interaction_hierarchy_age.png",
  type = "png", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 500,
  units = "cm" #you can change to pixels etc 
)

x <- pheatmap(mat_z_corr,show_colnames=FALSE,
              scale = "none",clustering_method="ward.D2",
              clustering_distance_cols="euclidean", clustering_distance_rows="euclidean",
              #cutree_rows = 4,cutree_cols = 4,
              #annotation_col = annotation_col, annotation_row = annotation_col, annotation_names_row = FALSE,
              legend = TRUE, legend_labels = "absulute correlation", fontsize_row = 12,
              cellwidth = 12, cellheight = 12,border_color = NA, color = col_dist[1:80],
              labels_row = annotation_col[,1]) #+

dev.off()
