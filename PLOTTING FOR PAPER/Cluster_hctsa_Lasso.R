# library to read matlab data formats into R
library(raveio)
library(tidyverse)
library(pheatmap)
library(dendextend)
library(scico)
library(RColorBrewer)
library(Cairo)

filepath = getwd()


# Load coefficients of lasso regression
coef <- read.csv(paste(filepath,'/DATA/lasso_coef.csv',sep=""), header = F)

# load features and numeric feature matrix and sort them according to this index
mat_z <- read.csv(paste(filepath,"/DATA/lasso_mat.csv" ,sep=""), header = F)

table <- read.csv(paste(filepath,'/DATA/lasso_feat.csv' ,sep=""), header = F)

# reshape mat into a long matrix where participants*parcels x features
df_corr <- abs(cor(mat_z,method = "spearman"))

# clustered heatmap
cols= str_c(table[,2])#,", coef = ", t(round(coef,4)))
rownames(df_corr) <- cols
colnames(df_corr) <- cols

annotation_col = data.frame(gsub(",.*", "", table[,4]))
rownames(annotation_col)=rownames(df_corr)
names(annotation_col) <- "Hctsa_group"

col_dist = rev(hcl.colors(100, "inferno"))
#col_label = list(Hctsa_group = c(correlation="#E78AC3",dynamicalSystem="#CAB2D6",forecasting="#66C2A5",FourierSpectrum="#BC80BD",
#                                 preprocessing="#8DA0CB",scaling="#BEAED4", stationarity="#B3E2CD",stepdetection="#FB9A99",wavelet="#FDDAEC"))

windowsFonts(A = windowsFont("Calibri"))

Cairo::Cairo(
  45, #length
  25, #width
  file = paste(filepath, "/PLOTS/Cluster_bestPerformingFeaturesLasso", ".png", sep = ""),
  type = "tiff", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 1000,
  units = "cm" #you can change to pixels etc 
)

x <- pheatmap(df_corr,show_colnames=FALSE,
              scale = "none",clustering_method="ward.D2",
              clustering_distance_cols="euclidean", clustering_distance_rows="euclidean",
              #cutree_rows = 4,cutree_cols = 4,
              annotation_col = annotation_col, annotation_row = annotation_col, annotation_names_row = FALSE,
              legend = TRUE, legend_labels = "absulute correlation", fontsize_row = 10,
              color = col_dist[1:80]) + # cellwidth = 12, cellheight = 12,border_color = NA, annotation_colors = col_label)+
  theme(text = element_text(family = "A")) 

dev.off()

