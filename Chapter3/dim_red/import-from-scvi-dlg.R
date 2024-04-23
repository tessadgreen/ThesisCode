library(Seurat)
library(DIALOGUE)


scvi <- read.csv("scvi-output/scvi_coords.csv", row.names = 1)

FullDset <- AddMetaData(FullDset, scvi)

lapply(0:9, FUN=function(x) {paste0('X',x)})