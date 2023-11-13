library(Seurat)
library(plyr)
setwd("~/Dropbox (HMS)/Tessa_DSM/Datasets/TNBC/Pal/17058077")


TNBC <- readRDS("SeuratObject_TNBC.rds")
# 0: epithelial tumor
# 1: immune
# 2: cycling epithelial tumor
# 3: immune
# 4: immune
# 5: immune
# 6: immune
# 7: immune
# 8: immune
# ID from fig 6a 

TNBCSub <- readRDS("SeuratObject_TNBCSub.rds")
# clusters 1 + 3 from original/full dataset
# non-epithelial TNBC see fig 7A 
# 0 T cells
# 1 TAMs
# 2 Plasma cells
# 3 CAFs
# 4 T cells
# 5 B cells
# 6 DCs
# 7 endothelial
# 8 pericytes
# 9 myeloid
TNBCSub$immune_label <- TNBCSub$seurat_clusters

TNBCSub$immune_label <- revalue(TNBCSub$immune_label, c(   "0" = "Tcells", 
                                   "1" = "TAMs",
                                   "2" = "plasma",
                                   "3" = "CAFs",
                                   "4" = "Tcells",
                                   "5" = "Bcells",
                                   "6" = "DCs", 
                                   "7" = "endothelial",
                                   "8" = "pericytes",
                                   "9" = "myeloid"
))



TNBCTC <- readRDS("SeuratObject_TNBCTC.rds")
TNBCTC$Tcell_type <- revalue(TNBCTC$seurat_clusters, c( "0" = "CD8effector",
                                                        "1" = "NaiveResting",
                                                        "2" = "Treg",
                                                        "3" = "cyclingT",
                                                        "4" = "TRMlike",
                                                        "5" = "NK",
                                                        "6" = "plasma"))
# see Fig EV4A
# TNBC T-cells
# 0 CD8 effector
# 1 naive/resting
# 2 Treg
# 3 cycling T 
# 4 TRM-like
# 5 NK 
# 6 plasma 

TNBCTum <- readRDS("SeuratObject_TNBCTum.rds")
# actually don't need this it's redundant (clusters already split)
subdf <- as.data.frame(TNBCSub$immune_label)
subdf[,1] <- as.character(subdf[,1])
fulldf <- as.data.frame(TNBC@meta.data)

#TNBC <- AddMetaData(TNBC, rep('other', length(colnames(TNBC))), col.name = "immune_label")
fulldf$immune_label <- as.character(fulldf$seurat_clusters)
fulldf[colnames(TNBCSub),'immune_label'] <-subdf[,1]

fulldf[fulldf$immune_label == "0", "immune_label"] <- "EpiTumor"
fulldf[fulldf$immune_label == "2", "immune_label"] <- "CyclingEpiTumor"


TNBC <- AddMetaData(TNBC, fulldf[,"immune_label", drop=F]) 
TNBC$immune_label <- as.factor(TNBC$immune_label)

subdf <- as.data.frame(TNBCTC$Tcell_type)
subdf[,1] <- as.character(subdf[,1])
fulldf[colnames(TNBCTC), 'immune_label'] <- subdf[,1]

TNBC <- AddMetaData(TNBC, fulldf[,"immune_label", drop=F], col.name = "immune_label_fine")
TNBC$immune_label_fine <- as.factor(TNBC$immune_label_fine)



# basal, LP, ML: label and sort
load("HumanBreast10X/Data/Human-PosSigGenes.RData")

# okay so this works, but we don't have any healthy epithelia/relevant tissue here
# so don't worry about it for this exercise.

#z <-TNBC@assays$RNA@data

#z_Basal <- z[rownames(z) %in% Basal, ]
#z_LP <- z[rownames(z) %in% LP, ]
#z_ML <- z[rownames(z) %in% ML, ]
#score_Basal <- colMeans(z_Basal)
#score_LP <- colMeans(z_LP)
#score_ML <- colMeans(z_ML)
#titles <- c("Basal", "LP", "ML")
#dat <- data.frame(score_Basal, score_LP, score_ML)
#TNBC@meta.data <- cbind(TNBC@meta.data, dat)

#DimPlot(TNBC[, colnames(TNBC) %in% colnames(TNBCSub)], label=T)

#HER2Sub <- readRDS("SeuratObject_HER2Sub.rds")
