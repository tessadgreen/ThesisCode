HER2 <- readRDS("SeuratObject_HER2.rds")
# 0: EpiTumor
# 1: 
# 2: 
# 3: cycling epithelial tumor CyclingEpiTumor
# 4: 
# 5:
# 6:
# 7: normal epithelial
# 8: 
# ID from fig 6B

HER2Sub <-readRDS("SeuratObject_HER2Sub.rds")
# 0: TAMs
# 1: T cells
# 2: plasma cells
# 3: CAFs
# 4: endothelial
# 5: B cells
# 6: T cells
# 7: pericytes
# 8: unknwon (ELf3, CDH1, KRT18, FOXA1, ERBB2, CD24)
# 9: unknown
# 10: unknown
# ID from fig 7B

HER2Sub$immune_label <- HER2Sub$seurat_clusters

HER2Sub$immune_label <- revalue(HER2Sub$immune_label, c(   "0" = "TAMs", 
                                                                 "1" = "Tcells",
                                                                 "2" = "plasma",
                                                                 "3" = "CAFs",
                                                                 "4" = "endothelial",
                                                                 "5" = "Bcells",
                                                                 "6" = "Tcells", 
                                                                 "7" = "pericytes",
                                                                 "8" = "unknown",
                                                                 "9" = "unknown",
                                                                 "10" = "unknown"))


HER2TC <- readRDS("SeuratObject_HER2TC.rds")
HER2TC$Tcell_type <- revalue(HER2TC$seurat_clusters, c("0" = "TEM",
                                                             "1"= "Treg",
                                                             "2" = "NaiveResting",
                                                             "3" = "cyclingTrmlike",
                                                             "4" = "NK",
                                                             "5" = "plasma",
                                                             "6" = "Trmlike",
                                                       "7" = "unknown"))

subdf <- as.data.frame(HER2Sub$immune_label)
subdf[,1] <- as.character(subdf[,1])
fulldf <- as.data.frame(HER2@meta.data)

fulldf$immune_label <- as.character(fulldf$seurat_clusters)
fulldf[colnames(HER2Sub),'immune_label'] <-subdf[,1]

fulldf[fulldf$immune_label == "0", "immune_label"] <- "EpiTumor"
fulldf[fulldf$immune_label == "3", "immune_label"] <- "CyclingEpiTumor"
fulldf[fulldf$immune_label == "7", "immune_label"] <- "EpiNormal"



HER2 <- AddMetaData(HER2, fulldf[,"immune_label", drop=F]) 
HER2$immune_label <- as.factor(HER2$immune_label)

subdf <- as.data.frame(HER2TC$Tcell_type)
subdf[,1] <- as.character(subdf[,1])
fulldf[colnames(HER2TC), 'immune_label'] <- subdf[,1]

HER2 <- AddMetaData(HER2, fulldf[,"immune_label", drop=F], col.name = "immune_label_fine")
HER2$immune_label_fine <- as.factor(HER2$immune_label_fine)

