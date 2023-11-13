ERTotal <- readRDS("SeuratObject_ERTotal.rds")
# 0: EpiTumor
# 1: 
# 2: 
# 3: 
# 4: cycling epithelial tumor CyclingEpiTumor
# 5: epithelial normal
# 6: epithelial normal
# 7: 
# 8: 
# ID from fig 6C

ERTotalSub <-readRDS("SeuratObject_ERTotalSub.rds")
# 0: T cells
# 1: TAMs
# 2: CAFs
# 3: Pericytes
# 4: unknown
# 5: endothelial
# 6: TAMs
# 7: B cells
# 8: Myeloid
# 9: CAFs
# 10: Plasma cells
# 11: unknown
# 12: unknown
# ID from fig 7C

ERTotalSub$immune_label <- ERTotalSub$seurat_clusters

ERTotalSub$immune_label <- revalue(ERTotalSub$immune_label, c(   "0" = "Tcells", 
                                                           "1" = "TAMs",
                                                           "2" = "CAFs",
                                                           "3" = "pericytes",
                                                           "4" = "unknown",
                                                           "5" = "endothelial",
                                                           "6" = "TAMs", 
                                                           "7" = "Bcells",
                                                           "8" = "myeloid",
                                                           "9" = "CAFs",
                                                           "10" = "unknown",
                                                           "11" = "unknown",
                                                           "12" = "unknown"
))


ERTotalTC <- readRDS("SeuratObject_ERTotalTC.rds")
ERTotalTC$Tcell_type <- revalue(ERTotalTC$seurat_clusters, c("0" = "CD8effector",
                                                             "1"= "NaiveResting",
                                                             "2" = "Treg",
                                                             "3" = "plasma",
                                                             "4" = "NK",
                                                             "5" = "unknown",
                                                             "6" = "unknown"))

subdf <- as.data.frame(ERTotalSub$immune_label)
subdf[,1] <- as.character(subdf[,1])
fulldf <- as.data.frame(ERTotal@meta.data)

fulldf$immune_label <- as.character(fulldf$seurat_clusters)
fulldf[colnames(ERTotalSub),'immune_label'] <-subdf[,1]

fulldf[fulldf$immune_label == "0", "immune_label"] <- "EpiTumor"
fulldf[fulldf$immune_label == "4", "immune_label"] <- "CyclingEpiTumor"
fulldf[fulldf$immune_label == "5", "immune_label"] <- "EpiNormal"
fulldf[fulldf$immune_label == "6", "immune_label"] <- "EpiNormal"



ERTotal <- AddMetaData(ERTotal, fulldf[,"immune_label", drop=F]) 
ERTotal$immune_label <- as.factor(ERTotal$immune_label)

subdf <- as.data.frame(ERTotalTC$Tcell_type)
subdf[,1] <- as.character(subdf[,1])
fulldf[colnames(ERTotalTC), 'immune_label'] <- subdf[,1]

ERTotal <- AddMetaData(ERTotal, fulldf[,"immune_label", drop=F], col.name = "immune_label_fine")
ERTotal$immune_label_fine <- as.factor(ERTotal$immune_label_fine)

