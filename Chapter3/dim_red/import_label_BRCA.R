BRCA1 <- readRDS("SeuratObject_BRCA1Tum.rds")
# only labeled in subsets
# see 4E 
# 0: cancer
# 5: LP (normal epi)
# 9: ML (normal epi)
# 6: basal
# 8: basal

BRCASub <-readRDS("SeuratObject_BRCA1TumSub.rds")
# 0: fibroblasts
# 1: T cells
# 2: vasc. endothelial
# 3: TAMs
# 4: plasma
# 5: pericytes
# 6: B cells
# 7: Myeloid
# 8: lymph endothelial

# ID from fig 5A

BRCASub$immune_label <- BRCASub$seurat_clusters

BRCASub$immune_label <- revalue(BRCASub$immune_label, c(   "0" = "fibroblasts", 
                                                                 "1" = "Tcells",
                                                                 "2" = "vascEndo",
                                                                 "3" = "TAMs",
                                                                 "4" = "plasma",
                                                                 "5" = "pericytes",
                                                                 "6" = "Bcells", 
                                                                 "7" = "myeloid",
                                                                 "8" = "lymphEndo"
))


# no T cell sub-labeling -- will need to DIY if I want to do this 

subdf <- as.data.frame(BRCASub$immune_label)
subdf[,1] <- as.character(subdf[,1])
fulldf <- as.data.frame(BRCA1@meta.data)

fulldf$immune_label <- as.character(fulldf$seurat_clusters)
fulldf[colnames(BRCASub),'immune_label'] <-subdf[,1]

# note: 0 also contains "CyclingEpiTumor" 
fulldf[fulldf$immune_label == "0", "immune_label"] <- "EpiTumor"
fulldf[fulldf$immune_label == "5", "immune_label"] <- "LP"
fulldf[fulldf$immune_label == "9", "immune_label"] <- "ML"
fulldf[fulldf$immune_label == "6", "immune_label"] <- "basal"
fulldf[fulldf$immune_label == "8", "immune_label"] <- "basal"



BRCA1 <- AddMetaData(BRCA1, fulldf[,"immune_label", drop=F]) 
BRCA1$immune_label <- as.factor(BRCA1$immune_label)

# for T cell type -- I will have to do the clustering myself
# same for cycling epithelial 

# keep only the TNBC samples
BRCA1 <- BRCA1[,startsWith(BRCA1$group, "T")]

# cluster the epithelial to extract cycling/non-cycling
epi1 <- BRCA1[,BRCA1$immune_label == "EpiTumor"]
#epi1 <- FindNeighbors(epi1, dims = 1:30)
#epi1 <- FindClusters(epi1, resolution= 0.1)

# complete separation by patient here 
# need anchors or similar to drive merger 
epi1 <- CreateSeuratObject(epi1@assays$RNA@counts, project = "epi", assay = "RNA", meta.data = epi1@meta.data)


epi.list <- SplitObject(epi1, split.by = "group")

# normalize and identify variable features for each dataset independently
epi.list <- lapply(X = epi.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = epi.list)
epi.anchors <- FindIntegrationAnchors(object.list = epi.list, anchor.features = features)
epi.combined <- IntegrateData(anchorset = epi.anchors)
DefaultAssay(epi.combined) <- "integrated"

epi.combined <- ScaleData(epi.combined)
epi.combined <- RunPCA(epi.combined, npcs = 20)
epi.combined <- RunUMAP(epi.combined, reduction = "pca", dims = 1:20)
epi.combined <- FindNeighbors(epi.combined, reduction = "pca", dims = 1:20)
epi.combined <- FindClusters(epi.combined, resolution = 0.1)

# 

epi.combined$cycling <- revalue(epi.combined$integrated_snn_res.0.1, c(   "0" = "EpiTumor", 
                                                           "1" = "CyclingEpiTumor",
                                                           "2" = "EpiTumor",
                                                           "3" = "EpiTumor"
))


#epi1 <- NormalizeData(epi1)
#epi1 <- FindVariableFeatures(epi1, selection.method = 'vst', nfeatures=1000)
#epi1 <- ScaleData(epi1, features=rownames(epi1))
#epi1 <- RunPCA(epi1, features = VariableFeatures(object = epi1))
#epi1 <- FindNeighbors(epi1, dims = 1:20)
#epi1 <- FindClusters(epi1, resolution = 0.1)
#epi1 <- RunUMAP(epi1, dims = 1:20)
# update the clusters are poorly aligned with MKI67 expression
# maybe compute PAM50 score and see what pops out? 
# variable genes as PAM50 or just make sure MKI67 included??
# 4 is cycling and 1 is not, but I don't know what makes 1/4 distinct from 0/3/2


subdf <- as.data.frame(epi.combined$cycling)
subdf[,1] <- as.character(subdf[,1])
fulldf <- as.data.frame(BRCA1@meta.data)
fulldf$immune_label <- as.character(fulldf$immune_label)
fulldf[colnames(epi.combined),'immune_label'] <-subdf[,1]
fulldf$immune_label <- as.factor(fulldf$immune_label)

BRCA1$immune_label <- fulldf$immune_label

