# source before this: 
# import_label_TNBC.R
# import_label_ER.R
# import_label_BRCA.R
# import_label_HER2.R
# (haven't done: import_label_pairedER)
library(Seurat)

FullDset <- readRDS("labeled_filtered_pal.rds")

FullDset <- CreateSeuratObject(FullDset@assays$RNA@counts, meta.data = FullDset@meta.data[,c("orig.ident", "immune_label", "immune_label_fine", "nCount_RNA", "nFeature_RNA", "group")])



# split the dataset into a list of two seurat objects (stim and CTRL)
pt.list <- SplitObject(FullDset, split.by = "group")

# normalize and identify variable features for each dataset independently
pt.list <- lapply(X = pt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = pt.list)
pt.anchors <- FindIntegrationAnchors(object.list = pt.list, anchor.features = features)

pt.combined <- IntegrateData(anchorset = pt.anchors)

saveRDS(pt.combined, file = "pt-integrated-pal.rds")

###################################################

# debug -- try doing the dim reduction all together??
FullDset <- FindVariableFeatures(FullDset, selection.method = "vst", nfeatures=1000)
FullDset <- NormalizeData(FullDset)
FullDset <- ScaleData(FullDset, features = VariableFeatures(FullDset))
FullDset <- RunPCA(FullDset, features = VariableFeatures(FullDset))


# Determine percent of variation associated with each PC
pct <- FullDset[["pca"]]@stdev / sum(FullDset[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
pcs <- min(co1, co2)
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

FullDset <- RunUMAP(FullDset, dims = 1:17)


celltypelist <- list()
for (celltype in levels(factor(levels(factor(FullDset$immune_label))))){
  celltypelist <- append(celltypelist, make.cell.type(name = celltype, tpm = as.matrix(FullDset@assays$RNA@data[,FullDset$immune_label == celltype]), 
                                                      samples = as.matrix(FullDset$group[FullDset$immune_label == celltype]), 
                                                      X =as.matrix(FullDset@reductions$pca@cell.embeddings[FullDset$immune_label == celltype,1:20]),
                                                      metadata = FullDset@meta.data[FullDset$immune_label == celltype,c("orig.ident", "immune_label_fine")],
                                                      cellQ =FullDset@meta.data[FullDset$immune_label==celltype,"nFeature_RNA"] ))
}

names(celltypelist) <- levels(factor(FullDset$immune_label))


R_20PCs_matched <- DIALOGUE.run(rA=celltypelist
                                , main = 'pal', k=10, results.dir = "out_20PCs/", 
                                plot.flag = F, conf = "cellQ")


# repeat for 17 PCs
celltypelist <- list()
for (celltype in levels(factor(levels(factor(FullDset$immune_label))))){
  celltypelist <- append(celltypelist, make.cell.type(name = celltype, tpm = as.matrix(FullDset@assays$RNA@data[,FullDset$immune_label == celltype]), 
                                                      samples = as.matrix(FullDset$group[FullDset$immune_label == celltype]), 
                                                      X =as.matrix(FullDset@reductions$pca@cell.embeddings[FullDset$immune_label == celltype,1:17]),
                                                      metadata = FullDset@meta.data[FullDset$immune_label == celltype,c("orig.ident", "immune_label_fine")],
                                                      cellQ =FullDset@meta.data[FullDset$immune_label==celltype,"nFeature_RNA"] ))
}

names(celltypelist) <- levels(factor(FullDset$immune_label))


R_17PCs_matched <- DIALOGUE.run(rA=celltypelist
                                , main = 'pal', k=10, results.dir = "out_17PCs/", 
                                plot.flag = F, conf = "cellQ")



# for unknown reasons, the for loopp works and the lapply doesn't
# so I guess we're just using the for loop from here on out! 

# let's first try splitting by cell type and doing PCA for each
pal.list <- SplitObject(FullDset, split.by = "immune_label")

# normalize and identify variable features for each dataset independently
# skip this during current debug

# this is giving the error below and I'm not sure why -- could be a weird formatting thing?
# rerun using code that looks more like the code from previous runs

pal.list <- lapply(X = pal.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
  x <- NormalizeData(x)
  x <- ScaleData(x, features=rownames(x))
  x <- RunPCA(x, features = VariableFeatures(object = x))
})


celltypelist <- list()
for (celltype in names(pal.list)){
  pct <- pal.list[[celltype]][["pca"]]@stdev / sum(pal.list[[celltype]][["pca"]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  # needs to be appended and not a lapply -- I don't know why, something to do with named variables?
  celltypelist <- append(celltypelist, make.cell.type(name = celltype, tpm = as.matrix(pal.list[[celltype]]@assays$RNA@data), 
                                                      samples = as.matrix(pal.list[[celltype]]$group), 
                                                      X =as.matrix(pal.list[[celltype]]@reductions$pca@cell.embeddings[,1:pcs]),
                                                      metadata = pal.list[[celltype]]@meta.data[,c("orig.ident", "immune_label_fine")],
                                                      cellQ =pal.list[[celltype]]@meta.data[,"nFeature_RNA"] ))
}

names(celltypelist) <- names(pal.list)

R_mixed <- DIALOGUE.run(rA=celltypelist
                                , main = 'pal', k=10, results.dir = "out_mixed/", 
                                plot.flag = F, conf = "cellQ")


saveRDS(celltypelist, file = "celltypelist_PCA_mixed_121721.rds")
  #celltypelist <- lapply( X = pal.list, FUN = function(x) {
#   x <- make.cell.type(name = names(x), tpm = as.matrix(x@assays$RNA@data), 
#                                                                  samples = x$group, 
#                                                                  X =as.matrix(x@reductions$pca@cell.embeddings[,1:18]),
#                                                                  metadata = x@meta.data[,c("orig.ident", "immune_label_fine")],
#                                                                  cellQ =x$nFeature_RNA)
#})

# possible issue: maybe if the expression for a give gene is 0 across a whole cell type?
# "Error in R$gene.pval[[r1@name]]

