library(Seurat)
library(DIALOGUE)

npcs <- function(seuratobj) {
	pct <- seuratobj[["pca"]]@stdev / sum(seuratobj[["pca"]]@stdev) * 100
	# Calculate cumulative percents for each PC
	cumu <- cumsum(pct)
	# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
	co1 <- which(cumu > 90 & pct < 5)[1]
	# Determine the difference between variation of PC and subsequent PC
	co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
	pcs <- min(co1, co2)
	pcs
}


FullDset <- readRDS("pt-integrated-pal.rds")

DefaultAssay(FullDset) <- "integrated"

pal.list <- SplitObject(FullDset, split.by = "immune_label")

pal.list <- lapply(X = pal.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
  x <- ScaleData(x, features=rownames(x))
  x <- RunPCA(x, features = VariableFeatures(object = x))
})

celltypelist <- list()
for (celltype in names(pal.list)){
  pcs <- npcs(pal.list[[celltype]])
  
  # needs to be appended and not a lapply -- I don't know why, something to do with named variables?
  celltypelist <- append(celltypelist, make.cell.type(name = celltype, tpm = as.matrix(pal.list[[celltype]]@assays$integrated@data), 
                                                      samples = as.matrix(pal.list[[celltype]]$group), 
                                                      X =as.matrix(pal.list[[celltype]]@reductions$pca@cell.embeddings[,1:pcs]),
                                                      metadata = pal.list[[celltype]]@meta.data[,c("orig.ident", "immune_label_fine")],
                                                      cellQ =pal.list[[celltype]]@meta.data[,"nFeature_RNA"] ))
}

names(celltypelist) <- names(pal.list)

R_PCs_mixed <- DIALOGUE.run(rA=celltypelist
                                , main = 'pal', k=10, results.dir = "out_PCs_mixed/", 
                                plot.flag = F, conf = "cellQ")

saveRDS(R_PCs_mixed, file = "out_PCs_mixed/output_object.rds")

