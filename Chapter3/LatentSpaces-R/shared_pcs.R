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

FullDset <- ScaleData(FullDset, verbose = FALSE)
FullDset <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)

use.pcs <- npcs(FullDset)

celltypelist <- list()
for (celltype in levels(factor(FullDset$immune_label))){
  celltypelist <- append(celltypelist, make.cell.type(name = celltype, tpm = as.matrix(FullDset@assays$integrated@data[,FullDset$immune_label == celltype]), 
                                                      samples = as.matrix(FullDset$group[FullDset$immune_label == celltype]), 
                                                      X =as.matrix(FullDset@reductions$pca@cell.embeddings[FullDset$immune_label == celltype,1:use.pcs]),
                                                      metadata = FullDset@meta.data[FullDset$immune_label == celltype,c("orig.ident", "immune_label_fine")],
                                                      cellQ =FullDset@meta.data[FullDset$immune_label==celltype,"nFeature_RNA"] ))
}

names(celltypelist) <- levels(factor(FullDset$immune_label))

R_PCs_matched <- DIALOGUE.run(rA=celltypelist
                                , main = 'pal', k=10, results.dir = paste0("out_",use.pcs,"_PCs_matched/"), 
                                plot.flag = F, conf = "cellQ")

saveRDS(R_PCs_matched, file = paste0("out_",use.pcs,"_PCs_matched/output_object.rds"))

