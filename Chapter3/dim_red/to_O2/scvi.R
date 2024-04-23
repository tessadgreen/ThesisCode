library(Seurat)
library(DIALOGUE)


# pull in filttered, nott this
FullDset <- readRDS("labeled_filtered_pal.rds")

scvi <- read.csv("scvi-output/scvi_coords.csv", row.names = 1)

DefaultAssay(FullDset) <- "integrated"

FullDset <- ScaleData(FullDset, verbose = FALSE)
FullDset <- RunPCA(FullDset, npcs = 50, verbose = FALSE)

FullDset <- AddMetaData(FullDset, scvi)

celltypelist <- list()
for (celltype in levels(factor(levels(factor(FullDset$immune_label))))){
  celltypelist <- append(celltypelist, make.cell.type(name = celltype, tpm = as.matrix(FullDset@assays$RNA@data[,FullDset$immune_label == celltype]), 
                                                      samples = as.matrix(FullDset$group[FullDset$immune_label == celltype]), 
                                                      X =as.matrix(FullDset@meta.data[FullDset$immune_label == celltype, unlist(lapply(0:9, FUN=function(x) {paste0('X',x)}))]),
                                                      metadata = FullDset@meta.data[FullDset$immune_label == celltype,c("orig.ident", "immune_label_fine")],
                                                      cellQ =FullDset@meta.data[FullDset$immune_label==celltype,"nFeature_RNA"] ))
}




names(celltypelist) <- levels(factor(FullDset$immune_label))

R_scvi_matched <- DIALOGUE.run(rA=celltypelist
                                , main = 'pal', k=10, results.dir = paste0("out_scvi_matched/"), 
                                plot.flag = F, conf = "cellQ")

saveRDS(R_scvi_matched, file = paste0("out_scvi_matched/output_object.rds"))