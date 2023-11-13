# source before this: 
# import_label_TNBC.R
# import_label_ER.R
# import_label_BRCA.R
# import_label_HER2.R
# (haven't done: import_label_pairedER)

library(Seurat)
library(DIALOGUE)
setwd("~/Dropbox (HMS)/Tessa_DSM/Datasets/TNBC/Pal/")

info <- read.csv("sample-info.csv")
info$Sample.Name
info$Patient <- as.character(info$Patient)
setwd("~/Dropbox (HMS)/Tessa_DSM/TNBC/pal_given_labels_120621")


FullDset <- merge(BRCA1, c(ERTotal, HER2, TNBC) )

FullDset$immune_label <- revalue(FullDset$immune_label, c("vascEndo" = "endothelial", "lymphEndo" = "endothelial"))
FullDset$immune_label <- revalue(FullDset$immune_label, c("DCs" = "myeloid"))
FullDset$immune_label <- revalue(FullDset$immune_label, c("fibroblasts" = "CAF"))
FullDset$immune_label <- revalue(FullDset$immune_label, c("plasma" = "Bcells"))
FullDset$immune_label <- revalue(FullDset$immune_label, c("CAF" = "CAFs"))

FullDset <- FullDset[ , FullDset$immune_label %in% c("Tcells", "Bcells", "CAFs", "CyclingEpiTumor", 
                                         "EpiTumor", "pericytes", "TAMs", "endothelial") ]

# then check sample membership of each
# 
isecs_pt <- table(FullDset$immune_label, FullDset$group)

FullDset <- FullDset[, FullDset$group != "ER_0001"]

# looks like at least a few cells per cell type in each sample
# better get going! 
saveRDS(FullDset, "labeled_filtered_pal.rds")


