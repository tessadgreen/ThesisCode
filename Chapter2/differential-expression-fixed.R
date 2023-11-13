library(SummarizedExperiment)
library(HDF5Array)
library(DESeq2)
library(biomaRt)
library(apeglm)
library(edgeR)
library(bigPint)
library(Glimma)
library(network)

setwd("~/Aspirin-rerun")
se <- loadHDF5SummarizedExperiment(dir="results-fixed/se-rerun-full")

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host= "grch37.ensembl.org")
name_table_all <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = rownames(se), bmHeader = T, mart = mart)

GetNames <- function(res, name_table=name_table_all) {
  res$symbol <- name_table[match(rownames(res), name_table[,1]),2]
  res[which(res$symbol == ''),]$symbol <- rownames(res[which(res$symbol == ''),])
  res
}
GetNamesBM <-function(res) {
  getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = rownames(res), bmHeader = T, mart = mart)
}

get.gene.name <- function(geneid, name_table = name_table_all){
  if (geneid %in% name_table_all$`Gene stable ID`){
    id <- which(name_table_all$`Gene stable ID`==geneid)
    return(name_table[id,2])
  }
  else{
    return(geneid)
  }
}

#convert counts to matrix in local memory
se@assays@data$counts<- matrix(se@assays@data$counts, nrow =nrow(se@assays@data$counts) )

# format conditions as factors
se$Collaborator.Participant.ID <- substr(se$Collaborator.Participant.ID,6,8)
se$condition <- factor(se$condition)
# A pre, B 60 mins, C visit 2
se$Collaborator.Participant.ID <- factor(se$Collaborator.Participant.ID)

ddsSE <- DESeqDataSet(se, ~condition + Collaborator.Participant.ID)

ddsSE <- ddsSE[rowSums(counts(ddsSE)) >= 50] #filter
ddsSE <-DESeq(ddsSE)
htmlwidgets::saveWidget(glimmaMDS(ddsSE), "./results-fixed/mds-plot-allsamples.html")

res_AvB <- results(ddsSE, contrast = c("condition", "A", "B"))
res_AvB <- GetNames(res_AvB)
res_AvB <- res_AvB[order(res_AvB$padj),]
res_AvB <- res_AvB[complete.cases(res_AvB),]
res_AvB <- results(ddsSE, contrast = c("condition", "A", "B"))

#res_AvB_cooks <- results(ddsSE, contrast = c("condition", "A", "B"))


write.csv(as.data.frame(res_AvB), file= './results-fixed/se_0_v_60_matchedAll_deSEQ.csv')
pdf("./results-fixed/se_0_v_60_matchedAll_deSEQ_plotMA.pdf")
DESeq2::plotMA(res_AvB)
dev.off()
sink("./results-fixed/se_0_v_60_matchedAll_deSEQ_resAvB_summary.txt")
summary(res_AvB)
sink()


res_AvC <- results(ddsSE, contrast = c("condition", "A", "C"))
res_AvC <- GetNames(res_AvC)
res_AvC <- res_AvC[complete.cases(res_AvC),]

#res_AvC_cooks <- results(ddsSE, contrast = c("condition", "A", "C"), cooksCutoff = TRUE)
#res_AvC$cooksCutoff <- is.na(res_AvC_cooks$padj)
res_AvC <- res_AvC[order(res_AvC$padj),]



write.csv(as.data.frame(res_AvC), file= './results-fixed/se_0_v_V2_matchedAll_deSEQ.csv')
summary(res_AvC)
pdf("./results-fixed/se_0_v_V2_matchedAll_deSEQ_plotMA.pdf")
DESeq2::plotMA(res_AvC)
dev.off()
sink('./results-fixed/se_0_v_V2_matchedAll_deSEQ_resAvC_summary.txt')
summary(res_AvC)
sink()

res_AvD <- results(ddsSE, contrast = c("condition", "A", "D"))
res_AvD <- GetNames(res_AvD)
res_AvD <- res_AvD[order(res_AvD$padj),]
res_AvD <- res_AvD[complete.cases(res_AvD),]
write.csv(as.data.frame(res_AvD), file= './results-fixed/se_0_v_v3_matchedAll_deSEQ.csv')
pdf("./results-fixed/se_0_v_V3_matchedAll_deSEQ_plotMA.pdf")
DESeq2::plotMA(res_AvD)
dev.off()
sink('./results-fixed/se_0_v_v3_matchedAll_deSEQ_resAvD_summary.txt')
summary(res_AvD)
sink()

anno <- as.data.frame(rownames(ddsSE))
anno$genename <- lapply(anno$`rownames(ddsSE)`, get.gene.name)
colnames(anno)[1] <- "GeneID"

glMDPlot(ddsSE, groups = ddsSE$condition, counts= ddsSE@assays@data$counts,
         anno=anno,samples=ddsSE$Collaborator.Sample.ID,
         sample.cols = as.color(ddsSE$Collaborator.Participant.ID),
         folder = "results-fixed/MDplot_allsamples", transform=T)

#deseq between time points with only >50M reads
#unmatched -- not aware of patient ID
#ddsSE_50M <- DESeqDataSet(se[, se$Sample.Is.On.Risk == "F"], ~condition)
#ddsSE_50M <- ddsSE_50M[rowSums(counts(ddsSE_50M)) >= 10] #filter
#ddsSE_50M <-DESeq(ddsSE_50M)#

#res_AvB_50M <- results(ddsSE_50M, contrast = c("condition", "A", "B"))
#res_AvC_50M <- results(ddsSE_50M, contrast = c("condition", "A", "C"))
#res_AvB_50M <- res_AvB_50M[order(res_AvB_50M$padj),]
#res_AvC_50M <- res_AvC_50M[order(res_AvC_50M$padj),]

#res_AvB_50M <- GetNames(res_AvB_50M)
#res_AvC_50M <- GetNames(res_AvC_50M)


#plotMA(res_AvB_50M)
#summary(res_AvB_50M)
#plotMA(res_AvC_50M)
#summary(res_AvC_50M)

#write.csv(as.data.frame(res_AvB), file= './results-fixed/se_0_v_60_50M_deSEQ.csv')

#write.csv(as.data.frame(res_AvC), file= './results-fixed/se_0_v_V2_50M_deSEQ.csv')


#SavePlots <- function(dds, res, rows) {
#  for (i in rows){
#    plotCounts(dds, rownames(res)[i], main=res$symbol[i])
#  }
#}

# OK so that's with all samples
# but probably I should check it with just the samples that 'worked'
#se_good <- se[, se$Sample.Is.On.Risk == "F"]
#se_good <- se[,se$Billed.=="Yes"]
se_good <- se[,se$Reads.Aligned.in.Pairs > 25000000]
se_good_AvB <- se_good[, se_good$condition %in% c("A","B")]
se_good_AvC <- se_good[, se_good$condition %in% c("A","C")]
se_good_AvD <- se_good[, se_good$condition %in% c("A","D")]

#se_good_AvB <- se_good_AvB[,duplicated(se_good_AvB$Collaborator.Participant.ID)]

keep_pts <- names(table(se_good_AvB$Collaborator.Participant.ID)[table(se_good_AvB$Collaborator.Participant.ID) == 2])
ddsSE_50M_AvB <- DESeqDataSet(se_good_AvB[,se_good_AvB$Collaborator.Participant.ID %in% keep_pts], ~condition + Collaborator.Participant.ID)
ddsSE_50M_AvB <- ddsSE_50M_AvB[rowSums(counts(ddsSE_50M_AvB)) >= 10] #filter
ddsSE_50M_AvB <-DESeq(ddsSE_50M_AvB)

res_AvB_50M <- results(ddsSE_50M_AvB, contrast = c("condition", "A", "B"))

keep_pts <- names(table(se_good_AvC$Collaborator.Participant.ID)[table(se_good_AvC$Collaborator.Participant.ID) == 2])
ddsSE_50M_AvC <- DESeqDataSet(se_good_AvC[,se_good_AvC$Collaborator.Participant.ID %in% keep_pts], ~condition + Collaborator.Participant.ID)
ddsSE_50M_AvC <- ddsSE_50M_AvC[rowSums(counts(ddsSE_50M_AvC)) >= 10] #filter
ddsSE_50M_AvC <-DESeq(ddsSE_50M_AvC)

res_AvC_50M <- results(ddsSE_50M_AvC, contrast = c("condition", "A", "C"))


keep_pts <- names(table(se_good_AvD$Collaborator.Participant.ID)[table(se_good_AvD$Collaborator.Participant.ID) == 2])
ddsSE_50M_AvD <- DESeqDataSet(se_good_AvD[,se_good_AvD$Collaborator.Participant.ID %in% keep_pts], ~condition + Collaborator.Participant.ID)
ddsSE_50M_AvD <- ddsSE_50M_AvD[rowSums(counts(ddsSE_50M_AvD)) >= 10] #filter
ddsSE_50M_AvD <-DESeq(ddsSE_50M_AvD)

res_AvD_50M <- results(ddsSE_50M_AvD, contrast = c("condition", "A", "D"))




res_AvB_50M <- res_AvB_50M[order(res_AvB_50M$padj),]
res_AvC_50M <- res_AvC_50M[order(res_AvC_50M$padj),]
res_AvD_50M <- res_AvD_50M[order(res_AvD_50M$padj),]

res_AvB_50M <- GetNames(res_AvB_50M)
res_AvC_50M <- GetNames(res_AvC_50M)
res_AvD_50M <- GetNames(res_AvD_50M)


#plotMA(res_AvB_50M)
sink("./results-fixed/ma-plot-50M_AvB_summary.txt")
summary(res_AvB_50M)
sink()
#plotMA(res_AvC_50M)
sink("./results-fixed/ma-plot-50M_AvC_summary.txt")
summary(res_AvC_50M)
sink()
#plotMA(res_AvD_50M)
sink("./results-fixed/ma-plot-50M_AvD_summary.txt")
summary(res_AvD_50M)
sink()

pdf("./results-fixed/ma-plot-50M_AvB.pdf")
DESeq2::plotMA(res_AvB_50M)
dev.off()


pdf("./results-fixed/ma-plot-50M_AvC.pdf")
DESeq2::plotMA(res_AvC_50M)
dev.off()

pdf("./results-fixed/ma-plot-50M_AvD.pdf")
DESeq2::plotMA(res_AvD_50M)
dev.off()


write.csv(as.data.frame(res_AvB_50M), file= './results-fixed/se_0_v_60_50M_deSEQ.csv')

write.csv(as.data.frame(res_AvC_50M), file= './results-fixed/se_0_v_V2_50M_deSEQ.csv')

write.csv(as.data.frame(res_AvD_50M), file= './results-fixed/se_0_v_V3_50M_deSEQ.csv')


anno <- as.data.frame(rownames(ddsSE_50M_AvB))
anno$genename <- lapply(anno$`rownames(ddsSE_50M_AvB)`, get.gene.name)
colnames(anno)[1] <- "GeneID"

glMDPlot(ddsSE_50M_AvB, groups = ddsSE_50M_AvB$condition, counts= ddsSE_50M_AvB@assays@data$counts,
         anno=anno,samples=ddsSE_50M_AvB$Collaborator.Sample.ID,
         sample.cols = as.color(ddsSE_50M_AvB$Collaborator.Participant.ID),
         folder = "results-fixed/MDplot_50M_AvB", transform=T)

anno <- as.data.frame(rownames(ddsSE_50M_AvC))
anno$genename <- lapply(anno$`rownames(ddsSE_50M_AvC)`, get.gene.name)
colnames(anno)[1] <- "GeneID"

glMDPlot(ddsSE_50M_AvC, groups = ddsSE_50M_AvC$condition, counts= ddsSE_50M_AvC@assays@data$counts,
         anno=anno,samples=ddsSE_50M_AvC$Collaborator.Sample.ID,
         sample.cols = as.color(ddsSE_50M_AvC$Collaborator.Participant.ID),
         folder = "results-fixed/MDplot_50M_AvC", transform=T)

anno <- as.data.frame(rownames(ddsSE_50M_AvD))
anno$genename <- lapply(anno$`rownames(ddsSE_50M_AvD)`, get.gene.name)
colnames(anno)[1] <- "GeneID"

glMDPlot(ddsSE_50M_AvD, groups = ddsSE_50M_AvD$condition, counts= ddsSE_50M_AvD@assays@data$counts,
         anno=anno,samples=ddsSE_50M_AvD$Collaborator.Sample.ID,
         sample.cols = as.color(ddsSE_50M_AvD$Collaborator.Participant.ID),
         folder = "results-fixed/MDplot_50M_AvD", transform=T)


#### CLICKABLE WITH JUST THE GOOD SAMPLES
se_good <- se[,se$Reads.Aligned.in.Pairs > 25000000]
# note ENSG00000107165 did not converge TYRP1
ddsSE_good <- DESeqDataSet(se_good, ~condition + Collaborator.Participant.ID)
ddsSE_good <- ddsSE_good[rowSums(counts(ddsSE_good)) >= 50] #filter
ddsSE_good <- estimateSizeFactors(ddsSE_good)
ddsSE_good <- estimateDispersions(ddsSE_good)
ddsSE_good <- nbinomWaldTest(ddsSE_good, maxit=1000)

#ddsSE_good <-DESeq(ddsSE_good)

res_AvB_25M <- results(ddsSE_good, contrast = c("condition", "A", "B"))

anno <- as.data.frame(rownames(ddsSE))
anno$genename <- lapply(anno$`rownames(ddsSE)`, get.gene.name)
colnames(anno)[1] <- "GeneID"



glMDPlot(ddsSE, groups = ddsSE$condition, counts= ddsSE@assays@data$counts,
         anno=anno,samples=ddsSE$Collaborator.Sample.ID,
         sample.cols = as.color(ddsSE$Collaborator.Participant.ID),
         folder = "results-fixed/MDplot_allsamples", transform=T)





# Then a third one: downsample the high-count ones to match the low-count ones
dsworkmat <- downsampleMatrix(workmat, min(colSums(workmat))/colSums(workmat), bycol = TRUE, sink = NULL)
colSums(dsworkmat)



se_subsampled <- se
se_subsampled@assays@data$counts<- as.matrix(dsworkmat)
saveHDF5SummarizedExperiment(se_subsampled, "se-subsampled")

ddsSE <- DESeqDataSet(se_subsampled, ~condition + Collaborator.Participant.ID)
ddsSE <- ddsSE[rowSums(counts(ddsSE)) >= 50] #filter
ddsSE <-DESeq(ddsSE)
htmlwidgets::saveWidget(glimmaMA(ddsSE), "./results-fixed/ma-plot-subsampled.html")

res_AvB <- results(ddsSE, contrast = c("condition", "A", "B"))
res_AvB <- GetNames(res_AvB)
res_AvB <- res_AvB[order(res_AvB$padj),]
write.csv(as.data.frame(res_AvB), file= './results-fixed/se_0_v_60_downsampled_deSEQ.csv')
sink("./results-fixed/downsampled_resAvB.txt")
summary(res_AvB)
sink()

pdf("./results-fixed/se_0_v_60_downsampled_deSEQ_plotMA.pdf")
DESeq2::plotMA(res_AvB)
dev.off()

res_AvC <- results(ddsSE, contrast = c("condition", "A", "C"))
res_AvC <- GetNames(res_AvC)
res_AvC <- res_AvC[order(res_AvC$padj),]

pdf("./results-fixed/se_0_v_V2_downsampled_deSEQ_plotMA.pdf")
DESeq2::plotMA(res_AvC)
dev.off()

write.csv(as.data.frame(res_AvC), file= './results-fixed/se2_0_v_V2_downsampled_deSEQ.csv')
plotMA(res_AvC)
sink("./results-fixed/downsampled_resAvC.txt")
summary(res_AvC)
sink()

res_AvD <- results(ddsSE, contrast = c("condition", "A", "D"))
res_AvD <- GetNames(res_AvD)
res_AvD <- res_AvD[order(res_AvD$padj),]
write.csv(as.data.frame(res_AvD), file= './results-fixed/se_0_v_v3_downsampled_deSEQ.csv')
pdf("./results-fixed/se_0_v_V3_downsampled_deSEQ_plotMA.pdf")
DESeq2::plotMA(res_AvD)
dev.off()

write.csv(as.data.frame(res_AvD), file= './results-fixed/se2_0_v_V3_downsampled_deSEQ.csv')
sink("./results-fixed/downsampled_resAvD.txt")
summary(res_AvD)
sink()

anno <- as.data.frame(rownames(ddsSE))
anno$genename <- lapply(anno$`rownames(ddsSE)`, get.gene.name)
colnames(anno)[1] <- "GeneID"

glMDPlot(ddsSE, groups = ddsSE$condition, counts= ddsSE@assays@data$counts,
         anno=anno,samples=ddsSE$Collaborator.Sample.ID,
         sample.cols = as.color(ddsSE$Collaborator.Participant.ID),
         folder = "results-fixed/MDplot_downsampled", transform=T)

