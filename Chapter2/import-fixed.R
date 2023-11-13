# running in bulk-rna environment
# imports from GoogleBucket
# calculates overlaps with Homo_sapiens.GRCh37 as it should
# saves hdf5 


library(Rsamtools)
library(readxl)
#library(GenomeInfoDb)
library(GenomicAlignments)
library(GenomicFeatures)
library(DESeq2)
library(SummarizedExperiment)
library(HDF5Array)
library(stringr)

setwd("~/Aspirin-rerun/")
# create list of all bamfiles in directory
bfl <- BamFileList(dir(".", "bam$", recursive = TRUE), yieldSize = 2e6)

#import xlsx of sample info 

#df[grep("eta", df$Title),]
#BAM names are in sample_metadata$`Collaborator Sample ID`

# starting comparisons
# xxx-V1-xxx-0-x vs xxx-V1-xxx-60-x
# xxx-V1-xxx-0-x vs xxx-V2-xxx-0-x

# also pull only ones that aren't classified as "On Risk" to start with

#tutorial: https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#locating-bam-files-and-the-sample-table
# https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#introduction
si1 <- seqinfo(bfl[1])
#okay, looks like chromosomes, mito, some poorly localized contigs, but essentially as expected

#define gene models

#apply GenomicAlignments
#pulls h sapiens by default
gtffile <- file.path("data","Homo_sapiens.GRCh37.87.gtf.gz")

#ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
txdb <- makeTxDbFromGFF(gtffile)


ebg <- exonsBy(txdb, by="gene")
#rerun to confirm that bfl table order maintained--wondering if rewriting added the wrong transcripts...

se <- summarizeOverlaps(features=ebg, 
                        reads=bfl,
                        mode="Union", #allows some sloppiness to gene models
                        singleEnd=FALSE,#paired-end reads
                        ignore.strand=FALSE, #strand specific
                        fragments=TRUE # count paired reads when only one of the pair aligns
                        )

saveHDF5SummarizedExperiment(se, dir = "se-rerun")
#se2 <- loadHDF5SummarizedExperiment( dir = "se-rerun")

colnames(se) <-  sub(".bam", "", colnames(se))
sample_metadata <- as.data.frame(read_excel("Marks (Harvard) - RNASeq Results - 2019.06.12.xlsx"))
sample_metadata <- DataFrame(sample_metadata)
rownames(sample_metadata) <- sample_metadata$Collaborator.Sample.ID
sample_metadata <- sample_metadata[colnames(se),]
colData(se) <- sample_metadata

se3 <- se


#set 1 metadata
#sample_metadata[grep("-V1-.*-0-", sample_metadata$`Collaborator Sample ID`),]
#sample_metadata[grep("-V1-.*-60-", sample_metadata$`Collaborator Sample ID`),]
#sample_metadata[grep("-V2-.*-0-", sample_metadata$`Collaborator Sample ID`),]


# label set indicates conditions for comparison
colData(se3)$condition <- rep('0', length(se3$Collaborator.Sample.ID))
colData(se3[,grep("-V1-.*-0-", se3$Collaborator.Sample.ID)])$condition <- 'A'
colData(se3[,grep("-V1-.*-60-", se3$Collaborator.Sample.ID)])$condition <- 'B'
colData(se3[,grep("-V2-.*-0-", se3$Collaborator.Sample.ID)])$condition <- 'C'
colData(se3[,grep("-V3-.*-0-", se3$Collaborator.Sample.ID)])$condition <- 'D'
#colData(se[,grep("-V3-.*-0-", sample_metadata$`Collaborator Sample ID`)])$condition <- 'D'

saveHDF5SummarizedExperiment(se3, dir = "se-rerun-full" )


#saveHDF5SummarizedExperiment(se, dir = "aspirin-hdf5" )
#rm(list=ls())
#se <- loadHDF5SummarizedExperiment(dir="aspirin-hdf5")
#se$bfl_id <- names(bfl)
#saveHDF5SummarizedExperiment(se, dir = "se-hdf5-bfl-fix", replace = TRUE )
#saveHDF5SummarizedExperiment(se, dir = "se-hdf5-bfl-fix2", replace = TRUE )
#colData(se) <- colData(se)[order(se$Collaborator.Sample.ID),]
