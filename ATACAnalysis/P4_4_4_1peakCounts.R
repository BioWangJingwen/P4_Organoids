rm(list = ls())

# library(ATACseqQC)
library(ChIPseeker)
library(ChIPpeakAnno)
library(idr)
library(reticulate)
library(DiffBind)
library(tidyverse)
library(purrr)
library(parallel)
library(GenomicRanges)
library(GenomicAlignments)

# use_condaenv(condaenv = "/chenfeilab/Avocado/Softwares/miniconda3/", required = NULL)

setwd("/chenfeilab/Avocado/P4_Organoids")
wd <- "/chenfeilab/Avocado/P4_Organoids/"

out.dir <- "P4_4ATAC/P4_4_4peakDiff/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.dir <- "P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_4_4_1_"

# All Info
samInfoATAC <- readRDS("P4_4ATAC/P4_4_1QC/P4_4_1_1_samInfoATAC.rds")
head(samInfoATAC)

sam.all <- table(samInfoATAC$SampleID)

samMergePeak <- readRDS("P4_4ATAC/P4_4_2Peaks/P4_4_2_1_samMergePeak.rds")
samInfoFilter <- readRDS("P4_4ATAC/P4_4_2Peaks/P4_4_2_1_samInfoFilter.rds")

# colors
source("P4_1Data/P4_1_0Colors/P4_1_0colors.R")

# Filter peaks ----------------------------------------------------------------

# get duplication of peaks in samples
peaks <- samMergePeak$peak

myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
# names(myPeaks) <- samMergePeak$SampleID
allPeaksSet_nR <- reduce(unlist(GRangesList(myPeaks)))
# allPeaksSet_nR <- GRangesList(myPeaks)
allPeaksSet_nR

overlap <- list()
for (i in 1:length(myPeaks)) {
  overlap[[i]] <- allPeaksSet_nR %over% myPeaks[[i]]
}
names(overlap) <- samMergePeak$SampleID

overlapMatrix <- do.call(cbind, overlap)
mcols(allPeaksSet_nR) <- overlapMatrix

allPeaksSet_nR[1:2, ]
dim(allPeaksSet_nR)

# delete blacklist and ChrM
blklist <- rtracklayer::import.bed("P4_1Data/hg19-blacklist.v2.bed")
nrToCount <- allPeaksSet_nR[!allPeaksSet_nR %over% blklist & !seqnames(allPeaksSet_nR) %in% "chrM"]
nrToCount

# Reproducible peaks found in at least two of the samples
# library(Rsubread)
occurrences <- rowSums(as.data.frame(elementMetadata(nrToCount)))

nrToCount <- nrToCount[occurrences >= 2, ]
nrToCount <- subset(nrToCount, seqnames %in% paste0("chr", 1:22))
nrToCount[1:5, 1:5]
length(nrToCount@seqnames)

summary(colSums(as.matrix(nrToCount@elementMetadata)))

saveRDS(nrToCount, file = paste0(out.dir, out.prefix, "UsedPeaks.rds"))

# get read counts
# consensusToCount <- GPos(nrToCount)

bamsToCount <- gsub("peaks$", "final.bam", samMergePeak$peak)

# save.image(paste0(out.dir, out.prefix, "image.rds"))
# 
# load(paste0(out.dir, out.prefix, "image.rds"))

myCounts <- summarizeOverlaps(nrToCount, bamsToCount, singleEnd = FALSE)

saveRDS(myCounts, file = paste0(out.dir, out.prefix, "PeakCounts.rds"))

myCounts[1:5, 1:5]

# Diff --------------------------------------------------------------------

library(DESeq2)
colnames(myCounts) <- gsub("^P4_4_2_1_", "", gsub(".final.bam$", "", colnames(myCounts)))
# samples has zero expression
sam <- names(which(colSums(assay(myCounts)) == 0))
usedSam <- setdiff(colnames(myCounts), sam)

myCounts <- myCounts[, usedSam]

Group <- factor(mapvalues(colnames(myCounts), samInfoATAC$SampleID, samInfoATAC$Source), levels = c("Normal", "Primary", "Metastatic"))
metaData <- data.frame(Group, row.names = colnames(myCounts))

summary(colSums(assay(myCounts)))
sort(colSums(assay(myCounts)))[1:5]

# Deseq running
atacDDS <- DESeqDataSetFromMatrix(assay(myCounts), metaData, ~Group, rowRanges = rowRanges(myCounts))
atacDDS <- DESeq(atacDDS)

saveRDS(atacDDS, file = paste0(out.dir, out.prefix, "atacDDS.rds"))

# Get diff peaks of normal
metaDataN <- metaData
metaDataN$Group <- mapvalues(metaDataN$Group, c("Primary", "Metastatic"), c("Tumor", "Tumor"))
metaDataN$Group <- factor(metaDataN$Group, levels = c("Tumor", "Normal"))
atacDDSNormal <- DESeqDataSetFromMatrix(assay(myCounts), metaDataN, ~Group, rowRanges = rowRanges(myCounts))
atacDDSNormal <- DESeq(atacDDSNormal)

saveRDS(atacDDSNormal, file = paste0(out.dir, out.prefix, "atacDDSNormal.rds"))

resNormal <- results(atacDDSNormal, format = "GRanges")
sigNormal <- subset(resNormal, log2FoldChange > 1 & padj < 0.01)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tracktables)
toOverLap <-  promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, 500, 500)
sigBedNormal <- 
  sigNormal[(!is.na(sigNormal$padj) &
               sigNormal$padj < 0.05 & 
               sigNormal$log2FoldChange > 1 &
               sigNormal %over% toOverLap), ]
myReport <- makebedtable(sigBedNormal, paste0(out.prefix, "sigBedNormal.html"), out.dir)
browseURL(myReport)

# Get diff peaks of primary to normal
metaDataN <- subset(metaData, Group %in% c("Normal", "Primary"))
metaDataN$Group <- factor(metaDataN$Group, levels = c("Normal", "Primary"))
myCountsN <- myCounts[, rownames(metaDataN)]

atacDDSP2N <- DESeqDataSetFromMatrix(assay(myCountsN), metaDataN, ~Group, rowRanges = rowRanges(myCountsN))
atacDDSP2N <- DESeq(atacDDSP2N)

saveRDS(atacDDSP2N, file = paste0(out.dir, out.prefix, "atacDDSP2N.rds"))

# Get diff peaks of primary to Metastatic
metaDataN <- subset(metaData, Group %in% c("Metastatic", "Primary"))
metaDataN$Group <- factor(metaDataN$Group, levels = c("Metastatic", "Primary"))
myCountsN <- myCounts[, rownames(metaDataN)]

atacDDSP2M <- DESeqDataSetFromMatrix(assay(myCountsN), metaDataN, ~Group, rowRanges = rowRanges(myCountsN))
atacDDSP2M <- DESeq(atacDDSP2M)

saveRDS(atacDDSP2M, file = paste0(out.dir, out.prefix, "atacDDSP2M.rds"))

# Get diff peaks of Metastatic to normal
metaDataN <- subset(metaData, Group %in% c("Normal", "Metastatic"))
metaDataN$Group <- factor(metaDataN$Group, levels = c("Normal", "Metastatic"))
myCountsN <- myCounts[, rownames(metaDataN)]

atacDDSM2N <- DESeqDataSetFromMatrix(assay(myCountsN), metaDataN, ~Group, rowRanges = rowRanges(myCountsN))
atacDDSM2N <- DESeq(atacDDSM2N)

saveRDS(atacDDSM2N, file = paste0(out.dir, out.prefix, "atacDDSM2N.rds"))

# Get diff peaks of Metastatic to Primary
metaDataN <- subset(metaData, Group %in% c("Primary", "Metastatic"))
metaDataN$Group <- factor(metaDataN$Group, levels = c("Primary", "Metastatic"))
myCountsN <- myCounts[, rownames(metaDataN)]

atacDDSM2P <- DESeqDataSetFromMatrix(assay(myCountsN), metaDataN, ~Group, rowRanges = rowRanges(myCountsN))
atacDDSM2P <- DESeq(atacDDSM2P)

saveRDS(atacDDSM2P, file = paste0(out.dir, out.prefix, "atacDDSM2P.rds"))

