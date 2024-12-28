rm(list = ls())
setwd("/chenfeilab/Avocado/P4_Organoids/")
wd <- "/chenfeilab/Avocado/P4_Organoids/"
library(ggplot2)
library(ggthemes)
library(ggVennDiagram)
library(ComplexHeatmap)
library(Seurat)
library(reshape)
library(ConsensusClusterPlus)
library(GenomicRanges)
library(GenomicAlignments)
library(plyr)

out.dir <- "P4_7Clustering/P4_7_2ATACClusters/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_7_2_1_"

# All Info
samInfoATAC <- readRDS("P4_4ATAC/P4_4_1QC/P4_4_1_1_samInfoATAC.rds")
head(samInfoATAC)
cli <- read.table("P4_1Data/Copy of Table_clinical information_1_2.txt", sep ="\t", header = T, check.names = F)
cli$SampleID <- paste0(cli$`Patient ID`, "T")
head(cli)

samInfoATAC <- join(samInfoATAC, cli, by = "SampleID")

sam.all <- table(samInfoATAC$SampleID)

samMergePeak <- readRDS("P4_4ATAC/P4_4_2Peaks/P4_4_2_1_samMergePeak.rds")
samInfoFilter <- readRDS("P4_4ATAC/P4_4_2Peaks/P4_4_2_1_samInfoFilter.rds")

# colors
source("P4_1Data/P4_1_0Colors/P4_1_0colors.R")
source("/chenfeilab/Avocado/MyPipline/peakVenn.R")
source("P4_0Functions/P4_0_1functions.R")
source("P4_7Clustering/P4_7_2_0ClusterPlots.R")

# peak counts
myCounts <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/P4_4_4_1_PeakCounts.rds")
myCounts[1:5, 1:5]
colnames(myCounts) <- gsub("^P4_4_2_1_", "", gsub(".final.bam$", "", colnames(myCounts)))
# samples has zero expression
sam <- names(which(colSums(assay(myCounts)) == 0))
usedSam <- setdiff(colnames(myCounts), sam)
myCounts <- myCounts[, usedSam]

myPeaks <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_5RXX_peakView/P4_4_4_5RXX_sigPeaks.rds")

# Clustering P&M -----------------------------------------------------------

# Top peaks
topPeaks <- lapply(myPeaks, function(x){
  if (length(x) > 1000) {
    x <- x[order(x$log2FoldChange, decreasing = T)[1:1000], ]
  }
  return(x)
})

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
allPeaksSet_nR <- GenomicRanges::reduce(unlist(GRangesList(topPeaks)))
# allPeaksSet_nR <- reduce(GRangesList(myPeaks))
allPeaksSet_nR

overlap <- list()
for (i in 1:length(topPeaks)) {
  overlap[[i]] <- allPeaksSet_nR %over% topPeaks[[i]]
}
names(overlap) <- names(topPeaks)
topPeakName <-
  lapply(topPeaks, function(x)paste(x@seqnames, x@ranges, sep = "_"))

couPeakName <-
  paste(myCounts@rowRanges@seqnames, myCounts@rowRanges@ranges, sep = "_")

lapply(topPeakName, function(x)length(setdiff(x, couPeakName)))

mat.Peak <-
  assay(myCounts)
rownames(mat.Peak) <- couPeakName
mat.Peak[1:5, 1:5]

mat.Peak <-
  mat.Peak[unique(unlist(topPeakName)), ]
mat.Peak <- log2(mat.Peak + 1)

# saveRDS(mat.Peak, file = paste0(out.dir, out.prefix, "mat.Peak.rds"))

# Get Clusters P M -----------------------------------------------------------

samInfo <- subset(samInfoATAC, Source %in% c("Primary", "Metastatic"))

peak.use <- mat.Peak[, intersect(colnames(mat.Peak), samInfo$SampleID)]

setwd(out.dir)
rcc3 <- ConsensusClusterPlus(
  as.matrix(peak.use),
  maxK = 20,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  title = paste0(out.prefix, "Primary_Metastatic_Clusters"),
  plot = "pdf",
  seed = 123,
  writeTable = T,
  distance = "euclidean",
  clusterAlg = "km"
)
setwd(wd)

save(peak.use, rcc3, samInfo, file = paste0(out.dir, out.prefix, "res.Primary_Metastatic_Clusters.RData"))

# GetPlots(results = rcc3,
#          Kvec = 2:20,
#          intersection_geneExp = peak.use,
#          cli = data.frame(sam = samInfo$SampleID, cli = samInfo$`pN stage`),
#          out.prefix = paste0(out.dir, out.prefix, "PM_"))

# Clustering P -----------------------------------------------------------

samInfo <- subset(samInfoATAC, Source %in% c("Primary"))

peak.use <- mat.Peak[, intersect(colnames(mat.Peak), samInfo$SampleID)]

setwd(out.dir)
rcc4 <- ConsensusClusterPlus(
  as.matrix(peak.use),
  maxK = 20,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  title = paste0(out.prefix, "Primary_Clusters"),
  plot = "pdf",
  seed = 123,
  writeTable = T,
  distance = "euclidean",
  clusterAlg = "km"
)
setwd(wd)

save(peak.use, rcc4, samInfo, file = paste0(out.dir, out.prefix, "res.Primary_Clusters.RData"))

# GetPlots(results = rcc4,
#          Kvec = 2:20,
#          intersection_geneExp = peak.use,
#          cli = data.frame(sam = samInfo$SampleID, cli = samInfo$`pN stage`),
#          out.prefix = paste0(out.dir, out.prefix, "P_"))
