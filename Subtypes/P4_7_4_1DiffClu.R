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
library(plyr)
library(GSVA)
library(ggpubr)
library(gridExtra)
library(DESeq2)
library(crayon)

out.dir <- "P4_7Clustering/P4_7_4DiffClu/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}

out.prefix <- "P4_7_4_1_"

# load data
source("P4_1Data/P4_1_0Colors/P4_1_0colors.R")
source("P4_0Functions/P4_0_1functions.R")
source("P4_7Clustering/P4_7_2_0ClusterPlots.R")

# fpk <- readRDS("P4_2RNA/P4_2_2DEG/P4_2_2_1_tpm.used.rds")
# colnames(fpk) <- mapvalues(colnames(fpk), samInfoOrgSel$Sample.ID.for.analysis, samInfoOrgSel$RepID)
cnt <- readRDS("P4_2RNA/P4_2_2DEG/P4_2_2_1_counts.used.rds")
dim(cnt)
cnt[1:5, 1:5]
cnt.samid <- aggregate(t(cnt), by = list(substr(colnames(cnt), 1, 5)), FUN = mean)
rownames(cnt.samid) <- cnt.samid$Group.1
cnt.samid <- t(cnt.samid[, -1])

myCounts <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/P4_4_4_1_PeakCounts.rds")
colnames(myCounts) <- gsub("^P4_4_2_1_", "", gsub(".final.bam$", "", colnames(myCounts)))
# samples has zero expression
sam <- names(which(colSums(assay(myCounts)) == 0))
usedSam <- setdiff(colnames(myCounts), sam)
usedSam <- usedSam[-grep("N$", usedSam)]
myCounts <- myCounts[, usedSam]
myCounts[1:5, 1:5]

# Get significance
aa <- load("P4_7Clustering/P4_7_2ATACClusters/P4_7_2_1_res.Primary_Metastatic_Clusters.RData")
aa
optK <- 4
results <- rcc3
Cluster <- paste0("Cluster", results[[optK]]$consensusClass)
names(Cluster) <- names(results[[optK]]$consensusClass)

# Diff Peaks --------------------------------------------------------------

Group <- factor(mapvalues(colnames(myCounts), names(Cluster), Cluster))
metaData <- data.frame(Group, row.names = colnames(myCounts))

summary(colSums(assay(myCounts)))
sort(colSums(assay(myCounts)))[1:5]

# Deseq running
# atacDDS <- DESeqDataSetFromMatrix(assay(myCounts), metaData, ~Group, rowRanges = rowRanges(myCounts))
# atacDDS <- DESeq(atacDDS)
# 
# saveRDS(atacDDS, file = paste0(out.dir, out.prefix, "atacDDS.rds"))

for (clu in unique(Cluster)) {
  # Get diff peaks of Cluster1
  # clu <- "Cluster1"
  cat(yellow$bgMagenta$bold(paste0("Diff-running ", clu, " ......\n")))
  
  metaDataN <- metaData
  metaDataN$Group <- mapvalues(metaDataN$Group, setdiff(Cluster, clu), c("Other", "Other", "Other"))
  metaDataN$Group <- factor(metaDataN$Group, levels = c("Other", clu))
  atacDDS <- DESeqDataSetFromMatrix(assay(myCounts), metaDataN, ~Group, rowRanges = rowRanges(myCounts))
  atacDDS <- DESeq(atacDDS)
  
  saveRDS(atacDDS, file = paste0(out.dir, out.prefix, "atacDDS", clu, ".rds"))
}

# DESeq RNA-seq

sam.int <- intersect(names(Cluster), colnames(cnt.samid))
cnt.p <- cnt.samid[, sam.int]
dim(cnt.p)
cnt.p[1:5, 1:5]

lisDEG <- list()
for (clu in unique(Cluster)) {
  # print(clu)
  # clu <- "Cluster1"
  # see vignette for suggestions on generating
  cat(yellow$bgMagenta$bold(paste0("Diff-running ", clu, " ......\n")))
 
  Group <- mapvalues(colnames(cnt.p), names(Cluster), Cluster)
  Group <- mapvalues(Group, setdiff(Cluster, clu), c("Other", "Other", "Other"))
  Group <- factor(Group, levels = c("Other", clu))
  names(Group) <- colnames(cnt.p)
  
  # object construction
  dds <- DESeqDataSetFromMatrix(round(cnt.p), DataFrame(Group), ~ Group)
  
  # standard analysis
  dds <- DESeq(dds)
  res <- results(dds)
  print(head(res))
  saveRDS(res, file = paste0(out.dir, out.prefix, "Res", clu, ".rds"))
  
  pos <- subset(na.omit(res), log2FoldChange > 1 & padj < 0.01)
  print(dim(pos))
  
  lisDEG[[clu]] <-
    rownames(pos)
}

saveRDS(lisDEG, file = paste0(out.dir, out.prefix, "lisDEG.rds"))
saveRDS(cnt.p, file = paste0(out.dir, out.prefix, "cnt.merged.rds"))

fpk <- readRDS("P4_2RNA/P4_2_2DEG/P4_2_2_1_tpm.used.rds")
dim(fpk)
fpk[1:5, 1:5]
fpk.samid <- aggregate(t(fpk), by = list(substr(colnames(fpk), 1, 5)), FUN = mean)
rownames(fpk.samid) <- fpk.samid$Group.1
fpk.samid <- t(fpk.samid[, -1])
setdiff(colnames(cnt.p), colnames(fpk.samid))
fpk.samid <- fpk.samid[, colnames(cnt.p)]

saveRDS(fpk.samid, file = paste0(out.dir, out.prefix, "fpk.merged.rds"))
