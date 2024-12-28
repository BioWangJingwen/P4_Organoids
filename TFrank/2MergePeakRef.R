rm(list = ls())

# reference:  Therapy-Induced Evolution of Human Lung Cancer Revealed by Single-Cell RNA Sequencing 
# DOI:https://doi.org/10.1016/j.cell.2020.07.017
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(parallel)
library(doParallel)
# registerDoParallel(cores = 2) 
library(GenomicAlignments)
library(ComplexHeatmap)
library(plyr)
library(circlize)
library(crayon)
# Single-Cell Transcriptomic Analysis of Primary and Metastatic
# Tumor Ecosystems in Head and Neck Cancer,
library(NMF)
library(RColorBrewer)
library(DESeq2)
library(ChIPseeker)
library(IRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# TxDb.Hsapiens.UCSC.hg19.knownGene

# memory.limit(size=xxx)

setwd("/chenfeilab/Avocado/P4_Organoids/")
wd <- "/chenfeilab/Avocado/P4_Organoids/"

out.dir <- "P4_8TFNet/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.dir <- "P4_8TFNet/P4_8_1GetNet/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.dir <- "P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakEach/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_8_1_2_"

# read expression

motif2gene <- 
  read.table("/chenfeilab/Avocado/MyPipline/TOBIAS/test_data/motif2gene_mapping.txt", sep = "\t")
colnames(motif2gene) <- c("Motif", "ENSEMBL")
head(motif2gene)
# This is a file containing information of the origin gene for each TF. 
# By "origin", what is meant is the gene from which the transcription factor is coded. 
# The example-file 'motif2gene_mapping.txt' within the test data gives a good impression of how such a file should look.
library(org.Hs.eg.db)
trs <- clusterProfiler::bitr(motif2gene$ENSEMBL, OrgDb = "org.Hs.eg.db", fromType = "ENSEMBL", toType = "SYMBOL")
head(trs)
motif2tf <- join(motif2gene, trs, by = "ENSEMBL")
motif2tf <- na.omit(motif2tf)
head(motif2tf)

samMergePeak <- readRDS("P4_4ATAC/P4_4_2Peaks/P4_4_2_1_samMergePeak.rds")

fil.bed <- paste0(out.dir, out.prefix, "dat.bed.rds")
if (!file.exists(fil.bed)) {
  dat.bed <- 
    read.table("P4_5Int/P4_5_2TFnet/P4_5_2_1Peak2Gene/P4_5_2_1_dat.bed.txt", header = T, sep = "\t")
  po <- c()
  for (x in 1:length(dat.bed$BED)) {
    if (file.size(dat.bed$BED[x]) != 0) {
      po <- c(po, x)
    }
  }
  dat.bed <- dat.bed[po, ]
  dim(dat.bed)
  head(dat.bed)
} else {
  dat.bed <- readRDS(fil.bed)
}

# Get clusters
aa <- load("P4_7Clustering/P4_7_2ATACClusters/P4_7_2_1_res.Primary_Metastatic_Clusters.RData")
aa
optK <- 4
results <- rcc3
Cluster <- paste0("Cluster", results[[optK]]$consensusClass)
names(Cluster) <- names(results[[optK]]$consensusClass)

cnt <- readRDS("P4_5Int/P4_5_2TFnet/P4_5_2_1Peak2Gene/P4_5_2_1_cnt.rds") 

# Get bam cnt -------------------------------------------------------------

fil.cnt <- paste0(out.dir, out.prefix, "cnt.norm.rds")
if (!file.exists(fil.cnt)) {
  sam <- intersect(colnames(cnt), samMergePeak$SampleID)
  bamsToCount <- gsub("peaks$", "final.bam", samMergePeak$peak[match(sam, samMergePeak$SampleID)])
  cnt.norm <- vst(as.matrix(cnt[, sam]))
  
  saveRDS(bamsToCount, file = paste0(out.dir, out.prefix, "bamsToCount.rds"))
  saveRDS(cnt.norm, file = fil.cnt)
} else {
  bamsToCount <- readRDS(paste0(out.dir, out.prefix, "bamsToCount.rds"))
  cnt.norm <- readRDS(fil.cnt)
}


# Set Samples -------------------------------------------------------------

dat.bed <- subset(dat.bed, SampleID == "RefSam")
out.dir <- "P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakEach/P4_8_1_2_RefSam/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}

# Read bed and combine peaks ----------------------------------------------

print("Read peaks")
fil.p <- paste0(out.dir, out.prefix, "myPeaks.rds")
nam <- paste0(dat.bed$SampleID, "_", dat.bed$TF)
if (!file.exists(fil.p)) {
  # myPeaks <- lapply(dat.bed$BED, ChIPQC:::GetGRanges, simple = TRUE)
  myPeaks <- 
    foreach(i = 1:nrow(dat.bed)) %dopar% {
      if (i %% 1000 == 0) {
        print(paste0("Read peaks: ", i, " of ", nrow(dat.bed)))
      }
      
      ChIPQC:::GetGRanges(dat.bed$BED[i])
    }
  names(myPeaks) <- nam
  saveRDS(myPeaks, file = fil.p)
} else {
  myPeaks <- readRDS(fil.p)
}

fil.p <- paste0(out.dir, out.prefix, "overlap.rds")
if (!file.exists(fil.p)) {
  allPeaksSet_nR <- GenomicRanges::reduce(unlist(GRangesList(myPeaks)))
  # allPeaksSet_nR <- GRangesList(myPeaks)
  allPeaksSet_nR
  
  overlap <- 
    foreach(i = 1:length(myPeaks)) %dopar% {
      if (i %% 1000 == 0) {
        print(paste0("Overlap peaks: ", i, " of ", length(myPeaks)))
      }
      allPeaksSet_nR %over% myPeaks[[i]]
    }
  names(overlap) <- nam
  
  saveRDS(overlap, file = fil.p)
  saveRDS(allPeaksSet_nR, file = paste0(out.dir, out.prefix, "allPeaksSet_nR.rds"))
} else {
  overlap <- readRDS(fil.p)
  allPeaksSet_nR <- readRDS(paste0(out.dir, out.prefix, "allPeaksSet_nR.rds"))
}

