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
out.dir <- "P4_8TFNet/P4_8_1GetNet/P4_8_1_4_TF2Peak/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_8_1_4_"

in.dir <- "P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakEach/"
in.prefix <- "P4_8_1_2_"
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

# fil.bed <- paste0(in.dir, out.prefix, "dat.bed.rds")
# dat.bed <- readRDS(fil.bed)

# Get clusters
aa <- load("P4_7Clustering/P4_7_2ATACClusters/P4_7_2_1_res.Primary_Metastatic_Clusters.RData")
aa
optK <- 4
results <- rcc3
Cluster <- paste0("Cluster", results[[optK]]$consensusClass)
names(Cluster) <- names(results[[optK]]$consensusClass)

cnt <- readRDS("P4_5Int/P4_5_2TFnet/P4_5_2_1Peak2Gene/P4_5_2_1_cnt.rds") 

# Get bam cnt -------------------------------------------------------------

fil.cnt <- paste0(in.dir, in.prefix, "cnt.norm.rds")
bamsToCount <- readRDS(paste0(in.dir, in.prefix, "bamsToCount.rds"))
cnt.norm <- readRDS(fil.cnt)

allPeaksSet <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_3_allPeaksSet.rds")
seqAll <- paste0(allPeaksSet@seqnames, allPeaksSet@ranges)

# Read bed and combine peaks ----------------------------------------------

# sams <- read.table("P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakEach/P4_8_1_2_samples.txt")$V1
# sams <- sams[1:5]
s <- "RefSam"

fil.mat <- paste0(out.dir, out.prefix, s, "overlapMatrix.rds")
fil.con <- paste0(out.dir, out.prefix, s, "con.peak.rds")
if (!file.exists(fil.mat) | !file.exists(fil.con)) {
  myPeaks <- readRDS(paste0(in.dir, in.prefix, s, "/", in.prefix, "myPeaks.rds"))
  
  overlap <- list()
  for (i in 1:length(myPeaks)) {
    overlap[[i]] <- allPeaksSet %over% myPeaks[[i]]
  }
  names(overlap) <- names(myPeaks)
  
  overlapMatrix <- do.call(cbind, overlap)
  
  saveRDS(overlapMatrix, file = fil.mat)
  
  
  con.peak <- rowSums(overlapMatrix)
  con.peak <- as.numeric(con.peak > 0)
  names(con.peak) <- seqAll
  
  saveRDS(con.peak, file = fil.con)
} else {
  overlapMatrix <- readRDS(fil.mat)
}

# TF2Peak -----------------------------------------------------------------

rownames(overlapMatrix) <- seqAll

index <- 1
for(i in which(colSums(overlapMatrix) > 0)) {
  print(paste0("TF2Peak: ", i, " of ", ncol(overlapMatrix)))
  dat <-
    data.frame(Peak = rownames(overlapMatrix)[which(overlapMatrix[, i])],
               TF = gsub(".*_", "", colnames(overlapMatrix)[i]),
               SampleID = gsub("_.*", "", colnames(overlapMatrix)[i]))
  print(paste0(nrow(overlapMatrix), " ", nrow(dat)))
  if (index == 1) {
    index <- index + 1
    dat.TF2Peak <- dat
  } else {
    dat.TF2Peak <- rbind(dat, dat.TF2Peak)
  }
}
saveRDS(dat.TF2Peak, file =  paste0(out.dir, out.prefix, s, "TF2Peak.rds"))





