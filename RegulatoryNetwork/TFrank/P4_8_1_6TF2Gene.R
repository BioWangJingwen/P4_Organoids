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
# out.dir <- "P4_8TFNet/P4_8_1GetNet/P4_8_1_6_TF2Gene/"
# if (!dir.exists(out.dir)) {
#   dir.create(out.dir)
# }
out.prefix <- "P4_8_1_6_"

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

fil.cnt <- paste0(in.dir, in.prefix, "cnt.norm.rds")
bamsToCount <- readRDS(paste0(in.dir, in.prefix, "bamsToCount.rds"))
cnt.norm <- readRDS(fil.cnt)

peakBound <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_5_Peak2Gene/P4_8_1_5_peakBound.rds")
seqBound <- paste0(peakBound@seqnames, peakBound@ranges)
peak2gene <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_5_Peak2Gene/P4_8_1_5_peak2gene.rds")
peak2gene <- na.omit(peak2gene)

allPeaksSet <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_3_allPeaksSet.rds")
seqAll <- paste0(allPeaksSet@seqnames, allPeaksSet@ranges)

# Get TF 2 Peak -----------------------------------------------------------

sams <- read.table("P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakEach/P4_8_1_2_samples.txt")$V1
# sams <- sams[1:5]

TF2Peak.lis <- 
  foreach (s = sams) %dopar% {
    print(paste0(s, ": ", match(s, sams), " of ", length(sams)))
    # s = sams[1]
    dat <- readRDS(paste0("P4_8TFNet/P4_8_1GetNet/P4_8_1_4_TF2Peak/P4_8_1_4_", s, "TF2Peak.rds"))
    dat <- subset(dat, Peak %in% seqBound)
  }
names(TF2Peak.lis) <- sams
TF2Peak.dat <- do.call(rbind, TF2Peak.lis)

print("Get TF 2 Gene ... ")

reg.net <- join(TF2Peak.dat, peak2gene, by = "Peak")

saveRDS(reg.net, file =  paste0(out.dir, out.prefix, "reg.net.rds"))




