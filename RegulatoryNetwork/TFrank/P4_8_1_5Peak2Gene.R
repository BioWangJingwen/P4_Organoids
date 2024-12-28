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
out.dir <- "P4_8TFNet/P4_8_1GetNet/P4_8_1_5_Peak2Gene/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_8_1_5_"

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

fil.p <- paste0(out.dir, out.prefix, "peakBound.rds")
if (!file.exists(fil.p)) {
  allPeaksSet <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_3_allPeaksSet.rds")
  seqAll <- paste0(allPeaksSet@seqnames, allPeaksSet@ranges)
  
  sams <- read.table("P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakEach/P4_8_1_2_samples.txt")$V1
  # sams <- sams[1:5]
  
  # for (s in sams) {
  #   overlapMatrix <- readRDS(paste0("P4_8TFNet/P4_8_1GetNet/P4_8_1_4_OverlapPeak/P4_8_1_4_", s, "overlapMatrix.rds"))
  #   print(overlapMatrix[1:5, 1:5])
  # }
  # print("Read overlapMatrix")
  # overlapList <- readRDS(paste0("P4_8TFNet/P4_8_1GetNet/P4_8_1_4_OverlapPeak/P4_8_1_4_", s, "overlapMatrix.rds"))
  # # overlapMatrix <- do.call(cbind, overlapList)
  # mcols(allPeaksSet) <- overlapMatrix
  
  print("Read filter matrix")
  conList <- lapply(sams, function(s)readRDS(paste0("P4_8TFNet/P4_8_1GetNet/P4_8_1_4_OverlapPeak/P4_8_1_4_", s, "con.peak.rds")))
  names(conList) <- sams
  dat <- do.call(cbind, conList)
  dat[1:5, 1:5]

  print("Filter peaks")
  occurrences <- rowSums(dat)
  seq <- names(occurrences)[which(occurrences >= 2)]

  peakBound <- allPeaksSet[match(seq, seqAll), ]

  saveRDS(peakBound, file = fil.p)
} else {
  peakBound <- readRDS(fil.p)
}

### Annotate peaks
print("Annotate peaks")
peakMerge <- peakBound[, 0]
# i <- 1
fil.a <- paste0(out.dir, out.prefix, "peakAnn.rds")
if (!file.exists(fil.a)) {
  peakAnn <-
    annotatePeak(peakMerge,
                 TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                 annoDb = "org.Hs.eg.db",
                 tssRegion = c(-5000, 5000),
                 genomicAnnotationPriority = "Promoter",
                 sameStrand = T,
                 verbose = F
    )
  saveRDS(peakAnn, file = fil.a)
} else {
  peakAnn <- readRDS(fil.a)
}

### get peak counts
print("Get peak counts")
fil.c <- paste0(out.dir, out.prefix, "peakCounts.rds")
if (!file.exists(fil.c)) {
  peakCounts <- summarizeOverlaps(peakMerge, bamsToCount, singleEnd = FALSE)
  saveRDS(peakCounts, file = fil.c)
} else {
  peakCounts <- readRDS(fil.c)
}

print("Norm peak counts")
fil.n <- paste0(out.dir, out.prefix, "conATAC.rds")
if (!file.exists(fil.n)) {
  conATAC <- assay(peakCounts)
  strs <- paste0(peakCounts@rowRanges@seqnames, peakCounts@rowRanges@ranges)
  po <- which(!duplicated(strs))
  conATAC <- conATAC[po, ]
  rownames(conATAC) <- strs[po]
  colnames(conATAC) <- substr(colnames(conATAC), 10, 14)

  # summary(colSums(conATAC))
  # summary(rowSums(conATAC))

  sam <- colnames(conATAC)[which(colSums(conATAC) > 0)]
  sam <- intersect(sam, colnames(cnt.norm))

  print(dim(conATAC))
  # conATAC <- log10(conATAC[, sam] + 1)
  conATAC <- varianceStabilizingTransformation(as.matrix(conATAC[, sam])) # vst normalized readcounts from DESeq2_1.22.2
  dim(conATAC)
  conATAC[1:5, 1:5]
  saveRDS(conATAC, file = fil.n)
} else {
  conATAC <- readRDS(fil.n)
}


### peak-gene condidates
print("Get Correlations")
sam <- intersect(colnames(cnt.norm), colnames(conATAC))
fil.p <- paste0(out.dir, out.prefix, "peak2gene.rds")
if (!file.exists(fil.p)) {
  peak2gene <-
    data.frame(Peak = paste0(peakAnn@anno@seqnames, peakAnn@anno@ranges),
               Gene = peakAnn@anno$SYMBOL)
  # peak2gene
  print(dim(peak2gene))
  peaks <- intersect(rownames(conATAC), peak2gene$Peak)
  genes <- intersect(rownames(cnt.norm), peak2gene$Gene)
  peak2gene <- subset(peak2gene, Peak %in% peaks & Gene %in% genes)

  peak2gene$r <- NA
  peak2gene$p <- NA
  for(i in 1:nrow(peak2gene)) {
    if(i%%1000 == 0){
      print(paste0("Get cor: ", i, " of ", nrow(peak2gene)))
    }
    x <- conATAC[peak2gene$Peak[i], sam]
    y <- cnt.norm[peak2gene$Gene[i], sam]
    cor.res <- psych::corr.test(x, y)

    peak2gene$r[i] <- cor.res$r
    peak2gene$p[i] <- cor.res$p
  }

  saveRDS(peak2gene, file = fil.p)
} else {
  peak2gene <- readRDS(fil.p)
}




