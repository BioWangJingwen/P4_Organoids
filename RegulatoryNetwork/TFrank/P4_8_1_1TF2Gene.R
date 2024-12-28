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
out.dir <- "P4_8TFNet/P4_8_1GetNet/P4_8_1_1_TF2Gene/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_8_1_1_"

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


# for(f in bamsToCount) {
#   if (file.size(f) == 0) {
#     print(f)
#   }
# }

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

fil.p <- paste0(out.dir, out.prefix, "peakBound.rds")
if (!file.exists(fil.p)) {
  overlapMatrix <- do.call(cbind, overlap)
  mcols(allPeaksSet_nR) <- overlapMatrix
  
  seqAll <- paste0(allPeaksSet_nR@seqnames, allPeaksSet_nR@ranges)
  
  dat <- as.data.frame(elementMetadata(allPeaksSet_nR))
  rownames(dat) <- seqAll
  
  dat <- aggregate(t(dat), by = list(substr(colnames(dat), 1, 5)), FUN = function(x){y=sum(x);if (y>0){y=1}; return(y)})
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1]
  
  occurrences <- colSums(dat)
  seq <- names(occurrences)[which(occurrences >= 2)]
  
  peakBound <- allPeaksSet_nR[match(seq, seqAll), ]
  
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

print("Get TF 2 Peak ... ")
### Get peak2TF
index <- 1
# mat.sig <- subset(peak2gene, p < 0.01)
# peakSig <- peakBound[which(paste0(peakBound@seqnames, peakBound@ranges) %in% mat.sig$Peak), ]
peakTF <- peakBound@elementMetadata
rownames(peakTF) <- paste0(peakBound@seqnames, peakBound@ranges)
peakTF <- as.matrix(peakTF)
for(i in which(colSums(peakTF) > 0)) {
  dat <-
    data.frame(Peak = rownames(peakTF)[peakTF[, i]],
               TF = gsub(".*_", "", colnames(peakTF)[i]),
               SampleID = gsub("_.*", "", colnames(peakTF)[i]))
  if (index == 1) {
    index <- index + 1
    dat.TF2Peak <- dat
  } else {
    dat.TF2Peak <- rbind(dat, dat.TF2Peak)
  }
}
saveRDS(dat.TF2Peak, file =  paste0(out.dir, out.prefix, "dat.TF2Peak.rds"))

reg.net <- join(dat.TF2Peak, peak2gene, by = "Peak")

saveRDS(reg.net, file =  paste0(out.dir, out.prefix, "reg.net.rds"))
