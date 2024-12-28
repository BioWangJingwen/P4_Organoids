# JobID 33662
rm(list = ls())

library(chromVAR)
library(BiocParallel)
register(MulticoreParam(8)) # Use 8 cores
library(data.table)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("/chenfeilab/Avocado/P4_Organoids/")
wd <- "/chenfeilab/Avocado/P4_Organoids/"

out.dir <- "P4_8TFNet/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.dir <- "P4_8TFNet/P4_8_2TFRank/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.dir <- "P4_8TFNet/P4_8_2TFRank/P4_8_2_2Accessibility/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_8_2_2_"


# read data ---------------------------------------------------------------

reg.net <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_7_reg.netR0.6FDR0.01.rds")
head(reg.net)

# Get clusters
aa <- load("P4_7Clustering/P4_7_2ATACClusters/P4_7_2_1_res.Primary_Metastatic_Clusters.RData")
aa
optK <- 4
results <- rcc3
Cluster <- paste0("Cluster", results[[optK]]$consensusClass)
names(Cluster) <- names(results[[optK]]$consensusClass)
write.table(data.frame(PatientID = names(Cluster),
                       Cluster = Cluster
                       ),
            file = paste0(out.dir, out.prefix, "Cluster.csv"), sep = ",", quote = F, col.names = T, row.names = F)

peakBound <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_5_Peak2Gene/P4_8_1_5_peakBound.rds")
seqBound <- paste0(peakBound@seqnames, peakBound@ranges)
# peakBound <- peakBound[which(seqBound %in% reg.net$Peak), ]
# seqBound <- paste0(peakBound@seqnames, peakBound@ranges)

bamfiles <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakEach/P4_8_1_2_bamsToCount.rds")
bamnames <- gsub(".*_", "", gsub("\\..*", "", basename(bamfiles)))
# bamfiles <- bamfiles[which(bamnames %in% reg.net$SampleID)]
# bamnames <- gsub(".*_", "", gsub("\\..*", "", basename(bamfiles)))


# motif 2 tf --------------------------------------------------------------

motif2gene <- 
  read.table("/chenfeilab/Avocado/MyPipline/TOBIAS/test_data/motif2gene_mapping.txt", sep = "\t")
colnames(motif2gene) <- c("Motif", "ENSEMBL")
head(motif2gene)
# This is a file containing information of the origin gene for each TF. 
# By "origin", what is meant is the gene from which the transcription factor is coded. 
# The example-file 'motif2gene_mapping.txt' within the test data gives a good impression of how such a file should look.
library(org.Hs.eg.db)
library(plyr)
trs <- clusterProfiler::bitr(motif2gene$ENSEMBL, OrgDb = "org.Hs.eg.db", fromType = "ENSEMBL", toType = "SYMBOL")
head(trs)
motif2tf <- join(motif2gene, trs, by = "ENSEMBL")
motif2tf <- na.omit(motif2tf)
head(motif2tf)

# Get counts --------------------------------------------------------------

fil <- "P4_8TFNet/P4_8_2TFRank/P4_8_2_2Accessibility/P4_8_2_2_scores.rds"

if (!file.exists(fil)) {
  # dat.col <-
  #   DataFrame(Source = substr(bamnames, 5, 5),
  #             Cluster = "NotClu",
  #             row.names = bamnames
  #   )
  # intnames <- intersect(bamnames, names(Cluster))
  # dat.col[intnames, "Cluster"] <- Cluster[intnames]
  
  fil.cnt <- paste0(out.dir, out.prefix, "fragment_counts.rds")
  if (!file.exists(fil.cnt)) {
    fragment_counts <- 
      getCounts(bamfiles, peakBound, 
                paired =  TRUE, 
                by_rg = FALSE, 
                format = "bam"
                # colData = dat.col
      )
    
    saveRDS(fragment_counts, file = fil.cnt)
  } else {
    fragment_counts <- readRDS(fil.cnt)
  }
  
  
  # GC context
  GC_counts <- addGCBias(fragment_counts ,
                         genome = BSgenome.Hsapiens.UCSC.hg19)
  
  # 过滤peaks，也有其他option.
  counts_filtered <- filterSamples(GC_counts, min_depth = 1500,
                                   min_in_peaks = 0.15)
  counts_filtered <- filterPeaks(counts_filtered)
  
  # match motif
  # motifs <- getJasparMotifs()
  motifs <- TFBSTools::readJASPARMatrix("P4_1Data/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt")
  motifs <- motifs[which(names(motifs) %in% motif2tf$Motif)]
  saveRDS(motifs, file = paste0(out.dir, out.prefix, "JasparMotifs.rds"))
  
  motif_ix <- matchMotifs(motifs, counts_filtered,
                          genome = BSgenome.Hsapiens.UCSC.hg19)
  
  # computing deviations（这一步计算很久）
  # "z-score": z-score也称之为"deviation score"
  dev <- computeDeviations(object = counts_filtered,
                           annotations = motif_ix)
  
  saveRDS(dev, file = paste0(out.dir, out.prefix, "dev.rds"))
  
  
  scores <- deviationScores(dev)
  colnames(scores) <- gsub("P.*_", "", gsub("\\..*", "", colnames(scores)))
  saveRDS(scores, file = paste0(out.dir, out.prefix, "scores.rds"))
  write.table(scores, file = paste0(out.dir, out.prefix, "scores.csv"), sep = ",", quote = F, row.names = T, col.names = T)
} else {
  scores <- readRDS(fil)
}

rownames(scores) <- 
  mapvalues(rownames(scores), motif2tf$Motif, motif2tf$SYMBOL)
scores <- scores[intersect(rownames(scores), reg.net$TF), ]

# A_Diff NPM ----------------------------------------------------------------

samN <- grep("N$", colnames(scores), value = T)
samP <- grep("T$", colnames(scores), value = T)
samM <- grep("M$", colnames(scores), value = T)

A.diff.NPM.lis <- list()

### NvsPM
out.N <- rowSums(scores[, samN])/length(samN)
out.PM <- rowSums(scores[, c(samP, samM)])/length(c(samP, samM))
A.NvsPM <- out.N - out.PM
summary(A.NvsPM); sort(A.NvsPM, decreasing = T)[1:10]
A.diff.NPM.lis[["NvsPM"]] <- A.NvsPM
write.table(A.NvsPM, file = paste0(out.dir, out.prefix, "A.diffNvsPM.txt"), quote = F, row.names = T, col.names = F)

### PvsN
out.N <- rowSums(scores[, samN])/length(samN)
out.P <- rowSums(scores[, samP])/length(samP)
A.PvsN <- out.P - out.N
summary(A.PvsN); sort(A.PvsN, decreasing = T)[1:10]
A.diff.NPM.lis[["PvsN"]] <- A.PvsN
write.table(A.PvsN, file = paste0(out.dir, out.prefix, "A.diffPvsN.txt"), quote = F, row.names = T, col.names = F)

### MvsP
out.P <- rowSums(scores[, samP])/length(samP)
out.M <- rowSums(scores[, samM])/length(samM)
A.MvsP <- out.M - out.P
summary(A.MvsP); sort(A.MvsP, decreasing = T)[1:10]
A.diff.NPM.lis[["MvsP"]] <- A.MvsP
write.table(A.MvsP, file = paste0(out.dir, out.prefix, "A.diffMvsP.txt"), quote = F, row.names = T, col.names = F)

saveRDS(A.diff.NPM.lis, file = paste0(out.dir, out.prefix, "A.diff.NPM.lis.rds"))

# A_diff Clusters ---------------------------------------------------------

Cluster <- Cluster[intersect(names(Cluster), colnames(scores))]
A.diff.clu.lis <- list()
for (clu in unique(Cluster)) {
  # clu <- "Cluster4"
  sam1 <- names(Cluster)[which(Cluster == clu)]
  sam2 <- names(Cluster)[which(Cluster != clu)]
  out1 <- rowSums(scores[, sam1])/length(sam1)
  out2 <- rowSums(scores[, sam2])/length(sam2)
  A.diff <- out1 - out2
  A.diff.clu.lis[[clu]] <- A.diff
  
  write.table(A.diff, file = paste0(out.dir, out.prefix, "A.diff", clu, ".txt"), quote = F, row.names = T, col.names = F)
  
}

saveRDS(A.diff.clu.lis, file = paste0(out.dir, out.prefix, "A.diff.clu.lis.rds"))









