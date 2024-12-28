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
out.dir <- "P4_8TFNet/P4_8_2TFRank/P4_8_2_3Expression/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_8_2_3_"


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

cnt <- readRDS("P4_2RNA/P4_2_2DEG/P4_2_2_1_counts.used.rds")

# NPM dif -----------------------------------------------------------------

resN2PM <- readRDS("P4_2RNA/P4_2_2DEG/P4_2_2_1_NormalRes.rds")

tfs <- intersect(unique(reg.net$TF), rownames(resN2PM))
resN2PM <- as.data.frame(resN2PM)[tfs, ]
fc <- resN2PM$log2FoldChange
Ediff_N <- -log10(resN2PM$pvalue + 10^(-24))*(fc/abs(fc))
names(Ediff_N) <- tfs
write.table(Ediff_N, file = paste0(out.dir, out.prefix, "E.diffNvsPM.txt"), quote = F, row.names = T, col.names = F)

resP2N <- readRDS("P4_2RNA/P4_2_3DEG1v1/P4_2_3_1_ResP2N.rds")

tfs <- intersect(unique(reg.net$TF), rownames(resP2N))
resP2N <- as.data.frame(resP2N)[tfs, ]
fc <- resP2N$log2FoldChange
Ediff_P <- -log10(resP2N$pvalue + 10^(-24))*(fc/abs(fc))
names(Ediff_P) <- tfs
write.table(Ediff_P, file = paste0(out.dir, out.prefix, "E.diffPvsN.txt"), quote = F, row.names = T, col.names = F)

resM2P <- readRDS("P4_2RNA/P4_2_3DEG1v1/P4_2_3_1_ResM2P.rds")

tfs <- intersect(unique(reg.net$TF), rownames(resM2P))
resM2P <- as.data.frame(resM2P)[tfs, ]
fc <- resM2P$log2FoldChange
Ediff_M <- -log10(resM2P$pvalue + 10^(-24))*(fc/abs(fc))
names(Ediff_M) <- tfs
write.table(Ediff_M, file = paste0(out.dir, out.prefix, "E.diffMvsP.txt"), quote = F, row.names = T, col.names = F)

Ediff.lis.npm <- list(NvsPM = Ediff_N, PvsN = Ediff_P, MvsP = Ediff_M)

saveRDS(Ediff.lis.npm, file = paste0(out.dir, out.prefix, "Ediff.lis_NPM.rds"))

# Clu diff ----------------------------------------------------------------

Ediff.lis <- list()
for (clu in unique(Cluster)) {
  res <- readRDS(paste0("P4_7Clustering/P4_7_4DiffClu/P4_7_4_1_Res", clu, ".rds"))
  tfs <- intersect(unique(reg.net$TF), rownames(res))
  res <- as.data.frame(res)[tfs, ]
  fc <- res$log2FoldChange
  Ediff <- -log10(res$pvalue + 10^(-24))*(fc/abs(fc))
  names(Ediff) <- tfs
  write.table(Ediff, file = paste0(out.dir, out.prefix, "E.diff", clu, ".txt"), quote = F, row.names = T, col.names = F)
  
  Ediff.lis[[clu]] <- Ediff
}

saveRDS(Ediff.lis, file = paste0(out.dir, out.prefix, "Ediff.lis_clu.rds"))








