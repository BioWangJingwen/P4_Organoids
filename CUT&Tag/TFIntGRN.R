rm(list = ls())

setwd("/chenfeilab/Avocado/P4_Organoids/")
wd <- "/chenfeilab/Avocado/P4_Organoids/"

library(corrplot)
library(ggplot2)
library(GenomicRanges)
library(ChIPseeker)
library(ChIPpeakAnno)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ucsc.hg19.knownGene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(annotatr)
library(reticulate)
use_condaenv(condaenv = "/work/chenfeilab/Avocado/Softwares/miniconda3/", required = NULL)
library(plyr)

out.dir <- "P4_9CutNTag/P4_9_5EP/P4_9_5_4TF_GRN/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = T)
}
out.prefix <- "P4_9_5_4_"

source("/chenfeilab/Avocado/P4_Organoids/P4_1Data/P4_1_0Colors/P4_1_0colors.R")

# Get cluster
aa <- load("P4_7Clustering/P4_7_2ATACClusters/P4_7_2_1_res.Primary_Metastatic_Clusters.RData")
aa
optK <- 4
results <- rcc3
Cluster <- paste0("Cluster", results[[optK]]$consensusClass)
names(Cluster) <- names(results[[optK]]$consensusClass)

# get ann -----------------------------------------------------------------

anndir <- "P4_9CutNTag/P4_9_5EP/P4_9_5_2EP_Ana/"

# ann with ChIPSeeker
fil.ChIPSeeker <- list.files(anndir, pattern = "_Ann_ChIPSeeker.rds")
fil.ChIPSeeker <- grep("IgG", fil.ChIPSeeker, value = T, invert = T)
fil.ChIPSeeker <- grep("H3K27me3", fil.ChIPSeeker, value = T, invert = T)
length(fil.ChIPSeeker)
sam <- gsub("P4_9_5_2_", "", gsub("_Ann_ChIPSeeker.rds", "", fil.ChIPSeeker))
sam <- unique(gsub(".*_", "", sam))
Cluster[sam]

net <- readRDS("/chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_1GetNet/P4_8_1_7_reg.netR0FDR0.01.rds")
bg <- unique(net$Gene)
head(net)

anti <- c("EBF", "MIST1", "ZBP89")
gene <- c("EBF3", "BHLHA15", "ZNF148")
for (s in sam) {
  print(s)
  # s <- sam[1]

  fils <- grep(s, fil.ChIPSeeker, value = T)
  TFs <- gsub("P4_9_5_2_", "", gsub(paste0("_", s, "_Ann_ChIPSeeker.rds"), "", fils))
  TFs <- mapvalues(TFs, anti, gene)
  print(TFs[which(!(TFs %in% net$TF))])
  
  net.int <- subset(net, SampleID == s & TF %in% TFs & Gene %in% TFs)
  if (nrow(net.int) > 0) {
    write.table(net.int, file = paste0(out.dir, out.prefix, s, "net.txt"),
                sep = "\t", quote = F, col.names = T, row.names = F)
  }
  
  net.s <- subset(net, SampleID == s & TF %in% TFs)
  TF.use <- unique(net.s$TF)
  lis.net <- lapply(TF.use, function(tf){
    aa = subset(net.s, TF == tf)
    return(unique(aa$Gene))
  })
  names(lis.net) <- TF.use
  
  mat.ji <- 
    matrix(0, nrow = length(TF.use), ncol = length(TF.use),
           dimnames = list(TF.use, TF.use))
  diag(mat.ji) <- 1
  
  dat.int <- t(combn(TF.use, 2))
  colnames(dat.int) <- c("TF1", "TF2")
  # dat.int <- as.data.frame(dat.int)
  # dat.int$Jaccard <- NA
  # dat.int$p <- NA
  for (i in 1:nrow(dat.int)) {
    # i <- 1
    t1 <- dat.int[i, 1]
    t2 <- dat.int[i, 2]
    g1 <- lis.net[[t1]]
    g2 <- lis.net[[t2]]
    
    # dat.int$Jaccard[i] <-
    mat.ji[t1, t2] <- mat.ji[t2, t1] <- 
      length(intersect(g1, g2))/length(union(g1, g2))
    
    # dat.int$p[i] <-
    #   1 - phyper(length(intersect(g1, g2)) - 1, length(setdiff(bg, g1)), length(g1), length(g2))
  }
  
  mat.cor <- cor(mat.ji, method = "spearman")
  
  save(mat.ji, mat.cor, file = paste0(out.dir, out.prefix, s, "mat.int.RData"))
  
  # heatmap(mat.ji)
  # heatmap(mat.cor)
  
  mat <- mat.ji
  mat[upper.tri(mat)] <- mat.cor[upper.tri(mat.cor)]
  rownames(mat) <- mapvalues(rownames(mat), gene, anti)
  colnames(mat) <- mapvalues(colnames(mat), gene, anti)
  
  # lower triangle is jaccard; upper is spearman of jaccard
  pdf(paste0(out.dir, out.prefix, s, "Ji_cor.pdf"), width = 6, height = 6)
  corrplot.mixed(mat, lower = 'number', upper = 'circle', order = 'original',
                 lower.col = hcl.colors(52, "Peach", rev = TRUE), 
                 upper.col = hcl.colors(52, "Blue-Red2"), main = s,
                 mar = c(5, 4, 4, 2) + 1
                 )
  dev.off()
}











