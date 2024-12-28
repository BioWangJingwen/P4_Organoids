rm(list = ls())

# library(ATACseqQC)
library(ChIPseeker)
library(ChIPpeakAnno)
library(idr)
library(reticulate)
library(DiffBind)
library(tidyverse)
library(purrr)
library(parallel)
library(GenomicRanges)
library(GenomicAlignments)

# use_condaenv(condaenv = "/chenfeilab/Avocado/Softwares/miniconda3/", required = NULL)

setwd("/chenfeilab/Avocado/P4_Organoids")
wd <- "/chenfeilab/Avocado/P4_Organoids/"

out.prefix <- "P4_4_4_2_"
out.dir <- paste0("P4_4ATAC/P4_4_4peakDiff/", out.prefix, "peakDiff/")
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}

# All Info
samInfoATAC <- readRDS("P4_4ATAC/P4_4_1QC/P4_4_1_1_samInfoATAC.rds")
head(samInfoATAC)

sam.all <- table(samInfoATAC$SampleID)

samMergePeak <- readRDS("P4_4ATAC/P4_4_2Peaks/P4_4_2_1_samMergePeak.rds")
samInfoFilter <- readRDS("P4_4ATAC/P4_4_2Peaks/P4_4_2_1_samInfoFilter.rds")

# colors
source("P4_1Data/P4_1_0Colors/P4_1_0colors.R")

# peak counts
myCounts <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/P4_4_4_1_PeakCounts.rds")
myCounts[1:5, 1:5]
colnames(myCounts) <- gsub("^P4_4_2_1_", "", gsub(".final.bam$", "", colnames(myCounts)))
# samples has zero expression
sam <- names(which(colSums(assay(myCounts)) == 0))
usedSam <- setdiff(colnames(myCounts), sam)
myCounts <- myCounts[, usedSam]

# Diff peak ---------------------------------------------------------------

# get duplication of peaks of DESeq
N2PM <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/P4_4_4_1_atacDDSNormal.rds")
P2N <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/P4_4_4_1_atacDDSP2N.rds")
P2M <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/P4_4_4_1_atacDDSP2M.rds")
M2N <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/P4_4_4_1_atacDDSM2N.rds")
M2P <- readRDS("P4_4ATAC/P4_4_4peakDiff/P4_4_4_1peakCounts/P4_4_4_1_atacDDSM2P.rds")
myPeaks <- 
  list(NvsPM = N2PM,
       PvsN = P2N,
       PvsM = P2M,
       MvsN = M2N,
       MvsP = M2P
  )
myPeaks <- 
  lapply(myPeaks, function(x){
    sig <- results(x, format = "GRanges")
    subset(sig, log2FoldChange > 1 & padj < 0.01)
  })
names(myPeaks)

saveRDS(myPeaks, file = paste0(out.dir, out.prefix, "sigPeaks.rds"))

# lapply(peaks, ChIPQC::GetGRanges, simple = TRUE)
# names(myPeaks) <- samMergePeak$SampleID
allPeaksSet_nR <- reduce(unlist(GRangesList(myPeaks)))
# allPeaksSet_nR <- GRangesList(myPeaks)
allPeaksSet_nR

overlap <- list()
for (i in 1:length(myPeaks)) {
  overlap[[i]] <- allPeaksSet_nR %over% myPeaks[[i]]
}
names(overlap) <- names(myPeaks)

overlapMatrix <- do.call(cbind, overlap)
mcols(allPeaksSet_nR) <- overlapMatrix
colSums(overlapMatrix)

allPeaksSet_nR[1:2, ]

saveRDS(allPeaksSet_nR, file = paste0(out.dir, out.prefix, "SigPeak.rds"))

# Venn plot of sig. peaks
library(limma)
vennDiagram(mcols(allPeaksSet_nR), 
            circle.col = col.deg[colnames(mcols(allPeaksSet_nR))],
            
)

pdf(paste0(out.dir, out.prefix, "VennPeaks.pdf"), width = 4.5, height = 4.5)
vennDiagram(mcols(allPeaksSet_nR), 
            circle.col = col.deg[colnames(mcols(allPeaksSet_nR))],
            cex=c(1.2,1,0.7), lwd = 1.5
)
dev.off()


# peak annotation ---------------------------------------------------------

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
peakAnno <- # annoPeaks(myPeaks$NvsPM)
  lapply(myPeaks, function(x){
    annotatePeak(x, tssRegion = c(-500, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                 annoDb = "org.Hs.eg.db")
  })
saveRDS(peakAnno, file = paste0(out.dir, out.prefix, "peakAnnoPos.rds"))

for (i in 1:length(peakAnno)) {
  sta <- peakAnno[[i]]@annoStat
  sta$Group <- names(peakAnno)[i]
  if (i == 1) {
    peakAnno.dat <- sta
  } else {
    peakAnno.dat <- rbind(peakAnno.dat, sta)
  }
}
peakPlot <- 
  ggplot(peakAnno.dat, aes(Group, Frequency, fill = Feature)) +
  geom_col() +
  labs(x = "", y = "Percentage(%)") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  scale_fill_d3(alpha = 0.62) +
  # facet_grid(Source~., scales = "free", space = "free") +
  theme_classic() +
  theme(axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text.x = element_text(size = 18), axis.title = element_text(size = 26),
        axis.text.y = element_text(size = 18),
        # strip.text = element_text(size = 21), 
        legend.title = element_text(size = 13), legend.text = element_text(size = 9)
  ) 

pdf(paste0(out.dir, out.prefix, "PositionPeaks.pdf"), width = 7, height = 3.5)
peakPlot
dev.off()

# # add gene id
# peakAnno <- lapply(peakAnno, function(x){
#   addGeneIDs(x, "org.Hs.eg.db", IDs2Add = c("symbol", "entrez_id"))
# })

# EnrichGO
g.enGO <- list()
for (g in names(peakAnno)) {
  print(g)
  x <- unique(na.omit(peakAnno[[g]]@anno$ENSEMBL))
  
  resGO <-
    getEnrichedGO(
      x,
      orgAnn="org.Hs.eg.db", 
      maxP=0.01,
      minGOterm=10,
      multiAdjMethod= NULL
    )
  saveRDS(resGO, file = paste0(out.dir, out.prefix, "resGO", g, ".rds"))
  g.enGO[[g]] <-
    enrichmentPlot(resGO, n = 10) +
    scale_fill_continuous( low = alpha(col.deg[g], alpha = 0.05),
                           high = col.deg[g])
}

pdf(paste0(out.dir, out.prefix, "enrichGO.pdf"), width = 9, height = 12)
ggarrange(plotlist = g.enGO, 
          # common.legend = T, 
          nrow = 3, ncol = 2,
          align = "hv",
          # label.x = 0,
          # label.y = 1,
          # hjust = 0,
          # vjust = 0,
          labels = names(peakAnno)
          )
dev.off()

# EnrichKEGG
# #安装Y叔的包，
# #安装创建KEGG数据库的包的包
# remotes::install_github("YuLab-SMU/createKEGGdb")
# #创建自己的物种的包create_kegg_db，会自动创建名称为KEGG.db_1.0.tar,gz的包。物种名称的简写，在
# createKEGGdb::create_kegg_db('hsa')
# 
# #安装这个包(默认的包的路径在当前工作目录，根据实际情况修改路径)
# install.packages("/chenfeilab/Avocado/P4_Organoids/KEGG.db_1.0.tar.gz",repos=NULL,type="source")

library(KEGG.db)
library(org.Hs.eg.db)
library(ggplot2)
hsa_kegg <- clusterProfiler::download_KEGG("hsa")
names(hsa_kegg)
head(hsa_kegg$KEGGPATHID2NAME)

g.enKEGG <- list()
for (g in names(peakAnno)) {
  print(g)
  x <- unique(na.omit(peakAnno[[g]]@anno$ENSEMBL))
  
  resKEGG <-
    getEnrichedPATH(
      x,
      orgAnn="org.Hs.eg.db", 
      pathAnn = "KEGG.db",
      maxP=0.01,
      multiAdjMethod= NULL
    )
  resKEGG$path.term <- plyr::mapvalues(resKEGG$path.id, hsa_kegg$KEGGPATHID2NAME$from,  hsa_kegg$KEGGPATHID2NAME$to)
  
  saveRDS(resKEGG, file = paste0(out.dir, out.prefix, "resKEGG", g, ".rds"))
  g.enKEGG[[g]] <-
    enrichmentPlot(resKEGG, n = 10) +
    scale_fill_continuous( low = alpha(col.deg[g], alpha = 0.05),
                           high = col.deg[g])
}

length(g.enKEGG)
library(ggpubr)
pdf(paste0(out.dir, out.prefix, "enrichKEGG.pdf"), width = 5, height = 12)
ggarrange(plotlist = g.enKEGG, 
          # common.legend = T, 
          nrow = 3, ncol = 2,
          align = "hv",
          # label.x = 0,
          # label.y = 1,
          # hjust = 0,
          # vjust = 0,
          labels = names(peakAnno)
)
dev.off()

# reactome.db
library(reactome.db)
# hsa_REACTOME <- clusterProfiler::download_REACTOME("hsa")
# names(hsa_REACTOME)
# head(hsa_REACTOME$REACTOMEPATHID2NAME)

g.enREACTOME <- list()
for (g in names(peakAnno)) {
  print(g)
  x <- unique(na.omit(peakAnno[[g]]@anno$ENSEMBL))
  
  resREACTOME <-
    getEnrichedPATH(
      x,
      orgAnn="org.Hs.eg.db", 
      pathAnn = "reactome.db",
      maxP=0.01,
      multiAdjMethod= NULL
    )
  resREACTOME$path.term <- gsub("^Homo sapiens: ", "", resREACTOME$path.term)
  
  saveRDS(resREACTOME, file = paste0(out.dir, out.prefix, "resREACTOME", g, ".rds"))
  g.enREACTOME[[g]] <-
    enrichmentPlot(resREACTOME, n = 10) +
    scale_fill_continuous( low = alpha(col.deg[g], alpha = 0.05),
                           high = col.deg[g])
}

length(g.enREACTOME)
pdf(paste0(out.dir, out.prefix, "enrichREACTOME.pdf"), width = 5, height = 12)
ggarrange(plotlist = g.enREACTOME, 
          # common.legend = T, 
          nrow = 3, ncol = 2,
          align = "hv",
          # label.x = 0,
          # label.y = 1,
          # hjust = 0,
          # vjust = 0,
          labels = names(peakAnno)
)
dev.off()
