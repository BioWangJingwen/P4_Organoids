# JobID 33662
rm(list = ls())
library(ggrepel)
library(ggthemes)
library(ggpubr)

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

out.prefix <- "P4_8_2_5_"

source("P4_8TFNet/P4_8_2_4RankFUN.R")
source("P4_1Data/P4_1_0Colors/P4_1_0colors.R")

# Get clusters
aa <- load("P4_7Clustering/P4_7_2ATACClusters/P4_7_2_1_res.Primary_Metastatic_Clusters.RData")
aa
optK <- 4
results <- rcc3
Cluster <- paste0("Cluster", results[[optK]]$consensusClass)
names(Cluster) <- names(results[[optK]]$consensusClass)

SelTF <- c("MEF2A", "NR2F2", "BHLHA15", "NFATC4", "NR2F1", "ZNF148", "TFAP4")

# Rank Cluster ------------------------------------------------------------

O_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_1OutDegree/P4_8_2_1_O.diff.clu.lis.rds")
A_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_2Accessibility/P4_8_2_2_A.diff.clu.lis.rds")
E_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_3Expression/P4_8_2_3_Ediff.lis_clu.rds")

ord_mat <- read.table("P4_8TFNet/P4_8_2TFRank/P4_8_2_4_ord_mat_Clusters.txt", row.names = 1, header = T, sep = "\t")

lis.dif <- list(
  A_diff = A_diff,
  E_diff = E_diff,
  O_diff = O_diff
)

g.lis.clu <- list()
for (clu in paste0("Cluster", 4)) {
  tfVec <- ord_mat[, clu]
  names(tfVec) <- rownames(ord_mat)
  # tfR <- sort(tfVec, decreasing = T)[1:5]
  
  for (i in 1:length(lis.dif)) {
    I_diff <- lis.dif[[i]]
    
    O_dat <-
      data.frame(I_diff = I_diff[[clu]],
                 Gene = names(I_diff[[clu]]))
    O_dat$Gene <- factor(O_dat$Gene, levels = O_dat$Gene[order(O_dat$I_diff, decreasing = T)])
    O_dat$Label <- NA
    po <- which(O_dat$Gene %in% SelTF)
    O_dat$Label[po] <- as.character(O_dat$Gene[po])
    
    col <- rep(col.level[names(lis.dif)[i]], length(po))
    names(col) <- O_dat$Label[po]
    
    g.lis.clu[[paste0(clu, names(lis.dif)[i])]] <-
      ggplot(O_dat, aes(Gene, I_diff, label = Label, col = Label)) +
      geom_text_repel() +
      geom_point() +
      theme_base() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none") +
      scale_color_manual(values = col, na.value = "#cccccc66") +
      scale_x_discrete(expand = c(.05, .05)) +
      labs(y = names(lis.dif)[i], title = clu)
  }
}

gg <-
  ggarrange(plotlist = g.lis.clu, 
            ncol = 3, nrow = 1,
            align = "hv"
  )
gg

pdf(paste0(out.dir, out.prefix, "Index_Cluster4.pdf"), width = 7, height = 3.6)
print(gg)
dev.off()

g.lis.clu <- list()
for (clu in paste0("Cluster", 1:4)) {
  tfVec <- ord_mat[, clu]
  names(tfVec) <- rownames(ord_mat)
  tfR <- sort(tfVec, decreasing = T)[1:25]
  
  for (i in 1:length(lis.dif)) {
    I_diff <- lis.dif[[i]]
    
    O_dat <-
      data.frame(I_diff = I_diff[[clu]],
                 Gene = names(I_diff[[clu]]))
    O_dat$Gene <- factor(O_dat$Gene, levels = O_dat$Gene[order(O_dat$I_diff, decreasing = T)])
    O_dat$Color <- NA
    po <- which(O_dat$Gene %in% names(tfR))
    O_dat$Color[po] <- as.character(O_dat$Gene[po])
    
    col <- rep(col.level[names(lis.dif)[i]], length(po))
    names(col) <- O_dat$Color[po]
    
    O_dat$Label <- NA
    po <- which(O_dat$Gene %in% names(tfR)[1:10])
    O_dat$Label[po] <- as.character(O_dat$Gene[po])
    
    g.lis.clu[[paste0(clu, names(lis.dif)[i])]] <-
      ggplot(O_dat, aes(Gene, I_diff, label = Label, color = Color)) +
      geom_text_repel(max.overlaps = 25,
                      force_pull = 0,
                      # xlim = c(600, 700),
                      nudge_x = 605,
                      force = 0.5,
                      direction    = "y",
                      hjust        = 0,
                      segment.size = 0.2) +
      geom_point() +
      geom_point(subset(O_dat, !is.na(Color)), mapping = aes(Gene, I_diff, color = Color ))+
      theme_base() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none") +
      scale_color_manual(values = col, na.value = "#cccccc66") +
      scale_x_discrete(expand = c(.05, .05)) +
      labs(y = names(lis.dif)[i], title = clu)
  }
}

gg <-
  ggarrange(plotlist = g.lis.clu, 
            ncol = 6, nrow = 2,
            align = "hv"
  )
# gg

pdf(paste0(out.dir, out.prefix, "Index_Cluster.pdf"), width = 16, height = 6)
print(gg)
dev.off()

# Score -------------------------------------------------------------------

library(chromVAR)
library(BiocParallel)
register(MulticoreParam(8)) # Use 8 cores
library(data.table)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)

peakBound <- readRDS("P4_8TFNet/P4_8_1GetNet/P4_8_1_5_Peak2Gene/P4_8_1_5_peakBound.rds")
seqBound <- paste0(peakBound@seqnames, peakBound@ranges)
# peakBound <- peakBound[which(seqBound %in% reg.net$Peak), ]
# seqBound <- paste0(peakBound@seqnames, peakBound@ranges)

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

fil <- paste0(out.dir, out.prefix, "scores.rds")
samMergePeak <- readRDS("P4_4ATAC/P4_4_2Peaks/P4_4_2_1_samMergePeak.rds")
bamfiles <- gsub("peaks$", "final.bam", samMergePeak$peak)
bamfiles <- grep("N.final.bam", bamfiles, invert = T, value = T)

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
  
  fil.cf <- paste0(out.dir, out.prefix, "counts_filtered.rds")
  if (!file.exists(fil.cnt)) {
    # GC context
    GC_counts <- addGCBias(fragment_counts ,
                           genome = BSgenome.Hsapiens.UCSC.hg19)
    
    # 过滤peaks，也有其他option.
    counts_filtered <- filterSamples(GC_counts, min_depth = 1500,
                                     min_in_peaks = 0.15)
    counts_filtered <- filterPeaks(counts_filtered)
    
    saveRDS(counts_filtered, file = fil.cf)
  } else {
    counts_filtered <- readRDS(fil.cf)
  }

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
# scores <- scores[intersect(rownames(scores), reg.net$TF), ]

# CAP heat ----------------------------------------------------------------
### deviation score of TF motifs


reg.net <- readRDS("/chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_1GetNet/P4_8_1_7_reg.netR0.6FDR0.01.rds")
reg.net.sel <- subset(reg.net, TF %in% SelTF)
reg.net.sel$Index <- paste0(reg.net.sel$TF, reg.net.sel$Peak)
reg.net.sel <- subset(reg.net.sel, !duplicated(reg.net.sel$Index))
head(reg.net.sel)
dim(reg.net.sel)
length(reg.net.sel$Peak)
length(unique(reg.net.sel$Peak))

# scores <- readRDS("/chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_2TFRank/P4_8_2_2Accessibility/P4_8_2_2_scores.rds")
# rownames(scores) <-
#   mapvalues(rownames(scores), motif2tf$Motif, motif2tf$SYMBOL)
scores <- scores[intersect(rownames(scores), reg.net.sel$TF), ]
sam <- intersect(colnames(scores), names(Cluster))
conATAC.tf <- scores[, sam]

#
# conATAC <- #readRDS("/chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_2TFRank/P4_8_2_2Accessibility/P4_8_2_2_scores.rds")
#   readRDS("/chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_1GetNet/P4_8_1_5_Peak2Gene/P4_8_1_5_conATAC.rds")
# conATAC[1:5, 1:5]
#
# which(!(reg.net.sel$Peak %in% rownames(conATAC)))
#
# peak <- intersect(reg.net.sel$Peak, rownames(conATAC))
#
# reg.net.tf <- subset(reg.net.sel, Peak %in% peak)
# sam <- intersect(colnames(conATAC), names(Cluster))
# ll <- lapply(unique(reg.net.tf$TF), function(t){
#   p <- reg.net.sel$Peak[which(reg.net.sel$TF == t)]
#   m <- conATAC[p, sam]
#   rownames(m) <- paste0(t, p)
#   return(m)
# })
# mat.exp <- do.call(rbind, ll)
# mm <- rowMeans(mat.exp)
# ss <- rowSds(mat.exp)
# conATAC.tf <- (mat.exp - mm)/ss
# # mat.exp <- apply(mat.exp, 1, scales::rescale)
# # mat.exp[1:5, 1:5]
# # conATAC.tf <- t(conATAC.tf)
# conATAC.tf[1:5, 1:5]
# dim(conATAC.tf)

library(plyr)
dat.row <-
  data.frame(TFs = gsub("chr.*", "", rownames(conATAC.tf)),
             Peak = gsub(".*chr", "chr", rownames(conATAC.tf))
  )
rownames(dat.row) <- rownames(conATAC.tf)

library(ComplexHeatmap)
ann.row <- rowAnnotation(
  TFs = anno_block(gp = gpar(fill = col.TFSel))
  # Pathway = dat.row$Pathway
)

samples <- colnames(conATAC.tf)
annCol <- data.frame(Cluster = mapvalues(samples, names(Cluster), Cluster),
                     row.names = samples)

col_fun <- function(x, c) {
  circlize::colorRamp2(breaks = seq(from = min(x), to = max(x), length.out = 2), colors = c("white", c))
}

ann.col <-
  columnAnnotation(
    Cluster = anno_block(gp = gpar(fill = mycol))
  )
spl.row <- dat.row$TFs
spl.col <- annCol$Cluster

set.seed(123456)
h.DP <-
  Heatmap(conATAC.tf,
          col = circlize::colorRamp2(seq(1, 30, length.out = 9), pal_material("green")(9)),
          top_annotation = ann.col,
          # left_annotation = ann.row,
          name = "z-score",
          # row_title_rot = 0,
          show_row_names = T,
          show_row_dend = F,
          show_column_names = F,
          show_column_dend = F,
          # column_title_rot = 90,
          column_dend_reorder = F,
          cluster_column_slices = F,
          cluster_row_slices = F,
          # layer_fun = layer_fun,
          # row_split = spl.row,
          column_split = spl.col
  )

pdf(paste0(out.dir, out.prefix, "HeatClu.pdf"), width = 6, height = 2)
h.DP
dev.off()


# colorRamp2(seq(1, 30, length.out = 9), pal_material("green")(9))


# box ---------------------------------------------------------------------

source("/work/chenfeilab/Avocado/MyPipline/S12_ggplots.R")

lis.plot <- list()
for (t in rownames(conATAC.tf)) {
  print(t)
  # t <- "NR2F1"
  cl <- "Cluster4"
  
  dat.plot <- data.frame(
    x = factor(mapvalues(Cluster[colnames(conATAC.tf)], setdiff(unique(Cluster), "Cluster4"), c("Others", "Others", "Others")), 
               levels = c(cl, "Others")),
    y = as.numeric(conATAC.tf[t, ])
  )
  cc <- c(Cluster4 = unname(col.TFSel[t]), Others = "#eeeeee")
  # names(cc) <- c(t, "Others")
  lis.plot[[t]] <- 
    regBoxPlot(dat.plot, 
               col = cc, 
               title = t, ylab = "deviation Z-scores", xlab = NULL, limits = NULL)
  
}

ggarrange(plotlist = lis.plot, nrow = 2, ncol = 5)

pdf(paste0(out.dir, out.prefix, "BoxClu.pdf"), width = 9, height = 5)
ggarrange(plotlist = lis.plot, nrow = 2, ncol = 5)
dev.off()

# "#79af97"

dev <- readRDS("/chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_2TFRank/P4_8_2_2Accessibility/P4_8_2_2_dev.rds")

dev$Cluster <- Cluster[gsub("P4_4_2_1_", "", gsub(".final.bam", "", colnames(dev)))]
dev <- subset(dev, !is.na(Cluster))
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 10, 
                               shiny = FALSE)
# plotDeviationsTsne 绘图
tsne_plots <- plotDeviationsTsne(dev, tsne_results, #annotation = "TEAD3", 
                                 sample_column = "Cluster", shiny = FALSE)
#注意sample_column 名字要和colData中你的colname对应上。
tsne_plots[[1]] +
  scale_color_manual(values = mycol, na.value = "#FFFFFFFF")

# kmer_ix <- matchKmers(6, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)
# kmer_dev <- computeDeviations(counts_filtered, kmer_ix)
# kmer_cov <- deviationsCovariability(kmer_dev)
# plotKmerMismatch("CATTCC",kmer_cov)
# 
# de_novos <- assembleKmers(kmer_dev, progress = FALSE) #no progress bar
# de_novos
# 
# dist_to_known <- pwmDistance(de_novos, motifs)
# closest_match1 <- which.min(dist_to_known$dist[1,])
# dist_to_known$strand[1,closest_match1]
# 
# library(ggmotif) # Package on github at AliciaSchep/ggmotif. Can use seqLogo alternatively
# library(TFBSTools)
# # De novo motif
# ggmotif_plot(de_novos[[1]])
# # Closest matching known
# ggmotif_plot(toPWM(reverseComplement(motifs[[closest_match1]]),type = "prob"))



