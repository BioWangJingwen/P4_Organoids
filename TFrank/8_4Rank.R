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
library(ggsci)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

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

out.prefix <- "P4_8_2_4_"

source("P4_8TFNet/P4_8_2_4RankFUN.R")


# legand ------------------------------------------------------------------

pdf(paste0(out.dir, out.prefix, "Rank_Legends.pdf"), width = 2, height = 4.5)
# Chromatin\naccessibility-g Outdegree-r Expression-b
grid.newpage()
col_fun = colorRamp2(seq(1, 30, length.out = 9), pal_material("red")(9))
lgd = Legend(col_fun = col_fun, title = "Outdegree", labels = c("Low", "High"), at = c(0, 30))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.8, "npc"), just = c("left", "center"))

col_fun = colorRamp2(seq(1, 30, length.out = 9), pal_material("green")(9))
lgd = Legend(col_fun = col_fun, title = "Chromatin\naccessibility", labels = c("Low", "High"), at = c(0, 30))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("left", "center"))

col_fun = colorRamp2(seq(1, 30, length.out = 9), pal_material("blue")(9))
lgd = Legend(col_fun = col_fun, title = "Expression", labels = c("Low", "High"), at = c(0, 30))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.2, "npc"), just = c("left", "center"))
while(!is.null(dev.list())) dev.off()

# Rank NPM ----------------------------------------------------------------

O_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_1OutDegree/P4_8_2_1_O.diff.NPM.lis.rds")
A_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_2Accessibility/P4_8_2_2_A.diff.NPM.lis.rds")
E_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_3Expression/P4_8_2_3_Ediff.lis_NPM.rds")

names(O_diff)
names(A_diff)
names(E_diff)

O_gen <- unique(unlist(lapply(O_diff, names)))
A_gen <- unique(unlist(lapply(A_diff, names)))
E_gen <- unique(unlist(lapply(E_diff, names)))

# genes fixed
gen <- intersect(O_gen, intersect(A_gen, E_gen))

O_lis <- lapply(O_diff, function(x)x[gen])
O_mat <- do.call(rbind, O_lis)

A_lis <- lapply(A_diff, function(x)x[gen])
A_mat <- do.call(rbind, A_lis)

E_lis <- lapply(E_diff, function(x)x[gen])
E_mat <- do.call(rbind, E_lis)

# order matrix
O_ord_mat <- apply(O_mat, 1, rank)
# rownames(O_ord_mat) <- colnames(O_mat)

A_ord_mat <- apply(A_mat, 1, rank)
# rownames(A_ord_mat) <- colnames(A_mat)

E_ord_mat <- apply(E_mat, 1, rank)
# rownames(E_ord_mat) <- colnames(E_mat)

clu <- sort(colnames(O_ord_mat))
ord_mat <- O_ord_mat[gen, clu] + A_ord_mat[gen, clu] + E_ord_mat[gen, clu]
write.table(ord_mat, file = paste0(out.dir, out.prefix, "ord_mat_NPM.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

# NvsPM
s <- "NvsPM"
title <- "NvsPM"

tfR <- sort(ord_mat[, s], decreasing = T)[1:25]
dat <- data.frame(
  TF = names(tfR),
  TF_rank = tfR,
  O_diff = O_mat[s, names(tfR)],
  A_diff = A_mat[s, names(tfR)],
  E_diff = E_mat[s, names(tfR)],
  Type = "y"
)
dat$TF <- factor(dat$TF, levels = dat$TF)

g.N <-
  GetRankPlot(dat, title)

# PvsN
s <- "PvsN"
title <- "PvsN"

tfR <- sort(ord_mat[, s], decreasing = T)[1:25]
dat <- data.frame(
  TF = names(tfR),
  TF_rank = tfR,
  O_diff = O_mat[s, names(tfR)],
  A_diff = A_mat[s, names(tfR)],
  E_diff = E_mat[s, names(tfR)],
  Type = "y"
)
dat$TF <- factor(dat$TF, levels = dat$TF)

g.P <-
  GetRankPlot(dat, title)


# MvsP_Pos
s <- "MvsP"
title <- "MvsP_Pos"

E_breaks <- c(min(E_mat[s, ]), max(E_mat[s, ]))
O_breaks <- c(min(O_mat[s, ]), max(O_mat[s, ]))
A_breaks <- c(min(A_mat[s, ]), max(A_mat[s, ]))

tfR <- sort(ord_mat[, s], decreasing = T)[1:25]
dat <- data.frame(
  TF = names(tfR),
  TF_rank = tfR,
  O_diff = O_mat[s, names(tfR)],
  A_diff = A_mat[s, names(tfR)],
  E_diff = E_mat[s, names(tfR)],
  Type = "y"
)
dat$TF <- factor(dat$TF, levels = dat$TF)

g.M_Pos <-
  GetRankPlot(dat, title) # , E_breaks = E_breaks, O_breaks = O_breaks, A_breaks = A_breaks


# MvsP_Neg
s <- "MvsP"
title <- "MvsP_Neg"

tfR <- sort(ord_mat[, s], decreasing = F)[1:25]
dat <- data.frame(
  TF = names(tfR),
  TF_rank = tfR,
  O_diff = O_mat[s, names(tfR)],
  A_diff = A_mat[s, names(tfR)],
  E_diff = E_mat[s, names(tfR)],
  Type = "y"
)
dat$TF <- factor(dat$TF, levels = dat$TF)

g.M_Neg <-
  GetRankPlot(dat, title) #, E_breaks = E_breaks, O_breaks = O_breaks, A_breaks = A_breaks

pdf(paste0(out.dir, out.prefix, "Rank_NPM.pdf"), width = 9, height = 7)
ggarrange(g.N, g.P, g.M_Neg, g.M_Pos, 
          ncol = 2, nrow = 2,
          align = "hv"
)
while(!is.null(dev.list())) dev.off()


# Rank Cluster ------------------------------------------------------------

O_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_1OutDegree/P4_8_2_1_O.diff.clu.lis.rds")
A_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_2Accessibility/P4_8_2_2_A.diff.clu.lis.rds")
E_diff <- readRDS("P4_8TFNet/P4_8_2TFRank/P4_8_2_3Expression/P4_8_2_3_Ediff.lis_clu.rds")

names(O_diff)
names(A_diff)
names(E_diff)

O_gen <- unique(unlist(lapply(O_diff, names)))
A_gen <- unique(unlist(lapply(A_diff, names)))
E_gen <- unique(unlist(lapply(E_diff, names)))

# genes fixed
gen <- intersect(O_gen, intersect(A_gen, E_gen))

O_lis <- lapply(O_diff, function(x)x[gen])
O_mat <- do.call(rbind, O_lis)

A_lis <- lapply(A_diff, function(x)x[gen])
A_mat <- do.call(rbind, A_lis)

E_lis <- lapply(E_diff, function(x)x[gen])
E_mat <- do.call(rbind, E_lis)

# order matrix
O_ord_mat <- apply(O_mat, 1, rank)
rownames(O_ord_mat) <- colnames(O_mat)

A_ord_mat <- apply(A_mat, 1, rank)
rownames(A_ord_mat) <- colnames(A_mat)

E_ord_mat <- apply(E_mat, 1, rank)
rownames(E_ord_mat) <- colnames(E_mat)

clu <- sort(colnames(O_ord_mat))
ord_mat <- O_ord_mat[gen, clu] + A_ord_mat[gen, clu] + E_ord_mat[gen, clu]
head(ord_mat)

write.table(ord_mat, file = paste0(out.dir, out.prefix, "ord_mat_Clusters.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

# Clusters
g.lis <- list()
dat.top <- data.frame(row.names = 1:25)
for (clu in paste0("Cluster", 1:4)) {
  print(clu)
  tfR <- sort(ord_mat[, clu], decreasing = T)[1:25]
  dat <- data.frame(
    TF = names(tfR),
    TF_rank = tfR,
    O_diff = O_mat[clu, names(tfR)],
    A_diff = A_mat[clu, names(tfR)],
    E_diff = E_mat[clu, names(tfR)],
    Type = "y"
  )
  dat.top[, clu] <- dat$TF
  dat$TF <- factor(dat$TF, levels = dat$TF)
   
  g.lis[[clu]] <-
    GetRankPlot(dat, clu)
}

saveRDS(dat.top, file = paste0(out.dir, out.prefix, "Rank_TF.rds"))
pdf(paste0(out.dir, out.prefix, "Rank_Clu.pdf"), width = 9, height = 7)
ggarrange(plotlist = g.lis, 
          ncol = 2, nrow = 2,
          align = "hv"
)
dev.off()

g <- c("MAZ", "ZNF148")
c <- "Cluster4"
O_mat[c, g]
E_mat[c, g]
A_mat[c, g]

O_ord_mat[g, c]
E_ord_mat[g, c]
A_ord_mat[g, c]

ord_mat[g, c]
