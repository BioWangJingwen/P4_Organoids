rm(list = ls())
setwd("/chenfeilab/Avocado/P4_Organoids/")
library(ggplot2)
library(ggthemes)
library(DESeq2)
library(plyr)

out.dir <- "P4_2RNA/P4_2_3DEG1v1/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_2_3_1_"

# load data
source("P4_1Data/P4_1_0Colors/P4_1_0colors.R")
samInfoOrgSel <- readRDS("P4_1Data/P4_1_1StaInfo/P4_1_1_2_2samInfoOrgSel.rds")
samInfoOrgSel <- subset(samInfoOrgSel, Sample.type == "RNA-Seq")
dim(samInfoOrgSel)

cnts <- readRDS("P4_2RNA/P4_2_2DEG/P4_2_2_1_counts.used.rds")

table(samInfoOrgSel$Source)

# Get N2PM ----------------------------------------------------------------

resN2PM <- readRDS("P4_2RNA/P4_2_2DEG/P4_2_2_1_NormalRes.rds")
print(head(resN2PM))

pos <- subset(na.omit(resN2PM), log2FoldChange > 1 & padj < 0.01)
dim(pos)
gN2PM <- rownames(pos)

# Get P2M -----------------------------------------------------------------

samInfoP2M <- subset(samInfoOrgSel, Source %in% c("Primary", "Metastatic"))
cntsP2M <- cnts[, samInfoP2M$RepID]

cond <- samInfoP2M[match(colnames(cntsP2M), samInfoP2M$RepID), "Source"]
# cond[which(cond != gp1)] <- "Others"
cond <- factor(cond, levels = c("Metastatic", "Primary"))

# object construction
dds <- DESeqDataSetFromMatrix(cntsP2M, DataFrame(cond), ~ cond)

# standard analysis
dds <- DESeq(dds)
resP2M <- results(dds)
print(head(resP2M))
saveRDS(resP2M, file = paste0(out.dir, out.prefix, "ResP2M.rds"))

pos <- subset(na.omit(resP2M), log2FoldChange > 1 & padj < 0.01)
dim(pos)
gP2M <- rownames(pos)

pos <- subset(na.omit(resP2M), log2FoldChange < -1 & padj < 0.01)
dim(pos)
gM2P <- rownames(pos)


# Get M2P -----------------------------------------------------------------

samInfoP2M <- subset(samInfoOrgSel, Source %in% c("Primary", "Metastatic"))
cntsP2M <- cnts[, samInfoP2M$RepID]

cond <- samInfoP2M[match(colnames(cntsP2M), samInfoP2M$RepID), "Source"]
# cond[which(cond != gp1)] <- "Others"
cond <- factor(cond, levels = c("Primary", "Metastatic"))

# object construction
dds <- DESeqDataSetFromMatrix(cntsP2M, DataFrame(cond), ~ cond)

# standard analysis
dds <- DESeq(dds)
resP2M <- results(dds)
print(head(resP2M))
saveRDS(resP2M, file = paste0(out.dir, out.prefix, "ResM2P.rds"))


# Get P2N -----------------------------------------------------------------
samInfoP2N <- subset(samInfoOrgSel, Source %in% c("Primary", "Normal"))
cntsP2N <- cnts[, samInfoP2N$RepID]

cond <- samInfoP2N[match(colnames(cntsP2N), samInfoP2N$RepID), "Source"]
# cond[which(cond != gp1)] <- "Others"
cond <- factor(cond, levels = c("Normal", "Primary"))

# object construction
dds <- DESeqDataSetFromMatrix(cntsP2N, DataFrame(cond), ~ cond)

# standard analysis
dds <- DESeq(dds)
resP2N <- results(dds)
print(head(resP2N))
saveRDS(resP2N, file = paste0(out.dir, out.prefix, "ResP2N.rds"))

pos <- subset(na.omit(resP2N), log2FoldChange > 1 & padj < 0.01)
dim(pos)
gP2N <- rownames(pos)

# Get M2N -----------------------------------------------------------------

samInfoM2N <- subset(samInfoOrgSel, Source %in% c("Metastatic", "Normal"))
cntsM2N <- cnts[, samInfoM2N$RepID]

cond <- samInfoM2N[match(colnames(cntsM2N), samInfoM2N$RepID), "Source"]
# cond[which(cond != gp1)] <- "Others"
cond <- factor(cond, levels = c("Normal", "Metastatic"))

# object construction
dds <- DESeqDataSetFromMatrix(cntsM2N, DataFrame(cond), ~ cond)

# standard analysis
dds <- DESeq(dds)
resM2N <- results(dds)
print(head(resM2N))
saveRDS(resM2N, file = paste0(out.dir, out.prefix, "ResM2N.rds"))

pos <- subset(na.omit(resM2N), log2FoldChange > 1 & padj < 0.01)
dim(pos)
gM2N <- rownames(pos)

lisDEG.Pos <-
  list(Normal = gN2PM,
       Primary = union(gP2N, gP2M),
       Metastatic = union(gM2N, gM2P)
  )
lisDEG.PosALL <-
  list(NvsPM = gN2PM,
       PvsN = gP2N,
       PvsM = gP2M,
       MvsN = gM2N,
       MvsP = gM2P
  )

saveRDS(lisDEG.Pos, file = paste0(out.dir, out.prefix, "lisDEG.Pos.rds"))
saveRDS(lisDEG.PosALL, file = paste0(out.dir, out.prefix, "lisDEG.PosALL.rds"))
