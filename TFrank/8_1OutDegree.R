rm(list = ls())

library(igraph)

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
out.dir <- "P4_8TFNet/P4_8_2TFRank/P4_8_2_1OutDegree/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.prefix <- "P4_8_2_1_"


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

# out degree --------------------------------------------------------------

fil.deg <- paste0(out.dir, out.prefix, "mat.deg.rds")
if (!file.exists(fil.deg)) {
  tfs <- unique(reg.net$TF)
  sams <- unique(reg.net$SampleID)
  mat.deg <- 
    matrix(0, nrow = length(tfs), ncol = length(sams), 
           dimnames = list(tfs, sams))
  
  out.graph <- paste0(out.dir, out.prefix, "Graphs/")
  if (!dir.exists(out.graph)) {
    dir.create(out.graph)
  }
  
  for (s in unique(reg.net$SampleID)) {
    print(s)
    # s <- "P001M"
    reg.net.s <- subset(reg.net, SampleID == s)
    # dim(reg.net.s)
    # head(reg.net.s)
    g <- graph_from_data_frame(reg.net.s[, c("TF", "Gene")], directed = TRUE)
    g <- simplify(g)
    saveRDS(g, file = paste0(out.graph, out.prefix, s, "_Graph.rds"))
    
    aa <- get.edgelist(g)
    deg <- 
      degree(
        g,
        v = unique(aa[, 1]),
        mode = "out",
        loops = F,
        normalized = FALSE
      )
    mat.deg[names(deg), s] <- deg
  }
  
  saveRDS(mat.deg, file = fil.deg)
} else {
  mat.deg <- readRDS(fil.deg)
}

# O_each sample -----------------------------------------------------------

out.deg <- paste0(out.dir, out.prefix, "OutdegreeEach/")
if (!dir.exists(out.deg)) {
  dir.create(out.deg)
}
# M  N  T 
# 39 67 75 

for (index in c("M",  "N",  "T")) {
  pdf(paste0(out.deg, out.prefix, index, "Outdegree.pdf"), width = 18, height = 21)
  par(mfrow = c(10, 8))
  for (s in grep(index, colnames(mat.deg), value = T)) {
    # s <- sam.ord[i]
    plot(density(mat.deg[, s], n = 30), 
         type = "b", xlab = "Number of target genes", 
         col = "royalblue", main = s)
  }
  while(!is.null(dev.list())) dev.off()
}


# O_diff NPM ------------------------------------------------------------------

samN <- grep("N$", colnames(mat.deg), value = T)
samP <- grep("T$", colnames(mat.deg), value = T)
samM <- grep("M$", colnames(mat.deg), value = T)

O.diff.NPM.lis <- list()

### NvsPM
out.sum <- rowSums(mat.deg)
out.N <- rowSums(log2(mat.deg[, samN]/out.sum + 10^(-24)))/length(samN)
out.PM <- rowSums(log2(mat.deg[, c(samP, samM)]/out.sum + 10^(-24)))/length(c(samP, samM))
O.NvsPM <- out.N - out.PM
summary(O.NvsPM); sort(O.NvsPM, decreasing = T)[1:10]
O.diff.NPM.lis[["NvsPM"]] <- O.NvsPM
write.table(O.NvsPM, file = paste0(out.dir, out.prefix, "O.diffNvsPM.txt"), quote = F, row.names = T, col.names = F)

### PvsN
out.sumPN <- rowSums(mat.deg[, c(samN, samP)])
out.N <- rowSums(log2(mat.deg[, samN]/out.sumPN + 10^(-24)))/length(samN)
out.P <- rowSums(log2(mat.deg[, samP]/out.sumPN + 10^(-24)))/length(samP)
O.PvsN <- out.P - out.N
summary(O.PvsN); sort(O.PvsN, decreasing = T)[1:10]
O.diff.NPM.lis[["PvsN"]] <- O.PvsN
write.table(O.PvsN, file = paste0(out.dir, out.prefix, "O.diffPvsN.txt"), quote = F, row.names = T, col.names = F)

### MvsP
out.sumPM <- rowSums(mat.deg[, c(samM, samP)])
out.M <- rowSums(log2(mat.deg[, samM]/out.sumPM + 10^(-24)))/length(samM)
out.P <- rowSums(log2(mat.deg[, samP]/out.sumPM + 10^(-24)))/length(samP)
O.MvsP <- out.M - out.P
summary(O.MvsP); sort(O.MvsP, decreasing = T)[1:10]
O.diff.NPM.lis[["MvsP"]] <- O.MvsP
write.table(O.MvsP, file = paste0(out.dir, out.prefix, "O.diffMvsP.txt"), quote = F, row.names = T, col.names = F)

saveRDS(O.diff.NPM.lis, file = paste0(out.dir, out.prefix, "O.diff.NPM.lis.rds"))

# O_diff Clusters ---------------------------------------------------------

out.sum <- rowSums(mat.deg)
O.diff.clu.lis <- list()
for (clu in unique(Cluster)) {
  # clu <- "Cluster4"
  sam1 <- names(Cluster)[which(Cluster == clu)]
  sam2 <- names(Cluster)[which(Cluster != clu)]
  out1 <- rowSums(log2(mat.deg[, sam1]/out.sum + 10^(-24)))/length(sam1)
  out2 <- rowSums(log2(mat.deg[, sam2]/out.sum + 10^(-24)))/length(sam2)
  
  O.diff <- out1 - out2
  O.diff.clu.lis[[clu]] <- O.diff
  
  write.table(O.diff, file = paste0(out.dir, out.prefix, "O.diff", clu, ".txt"), quote = F, row.names = T, col.names = F)
  
}

saveRDS(O.diff.clu.lis, file = paste0(out.dir, out.prefix, "O.diff.clu.lis.rds"))













