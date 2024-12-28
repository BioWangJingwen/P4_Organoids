rm(list = ls())

setwd("/chenfeilab/Avocado/P4_Organoids/")
wd <- "/chenfeilab/Avocado/P4_Organoids/"

library(DESeq2)
library(GenomicRanges)
library(magrittr)
library(rangr)
library(dplyr)
library(ChIPseeker)

out.dir <- "P4_9CutNTag/P4_9_4Signals/P4_9_4_3CrossSignal/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = T)
}
out.prefix <- "P4_9_4_3Cross_"

# read --------------------------------------------------------------------

samInfoATAC <- readRDS("P4_4ATAC/P4_4_1QC/P4_4_1_1_samInfoATAC.rds")

# Get cluster
aa <- load("P4_7Clustering/P4_7_2ATACClusters/P4_7_2_1_res.Primary_Metastatic_Clusters.RData")
aa
optK <- 4
results <- rcc3
Cluster <- paste0("Cluster", results[[optK]]$consensusClass)
names(Cluster) <- names(results[[optK]]$consensusClass)
sam.C4 <- names(Cluster)[which(Cluster == "Cluster4")]

# Get bed -----------------------------------------------------------------

fil.peak <- list(
  MIST1 = list.files("P4_9CutNTag/P4_9_2Motifs/P4_9_2_3MergePeak", full.names = T, pattern = "MIST1.*.merge.narrowPeak"),
  TFAP4 = list.files("P4_9CutNTag/P4_9_2Motifs/P4_9_2_3MergePeak", full.names = T, pattern = "TFAP4.*.merge.narrowPeak")
)

out.bed <- paste0(out.dir, out.prefix, "bed/")
if (!dir.exists(out.bed)) {
  dir.create(out.bed, recursive = T)
}

fil.bed <- list(
  MIST1 = c(),
  TFAP4 = c()
)
for (t in names(fil.peak)) {
  # t <- names(fil.peak)[1]
  for (f in fil.peak[[t]]) {
    # f <- fil.peak[[t]][1]
    peak <- readPeakFile(f)
    s <- gsub("P4_9_2_3_", "", gsub(".merge.narrowPeak", "", basename(f)))
    f.p <- paste0(out.bed, out.prefix, s, "_sort_peakAll.bed")
    fil.bed[[t]] <- c(fil.bed[[t]], f.p)
    if (!file.exists(f.p)) {
      tbl <- as.data.frame(peak)
      tbl <- subset(tbl, seqnames %in% paste0("chr", 1:22))
      # po <- mixedorder(as.character(tbl$seqnames), decreasing = F)
      fb <- paste0(out.bed, out.prefix, s, "_peakAll.bed")
      write.table(tbl[, c("seqnames", "start", "end", "strand", "width")], 
                  file = fb, quote = F, sep = "\t",
                  row.names = F, col.names = F
      )
      system(paste0("sort -k1,1 -k2,2n ", fb, " > ", f.p), wait = T)
    }
  }
}

# fil bed -------------------------------------------------------------------

# heat of each sample 
fil.bw.all <- list.files(
  "/chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_2Motifs/P4_9_2_3MergePeak", 
  recursive = T, full.names = T, pattern = "merge.CPM.bw"
)
sam.all <- gsub(".merge.CPM.bw", "", gsub(".*_", "", basename(fil.bw.all)))
fil.bw.all <- fil.bw.all[which(sam.all %in% sam.C4)]

run <- c()
sam.bed <- unique(gsub(".merge.CPM.bw", "", gsub(".*_", "", basename(fil.bw.all))))
for (s in sam.bed) {
  # TFAP4
  fil.bed.use <- paste0(
    "/chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_4Signals/P4_9_4_3CrossSignal/P4_9_4_3Cross_bed/P4_9_4_3Cross_MIST1_", s, "_sort_peakAll.bed"
  )
  fils.bw.use <- paste0(
    "/chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_2Motifs/P4_9_2_3MergePeak/P4_9_2_3_TFAP4_", s, ".merge.CPM.bw"
  )
  f.mat <- paste0(wd, out.dir, out.prefix, "TFAP4_", s, ".mat.gz")
  if (!file.exists(f.mat)) {
    print("TFAP4")
    ### plot heatmap
    cmd <- paste0(
      "computeMatrix reference-point -S ",
      fils.bw.use, 
      " -R ",
      fil.bed.use,
      " --skipZeros -o ",
      f.mat, 
      " --samplesLabel ",
      "TFAP4", 
      " -p 8 -a 2000 -b 2000 --referencePoint center --missingDataAsZero --quiet"
    )
    run <- c(run, cmd)
  }
  f.plot <- paste0(wd, out.dir, out.prefix, "TFAP4_", s, ".pdf")
  if (!file.exists(f.plot)) {
    cmd <- paste0(
      "plotHeatmap -m ",
      f.mat, " -out ", f.plot,
      " --regionsLabel ",
      "MIST1_peak", 
      " --plotTitle ", s,
      " --colorList ", " '#ffffff, #374e55' ", " --legendLocation none --heatmapWidth 6 &"
    )
    run <- c(run, cmd)
  }
  #  --heatmapHeight 12 --heatmapWidth 3
  
  # MIST1
  fil.bed.use <- paste0(
    "/chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_4Signals/P4_9_4_3CrossSignal/P4_9_4_3Cross_bed/P4_9_4_3Cross_TFAP4_", s, "_sort_peakAll.bed"
  )
  fils.bw.use <- paste0(
    "/chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_2Motifs/P4_9_2_3MergePeak/P4_9_2_3_MIST1_", s, ".merge.CPM.bw"
  )
  f.mat <- paste0(wd, out.dir, out.prefix, "MIST1_", s, ".mat.gz")
  if (!file.exists(f.mat)) {
    print("MIST1")
    ### plot heatmap
    cmd <- paste0(
      "computeMatrix reference-point -S ",
      fils.bw.use, 
      " -R ",
      fil.bed.use,
      " --skipZeros -o ",
      f.mat, 
      " --samplesLabel ",
      "MIST1", 
      " -p 8 -a 2000 -b 2000 --referencePoint center --missingDataAsZero --quiet"
    )
    run <- c(run, cmd)
  }
  f.plot <- paste0(wd, out.dir, out.prefix, "MIST1_", s, ".pdf")
  if (!file.exists(f.plot)) {
    
    cmd <- paste0(
      "plotHeatmap -m ",
      f.mat, " -out ", f.plot,
      " --regionsLabel ",
      "TFAP4_peak", 
      " --plotTitle ", s,
      " --colorList ", " '#ffffff, #00a1d5' ", " --legendLocation none --heatmapWidth 6 &"
    )
    run <- c(run, cmd)
  }
  
}

write.table(run, file = paste0(out.dir, out.prefix, "run.sh"), quote = F, row.names = F, col.names = F)
run
length(run)


#  --zMin -3 --zMax 3






















