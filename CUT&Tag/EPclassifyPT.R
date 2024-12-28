rm(list = ls())

setwd("/chenfeilab/Avocado/P4_Organoids/")
wd <- "/chenfeilab/Avocado/P4_Organoids/"

# library(DESeq2)
library(GenomicRanges)
library(ChIPseeker)
library(ChIPpeakAnno)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ucsc.hg19.knownGene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(annotatr)

out.dir <- "P4_9CutNTag/P4_9_5EP/P4_9_5_1EPclassify/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = T)
}
# out.prefix <- "P4_9_4_1less_"
out.prefix <- "P4_9_5_1_"

# peakTool <- "/work/chenfeilab/Avocado/Software/peak-tool/human_hg19/peak_tool_multi"
# conda activate EP

# class -------------------------------------------------------------------

fil.bw.all <- list.files(
  "/chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_2Motifs/P4_9_2_3MergePeak", 
  recursive = T, full.names = T, pattern = "merge.CPM.bw"
)

fil.peak.all <- list.files(
  "/chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_2Motifs/P4_9_2_3MergePeak", 
  full.names = T, pattern = "merge.narrowPeak"
)
fil.peak.all <- grep("H3K27me3", fil.peak.all, invert = T, value = T)
fil.peak.all <- grep("IgG", fil.peak.all, invert = T, value = T)
fil.peak.all <- grep("EBF", fil.peak.all, invert = T, value = T)


# setwd("/work/chenfeilab/Avocado/Software/peak-tool/human_hg19")
run <- c()
for (f in fil.peak.all) {
  # f <- fil.peak.all[56]
  sam <- gsub(".merge.narrowPeak", "", gsub("P4_9_2_3_*", "", basename(f)))
  
  peak <- readPeakFile(f)
  
  # ann with ChIPSeeker
  f.cs <- paste0(out.dir, out.prefix, "ChIPSeeker", sam, ".rds")
  if (!file.exists(f.cs)) {
    ann.cs <- annotatePeak(peak, TxDb=txdb,
                           tssRegion=c(-1000, 1000),
                           verbose=FALSE,
                           addFlankGeneInfo=TRUE,
                           annoDb="org.Hs.eg.db"
    )
    saveRDS(ann.cs, file = f.cs)
    write.table(as.data.frame(ann.cs)[, c(1:6, 13, 21, 23)], file = paste0(out.dir, out.prefix, "ChIPSeeker", sam, ".txt"),
                sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(ann.cs@annoStat, file = paste0(out.dir, out.prefix, "ChIPSeekerTbl", sam, ".txt"),
                sep = "\t", col.names = T, row.names = F, quote = F)
  }
  
  # ann with ChIPpeakAnno
  f.ca <- paste0(out.dir, out.prefix, "ChIPpeakAnno", sam, ".rds")
  if (!file.exists(f.ca)) {
    ann.ca <- annotatePeakInBatch(peak, 
                                  AnnotationData = ucsc.hg19.knownGene)
    ann.ca <- addGeneIDs(annotatedPeak=ann.ca, 
                         orgAnn="org.Hs.eg.db", 
                         feature_id_type="entrez_id",
                         IDs2Add="symbol")
    
    saveRDS(ann.ca, file = f.ca)
    write.table(as.data.frame(ann.ca), file = paste0(out.dir, out.prefix, "ChIPpeakAnno", sam, ".txt"),
                sep = "\t", col.names = T, row.names = F, quote = F)
   
  }
  
  # ann with annotatr
  
  # ann with homer
  # annotatePeaks.pl /chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_2Motifs/P4_9_2_3MergePeak/P4_9_2_3_ZBP89_P084T.merge.narrowPeak hg19 > peak.annotation.xls
  f.hom <- paste0(out.dir, out.prefix, "homer", sam, ".xls")
  if (!file.exists(f.hom)) {
    print(sam)
    
    cmd <- paste0(
      "annotatePeaks.pl ",
      f, 
      " hg19 > ",
      wd, f.hom
    )
    system(cmd, wait = F)
  }
  
  # ann with peakTool
  f.bed <- paste0(out.dir, out.prefix, sam, ".bed")
  system(paste0(
    "cat ",
    f,
    " > ",
    f.bed
  ))
  
  
  f.mat <- paste0(out.dir, out.prefix, "peakTool", sam, ".txt")

  if (!file.exists(f.mat)) {
    print(sam)
    
    run <- c(run, 
             paste0(
               "./peak_tool_multi ",
               wd, f.bed, 
               " > ",
               wd, f.mat, "&"
             ))
    
  }
}

write.table(run, file = "/work/chenfeilab/Avocado/Softwares/peak-tool/human_hg19/run.sh", quote = F, row.names = F, col.names = F)
run
