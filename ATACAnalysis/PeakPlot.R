# P4_4ATAC/P4_4_4peakDiff/P4_4_4_5RXX_peakView/P4_4_4_5RXX_VennPeaks
# data: P4_4ATAC/P4_4_4peakDiff/P4_4_4_5RXX_peakView/P4_4_4_5RXX_SigPeak.rds

source("/chenfeilab/Avocado/MyPipline/peakVenn.R")

pdf(paste0(out.dir, out.prefix, "VennPeaks.pdf"), width = 4.5, height = 4.5)
myVenn(mcols(allPeaksSet_nR), 
       circle.col = col.deg[colnames(mcols(allPeaksSet_nR))],
       alpha = 0.7
)
dev.off()

# P4_4ATAC/P4_4_4peakDiff/P4_4_4_5RXX_peakView/P4_4_4_5RXX_PositionPeaks
# data: P4_4ATAC/P4_4_4peakDiff/P4_4_4_5RXX_peakView/P4_4_4_5RXX_peakAnnoPos.rds
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

pdf(paste0(out.dir, out.prefix, "PositionPeaks.pdf"), width = 7, height = 2.5)
peakPlot
dev.off()


# P4_4ATAC/P4_4_4peakDiff/P4_4_4_5RXX_peakView/P4_4_4_5RXX_HeatDP
ann.row <- rowAnnotation(
  NvsPM = dat.row$NvsPM,
  PvsN = dat.row$PvsN,
  MvsP = dat.row$MvsP,
  col = list(NvsPM = c(Up = unname(col.deg["NvsPM"]), Non_DCAP = "white"),
             PvsN = c(Up = unname(col.deg["PvsN"]), Non_DCAP = "white"),
             MvsP = c(Up = unname(col.deg["MvsP_up"]), Down = unname(col.deg["MvsP_down"]),  Non_DCAP = "white")
  )
)

po <- match(colnames(mat.exp), samInfoATAC$SampleID)

ann.col <- 
  columnAnnotation(
    Source = anno_block(gp = gpar(fill = col.source)),
    # PAM50 = samInfoATAC$PAM50[po],
    # AIMS = samInfoATAC$AIMS[po],
    # claudinLow = samInfoATAC$claudinLow[po],
    `menopausal status` = samInfoATAC$`menopausal status`[po],
    Stage = samInfoATAC$Stage[po],
    `Molecular subtype1` = samInfoATAC$`Molecular subtype1`[po],
    col = list(IHC = col.IHC, 
               PAM50 = col.PAM50,
               AIMS = col.AIMS,
               claudinLow = col.claudinLow,
               `menopausal status` = col.status,
               Stage = col.stage,
               `Molecular subtype1` = col.sub1
    )
  )

# spl.row <- 
#   factor(c(rep("Normal", length(gen.nor)), rep("Primary", length(gen.pri)), rep("Metastatic", length(gen.met))),
#          levels = c("Normal", "Primary", "Metastatic"))
spl.col <- factor(samInfoATAC$Source[po], levels = c("Normal", "Primary", "Metastatic"))

set.seed(123456)
h.DP <-
  Heatmap(mat.exp,
          col = circlize::colorRamp2(breaks = c(0, 5, 10), colors = col.exp5),
          top_annotation = ann.col,
          left_annotation = ann.row,
          name = "Scaled_Exp",
          # row_title_rot = 0,
          show_row_names = F, 
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

Cairo::CairoPNG(paste0(out.dir, out.prefix, "HeatDP.png"), width = 9, height = 9, units = "in", res = 108)
h.DP
dev.off()
