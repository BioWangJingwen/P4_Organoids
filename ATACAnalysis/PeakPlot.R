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

# data: P4_4ATAC/P4_4_4peakDiff/P4_4_4_5RXX_peakView/P4_4_4_5RXX_peakAnnoPos.rds
