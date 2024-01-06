# P4_2RNA/P4_2_3DEG1v1/P4_2_3_2RXX_vennDEGALL
vennALL <-
  ggVennDiagram(lisDEG.PosALL, 
                set_color = col.deg[names(lisDEG.PosALL)]
  ) +
  scale_fill_material("light-blue") +
  scale_color_manual(values = col.deg[names(lisDEG.PosALL)]) +
  NoLegend()
vennALL

save(lisDEG.PosALL, col.deg, file = paste0(out.dir, out.prefix, "vennDEGALL.RData"))
pdf(paste0(out.dir, out.prefix, "vennDEGALL.pdf"), width = 9, height = 9)
vennALL
dev.off()
