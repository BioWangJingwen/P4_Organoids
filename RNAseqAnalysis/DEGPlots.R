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

# P4_2RNA/P4_2_3DEG1v1/P4_2_3_2RXX_vennDEGALL2
ggvenn_4 <-
  ggplot(dat.row, aes(A = `NvsPM`, B = `PvsN`, C = `MvsP_down`, D = `MvsP_up`)) +
  geom_venn(fill_color=col.deg[names(lisDEG.PosALL)],fill_alpha = .7,set_name_size = 8,
            set_name_color = col.deg[names(lisDEG.PosALL)],
            text_size=5) + 
  theme_void()+
  coord_fixed() 
  
  # labs(title = "Example of ggvenn:: geom_venn function",
  #      subtitle = "processed charts with geom_venn()",
  #      caption = "Visualization by DataCharm") +
  # theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black",face = 'bold',
  #                                 size = 20, margin = margin(t = 1, b = 12)),
  #       plot.subtitle = element_text(hjust = 0,vjust = .5,size=15),
  #       plot.caption = element_text(face = 'bold',size = 12))
pdf(paste0(out.dir, out.prefix, "vennDEGALL2.pdf"), width = 9, height = 9)
ggvenn_4
dev.off()
