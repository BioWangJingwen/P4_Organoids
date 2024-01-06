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
pdf(paste0(out.dir, out.prefix, "vennDEGALL2.pdf"), width = 9, height = 9)
ggvenn_4
dev.off()

# P4_2RNA/P4_2_3DEG1v1/P4_2_3_2RXX_HeatDEGALL
ann.row <- rowAnnotation(
  NvsPM = dat.row$NvsPM,
  PvsN = dat.row$PvsN,
  MvsP = dat.row$MvsP,
  col = list(NvsPM = c(Up = unname(col.deg["NvsPM"]), Non_DEG = "white"),
             PvsN = c(Up = unname(col.deg["PvsN"]), Non_DEG = "white"),
             MvsP = c(Up = unname(col.deg["MvsP_up"]), Down = unname(col.deg["MvsP_down"]),  Non_DEG = "white")
             )
)
po <- match(colnames(mat.exp), samInfoOrgSel$RepID)
samInfoOrgSel$Stage <- factor(samInfoOrgSel$Stage, levels = c("I", "IIA", "IIB", " III", "IV",  "Tis"))
ann.col <- 
  columnAnnotation(
    Source = anno_block(gp = gpar(fill = col.source)),
    # IHC = samInfoOrgSel$`Molecular subtype2`[po],
    PAM50 = samInfoOrgSel$PAM50[po],
    AIMS = samInfoOrgSel$AIMS[po],
    claudinLow = samInfoOrgSel$claudinLow[po],
    `menopausal status` = samInfoOrgSel$`menopausal status`[po],
    Stage = samInfoOrgSel$Stage[po],
    `Molecular subtype1` = samInfoOrgSel$`Molecular subtype1`[po],
    col = list(IHC = col.IHC, 
               PAM50 = col.PAM50,
               AIMS = col.AIMS,
               claudinLow = col.claudinLow,
               `menopausal status` = col.status,
               Stage = col.stage,
               `Molecular subtype1` = col.sub1
    )
  )
spl.col <- factor(samInfoOrgSel$Source[po], levels = c("Normal", "Primary", "Metastatic"))

set.seed(123456)
h.deg <-
  Heatmap(mat.exp,
          col = colorRamp2(breaks = c(-5, 0, 5), colors = col.exp5),
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
Cairo::CairoPDF(paste0(out.dir, out.prefix, "HeatDEGALL.pdf"), width = 9, height = 7)
h.deg
dev.off()
