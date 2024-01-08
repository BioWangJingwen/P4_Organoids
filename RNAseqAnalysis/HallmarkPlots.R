# P4_2RNA/P4_2_5PathAct/P4_2_5_1RXX_HeatGSEA

po <- match(rownames(mat.gsea), dat.row$Hallmark)
col.NvsPM <- colorRampPalette(c("white", col.deg["NvsPM"]))(5)
col.PvsN <- colorRampPalette(c("white", col.deg["PvsN"]))(5)
col.MvsP <- colorRampPalette(c("white", col.deg["MvsP_up"]))(5)
names(col.NvsPM) <- names(col.PvsN) <- names(col.MvsP) <- names(col.sig)

ann.lef <-
  rowAnnotation(
    NvsPM = dat.row$NvsPM[po],
    PvsN = dat.row$PvsN[po],
    MvsP = dat.row$MvsP[po],
    col = list(NvsPM = col.NvsPM,
               PvsN = col.PvsN,
               MvsP = col.MvsP
               )
  )

po <- match(colnames(mat.gsea), samInfoOrgSel$RepID)

ann.col <- 
  columnAnnotation(
    Source = anno_block(gp = gpar(fill = col.source))
  )

spl.col <- factor(samInfoOrgSel$Source[po], levels = c("Normal", "Primary", "Metastatic"))

col.rowname <- unlist(lapply(p.adj, function(p){
  if (p < 0.01) {
    c = "red"
  } else if (p < 0.05) {
    c = "orange"
  } else {
    c = 'black'
  }
  return(c)
}))

set.seed(123456)
h.deg <-
  Heatmap(mat.gsea,
          col = colorRamp2(breaks = c(-5, 0, 5), colors = col.exp1),
          top_annotation = ann.col,
          # left_annotation = ann.row, # used to check colors
          left_annotation = ann.lef, 
          name = "ssGSEA scores",
          # row_title_rot = 0,
          # show_row_names = F, 
          row_names_gp = gpar(col = col.rowname),
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
h.deg

Cairo::CairoPDF(paste0(out.dir, out.prefix, "HeatGSEA.pdf"), width = 12, height = 9)
h.deg
dev.off()

# P4_2RNA/P4_2_5PathAct/P4_2_5_5GeneProgram/P4_2_5_5RXX_HeatNMF

po <- match(rownames(mat.gsea), dat.row$Hallmark)
col.NvsPM <- colorRampPalette(c("white", col.deg["NvsPM"]))(5)
col.PvsN <- colorRampPalette(c("white", col.deg["PvsN"]))(5)
col.MvsP <- colorRampPalette(c("white", col.deg["MvsP_up"]))(5)
names(col.NvsPM) <- names(col.PvsN) <- names(col.MvsP) <- names(col.sig)

ann.lef <-
  rowAnnotation(
    NvsPM = dat.row$NvsPM[po],
    PvsN = dat.row$PvsN[po],
    MvsP = dat.row$MvsP[po],
    col = list(NvsPM = col.NvsPM,
               PvsN = col.PvsN,
               MvsP = col.MvsP
    )
  )

po <- match(colnames(mat.gsea), samInfoOrgSel$RepID)

ann.col <- 
  columnAnnotation(
    Source = anno_block(gp = gpar(fill = col.source))
  )

spl.col <- factor(samInfoOrgSel$Source[po], levels = c("Normal", "Primary", "Metastatic"))

col.rowname <- unlist(lapply(p.adj, function(p){
  if (p < 0.01) {
    c = "red"
  } else if (p < 0.05) {
    c = "orange"
  } else {
    c = 'black'
  }
  return(c)
}))

set.seed(123456)
h.deg <-
  Heatmap(mat.gsea,
          col = colorRamp2(breaks = c(-5, 0, 5), colors = col.exp1),
          top_annotation = ann.col,
          # left_annotation = ann.row, # used to check colors
          left_annotation = ann.lef, 
          name = "ssGSEA scores",
          # row_title_rot = 0,
          # show_row_names = F, 
          row_names_gp = gpar(col = col.rowname),
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
h.deg

# Cairo::CairoPNG(paste0(out.dir, out.prefix, "HeatGSEA.png"), width = 12, height = 9, units = "in", res = 108)
# h.deg
# dev.off()

Cairo::CairoPDF(paste0(out.dir, out.prefix, "HeatNMF.pdf"), width = 12, height = 9)
h.deg
dev.off()












