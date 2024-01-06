library(ggsci)
library(RColorBrewer)

col.seq <- c(ATAC = "turquoise", RNA = "seagreen2", WES = "tan1")
col.source <- c(Normal = "springgreen4", Primary = "orangered", Metastatic = "firebrick4")

col.exp1 <- c("#3C87C8", "#FFFFFF", "#EB5D1D")
col.exp2 <- c("darkgreen", "#C7E5C9", "darkred")
col.exp3 <- c("#3C87C8", "#FFD800", "#EB5D1D")
col.exp4 <- c("#253d24", "#ffffff", "#f34718")
col.exp5 <- c("#71ABB6", "#ffffff", "#C75C64")

col.IHC <- c(TNBC = "dodgerblue", `Luminal A` = "brown1", `Luminal B` = "lightseagreen", 'HER2+' = "maroon")

# col.IHC = pal_aaas()(5)
# names(col.IHC) <- c('TNBC', 'LuminalA', 'LuminalB', 'HER2+', 'LuminalA')
col.PAM50 = pal_cosmic()(5)
names(col.PAM50) <- c('Her2', 'Normal-like', 'Basal', 'LumA', 'LumB')
col.AIMS = pal_d3()(5)
names(col.AIMS) <- c('Her2', 'Normal-like', 'Basal', 'LumA', 'LumB')
col.claudinLow = pal_flatui()(2)
names(col.claudinLow) <- c("Others", "Claudin-low")


col.status <- c("0" = "grey60", "1" = "#5BAF4A")
col.stage <- c(pal_material("indigo")(5), "grey60")
# "#E7EAF6FF" "#C5CAE9FF" "#9FA7D9FF" "#7985CBFF" "#5B6BBFFF"    "grey60"
names(col.stage) <- c("I", "IIA", "IIB", " III", "IV",  "Tis")
col.sub1 <- c(TNBC = "dodgerblue", `HR+, HER2-` = "brown1", `HR+, HER2+` = "lightseagreen", 'HR-, HER2+' = "maroon")

# col.deg <-
#   pal_npg()(5)
# names(col.deg) <- c("NvsPM", "PvsN", "PvsM", "MvsP", "MvsN")
# col.deg <- col.deg[c("NvsPM", "PvsN",  "PvsM", "MvsP")]
# names(col.deg) <- mapvalues(names(col.deg), c("PvsM", "MvsP"), c("MvsP_down", "MvsP_up"))

col.deg <-
  c("#71ABB6", "#F0B57D", "#D3E1AE", "#C75C64")
names(col.deg) <- c("NvsPM", "PvsN", "MvsP_down", "MvsP_up")

col.muType <- c(Organoid = "#0099FF", Tissue = "#666699")
# col.mu <- c(c('Frame_Shift_Del' = "#45EDCC",'Missense_Mutation' = "#3EBDEE",
#               'Nonsense_Mutation' = "#C0BEBE", 'Multi_Hit' = "#ED45A5", 
#               'Frame_Shift_Ins' = "#45ED98",
#               'In_Frame_Ins' = "#9FC85E", 'Splice_Site' = "#E5291E", 'In_Frame_Del' = "#EE8044"))
col.mu <- brewer.pal(n = 8, name = "Set3") # Pastel1
names(col.mu) <-
  c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 
    'Multi_Hit', 'Frame_Shift_Ins',
    'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del'
    )

col.muSam <-
  brewer.pal(6, "Set1")
names(col.muSam) <- c('P030T', 'P110T', 'P068T', 'P078T', 'P049M', 'P049T')

col.Epi <-
  c(S01 = "#E41A71",
    S02 = "#379DB8",
    S03 = "#5BAF4A",
    S04 = "#7B4EA3",
    S05 = "#FF7600",
    S06 = "#FFC800")

col.EnDB <-
  c(GOBP = '#f9bdb0', GOCC = '#aacebc', GOMF = '#ead295', KEGG = '#8db0c5', REACTOME = '#ffddbb')




