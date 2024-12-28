rm(list = ls())

setwd("/chenfeilab/Avocado/P4_Organoids/")
wd <- "/chenfeilab/Avocado/P4_Organoids/"

out.dir <- "P4_8TFNet/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.dir <- "P4_8TFNet/P4_8_1GetNet/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.pat <- "P4_8TFNet/P4_8_1GetNet/P4_8_1_7RegNetSample/"
if (!dir.exists(out.pat)) {
  dir.create(out.pat)
}
out.prefix <- "P4_8_1_7_"


# read data ---------------------------------------------------------------

reg.net <- readRDS("/chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_1GetNet/P4_8_1_6_reg.net.rds")
reg.net <- na.omit(reg.net)

# Get regulation network --------------------------------------------------

reg.net$fdr <- p.adjust(reg.net$p, method = "fdr")
reg.net.sel1 <- subset(reg.net, r > 0.6 & fdr < 0.01)

saveRDS(reg.net.sel1, file = paste0(out.dir, out.prefix, "reg.netR0.6FDR0.01.rds"))

reg.net.sel2 <- subset(reg.net, r > 0 & fdr < 0.01)

saveRDS(reg.net.sel2, file = paste0(out.dir, out.prefix, "reg.netR0FDR0.01.rds"))

# write reg net -----------------------------------------------------------

head(reg.net.sel1)
reg.net.sel1 <- as.data.frame(reg.net.sel1)
for (p in unique(reg.net.sel1$SampleID)) {
  # p <- "P073M"
  print(p)
  reg.p <- reg.net.sel1[which(reg.net.sel1$SampleID == p), ]
  print(dim(reg.p))
  write.table(reg.p, file = paste0(out.pat, out.prefix, p, "_reg.netR0.6FDR0.01.txt"),
              sep = "\t", row.names = F, col.names = T, quote = F)
}

reg.p <- reg.net[which(reg.net$SampleID == p), ]



