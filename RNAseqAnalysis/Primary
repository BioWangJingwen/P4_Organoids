# P4_2RNA/P4_2_9Primary/P4_2_9_1Scores/P4_2_9_1_gbox.gsea


gg.gsea <-
  ggplot(mel.gsea, aes(pN_Stage, Score)) +
  facet_wrap(.~Signature, ncol = 5, nrow = 2, scales = "free") +
  # geom_boxplot(aes(fill = pN_Stage)) +
    geom_violin(aes(fill = pN_Stage), draw_quantiles = 0.5) +
  # geom_point() +
  stat_compare_means(comparisons = list(c("N0", "N3")),
                     label = "p.signif", vjust = 1.2,
                     method = "t.test"
                     ) +
  scale_fill_manual(values = alpha(col.pN, alpha = 0.38)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1, size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid = element_blank()
  ) +
  labs(x = "", y = "ssGSEA scores") +
  NoLegend()

gg.gsea

Cairo::CairoPDF(paste0(out.dir, out.prefix, "gbox.gsea.pdf"), width = 8, height = 6)
gg.gsea
dev.off()
