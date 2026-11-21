suppressPackageStartupMessages(library(ggplot2))

setwd("/Volumes/RES/ws/paper1_ACLF_silence")

box <- function(x, y, w, h, label, fill = "#EAF2F8", border = "#2C3E50",
                fontsize = 3.4, fontface = "plain") {
  list(
    geom_rect(aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
              fill = fill, color = border, linewidth = 0.7),
    geom_text(aes(x = x, y = y, label = label), size = fontsize,
              fontface = fontface, lineheight = 0.95)
  )
}

arr <- function(x1, y1, x2, y2) {
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.24, "cm"), type = "closed"),
               linewidth = 0.8, color = "#2C3E50")
}

group_cols <- c("water+sham" = "#D6EAF8", "water+CLP" = "#FDEBD0",
                "TAA+sham" = "#D5F5E3", "TAA+CLP" = "#FADBD8")

p <- ggplot() +
  geom_rect(aes(xmin = 8, xmax = 44, ymin = 46, ymax = 60),
            fill = NA, color = "#5D6D7E", linetype = "dashed", linewidth = 0.7) +
  geom_rect(aes(xmin = 46, xmax = 82, ymin = 46, ymax = 60),
            fill = NA, color = "#5D6D7E", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 26, y = 63, label = "Water background",
           size = 3.8, fontface = "bold") +
  annotate("text", x = 64, y = 63, label = "TAA background (cirrhosis)",
           size = 3.8, fontface = "bold") +
  box(17, 53, 15, 9, "water + sham\n(n = 3)", fill = group_cols["water+sham"]) +
  box(35, 53, 15, 9, "water + CLP\n(n = 3)", fill = group_cols["water+CLP"]) +
  box(55, 53, 15, 9, "TAA + sham\n(n = 3)", fill = group_cols["TAA+sham"]) +
  box(73, 53, 15, 9, "TAA + CLP\n(n = 3)", fill = group_cols["TAA+CLP"]) +
  arr(26, 46, 26, 41.5) +
  arr(64, 46, 64, 41.5) +
  box(45, 38, 62, 7,
      "Whole-liver RNA-seq, 24 h post-CLP\nGSE228610 - 12 samples - this study: reanalysis",
      fill = "#E8DAEF", fontsize = 3.6) +
  arr(45, 34.5, 45, 29) +
  box(10, 23, 17, 10, "Salmon 2.11.0\ntranscript\nquantification\n(Ensembl GRCm39)",
      fill = "#EBF5FB", fontsize = 3.2) +
  box(30, 23, 17, 10, "tximport 1.38.2\ngene-level\nsummarization",
      fill = "#EBF5FB", fontsize = 3.2) +
  box(50, 23, 19, 10, "DESeq2 1.50.2\ndifferential expression\n(~ group; primary,\nsecondary, interaction)",
      fill = "#FDEBD0", fontsize = 3.1) +
  box(70, 23, 17, 10, "fgsea 1.36.2\ngene set\nenrichment\n(CellDeath /\nLiverMetabolism)",
      fill = "#FDEBD0", fontsize = 3.1) +
  box(90, 23, 15, 10, "Figures\n2-4, S1\nAlb - Ripk3 -\nGSEA - PCA",
      fill = "#D5F5E3", fontsize = 3.2) +
  arr(18.5, 23, 21.5, 23) +
  arr(38.5, 23, 40.5, 23) +
  arr(59.5, 23, 61.5, 23) +
  arr(78.5, 23, 82.5, 23) +
  annotate("text", x = 50, y = 12, label = "Study design and analysis workflow",
           size = 4.6, fontface = "bold") +
  annotate("text", x = 50, y = 8,
           label = "CLP: cecal ligation and puncture; TAA: thioacetamide; GSEA: gene set enrichment analysis",
           size = 2.9, color = "#626567") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 66)) +
  theme_void()

ggsave("figures/Fig1_design.png", p, width = 10, height = 6.6, dpi = 300, bg = "white")
ggsave("figures/Fig1_design.pdf", p, width = 10, height = 6.6, bg = "white")

cat("DONE Fig1_design v2\n")
