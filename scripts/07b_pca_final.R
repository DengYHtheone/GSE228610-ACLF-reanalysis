suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

setwd("/Volumes/RES/ws/paper1_ACLF_silence")

dds <- readRDS("data/processed/dds_GSE228610_REAL.rds")
samples <- read.csv("data/processed/samplesheet_GSE228610_REAL.csv", stringsAsFactors = FALSE)

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
samples <- samples[match(colnames(mat), samples$srr), ]
samples$group_label <- gsub("_", "+", samples$group)
samples$group_label <- factor(samples$group_label,
                              levels = c("water+sham", "water+CLP", "TAA+sham", "TAA+CLP"))

pca <- prcomp(t(mat))
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cat("各主成分解释方差 (%):", paste(var_exp[1:4], collapse = ", "), "\n")

df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                 group = samples$group_label, srr = samples$srr)

cols <- c("water+sham" = "#1f77b4", "water+CLP" = "#ff7f0e",
          "TAA+sham" = "#2ca02c", "TAA+CLP" = "#d62728")

p <- ggplot(df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(aes(label = srr), size = 3, show.legend = FALSE,
                           max.overlaps = 20) +
  scale_color_manual(values = cols) +
  theme_classic() +
  labs(title = "PCA of VST-transformed gene expression",
       x = paste0("PC1 (", var_exp[1], "% variance)"),
       y = paste0("PC2 (", var_exp[2], "% variance)"),
       color = "")

ggsave("figures/FigS1_PCA.png", p, width = 7, height = 5.5, dpi = 200)
ggsave("figures/FigS1_PCA.pdf", p, width = 7, height = 5.5)

cat("✅ FigS1_PCA 最终版已生成\n")
