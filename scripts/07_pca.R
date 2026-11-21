suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))

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

df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = samples$group_label)
cat("各主成分解释方差 (%):", paste(var_exp[1:4], collapse = ", "), "\n")

cols <- c("water+sham" = "#1f77b4", "water+CLP" = "#ff7f0e",
          "TAA+sham" = "#2ca02c", "TAA+CLP" = "#d62728")

p <- ggplot(df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95, linewidth = 0.6) +
  scale_color_manual(values = cols) +
  theme_classic() +
  labs(title = "PCA of VST-transformed gene expression",
       x = paste0("PC1 (", var_exp[1], "% variance)"),
       y = paste0("PC2 (", var_exp[2], "% variance)"),
       color = "")

ggsave("figures/FigS1_PCA.png", p, width = 6.5, height = 5, dpi = 200)
ggsave("figures/FigS1_PCA.pdf", p, width = 6.5, height = 5)

cat("✅ FigS1_PCA 已生成\n")
