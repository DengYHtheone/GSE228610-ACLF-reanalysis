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

res <- read.csv("results/DE_interaction_final_REAL.csv", stringsAsFactors = FALSE)
sig <- subset(res, padj < 0.05)
sig <- sig[order(-abs(sig$log2FoldChange)), ]

detox <- c("Cyp2c37","Cyp2c23","Sult2a8","Akr1c19","Hsd3b3","Aox3",
           "Afmid","Mup10","Slc22a30","Bdh2","Tlcd2")

mat38 <- mat[rownames(mat) %in% sig$X, ]
z <- t(scale(t(mat38)))
z[z > 2] <- 2; z[z < -2] <- -2

df <- expand.grid(gene = rownames(z), srr = colnames(z))
df$z <- as.vector(z)
df$symbol <- sig$symbol[match(df$gene, sig$X)]
df$group <- samples$group_label[match(df$srr, samples$srr)]
df$module <- ifelse(df$symbol %in% detox, "drug/xenobiotic", "other")
df$symbol <- factor(df$symbol, levels = rev(sig$symbol))
df$srr <- factor(df$srr, levels = samples$srr[order(samples$group_label)])

p <- ggplot(df, aes(x = srr, y = symbol, fill = z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(module ~ group, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "row z-score") +
  theme_minimal() +
  labs(title = "38 genes with significant background x sepsis interaction",
       x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        strip.text.y = element_text(angle = 0, size = 7),
        panel.grid = element_blank())

ggsave("figures/Fig_heatmap38.png", p, width = 8, height = 9, dpi = 300, bg = "white")
ggsave("figures/Fig_heatmap38.pdf", p, width = 8, height = 9, bg = "white")
cat("DONE heatmap38\n")
