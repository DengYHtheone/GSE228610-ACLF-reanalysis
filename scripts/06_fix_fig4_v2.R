suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(ggplot2)
})

setwd("/Volumes/RES/ws/paper1_ACLF_silence")
dir.create("figures", showWarnings = FALSE)

dds <- readRDS("data/processed/dds_GSE228610_REAL.rds")
res <- results(dds, contrast = c("group", "TAA_CLP", "TAA_sham"))

symbol_map <- read.csv("results/DE_TAA_CLP_vs_TAA_sham_REAL.csv", stringsAsFactors = FALSE)
ens2sym <- symbol_map$symbol
names(ens2sym) <- symbol_map$X

gene_list <- res$log2FoldChange
names(gene_list) <- ens2sym[rownames(res)]
gene_list <- gene_list[!is.na(names(gene_list)) & !is.na(gene_list)]
gene_list <- tapply(gene_list, names(gene_list), mean)
gene_list <- sort(gene_list, decreasing = TRUE)

pathways <- list(
  LiverMetabolism = c("Alb","Ttr","Hnf4a","G6pc","Pck1","Cyp3a11",
                      "Cyp2e1","Sult2a1","Aox1","Cat","Tat","Arg1"),
  CellDeath = c("Ripk3","Mlkl","Casp3","Casp1","Gsdmd","Tnf",
                "Fas","Bax","Bak1","Cyba","Nox2")
)

set.seed(42)
fgsea_res <- fgsea(pathways, gene_list, minSize = 3, maxSize = 500)
print(fgsea_res[, .(pathway, NES, pval, padj, size)])

fgsea_res$label <- paste0("NES = ", round(fgsea_res$NES, 2),
                          "  (p = ", format.pval(fgsea_res$pval, digits = 2), ")")
fgsea_res$pathway <- factor(fgsea_res$pathway,
                            levels = fgsea_res$pathway[order(fgsea_res$NES)])
fgsea_res$signif <- ifelse(fgsea_res$padj < 0.05, "FDR < 0.05", "ns")

# 标签放条形外侧，黑色，不会溢出
fgsea_res$label_y <- fgsea_res$NES + ifelse(fgsea_res$NES > 0, 0.15, -0.15)
fgsea_res$label_h <- ifelse(fgsea_res$NES > 0, 0, 1)

p_nes <- ggplot(fgsea_res, aes(x = pathway, y = NES, fill = signif)) +
  geom_col(width = 0.55) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(aes(y = label_y, label = label, hjust = label_h), size = 3.8) +
  scale_fill_manual(values = c("FDR < 0.05" = "#d62728", "ns" = "#7f7f7f")) +
  scale_y_continuous(limits = c(-2.4, 2.4)) +
  coord_flip() +
  theme_classic() +
  labs(title = "GSEA: TAA+CLP vs TAA+sham",
       x = "", y = "Normalized Enrichment Score (NES)") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

ggsave("figures/Fig4_GSEA_NES_v2.png", p_nes, width = 7.5, height = 3.5, dpi = 200)
ggsave("figures/Fig4_GSEA_NES_v2.pdf", p_nes, width = 7.5, height = 3.5)

cat("✅ Fig4 v2 已生成\n")
