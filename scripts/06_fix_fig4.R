suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(ggplot2)
})

setwd("/Volumes/RES/ws/paper1_ACLF_silence")

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
fgsea_res$pathway <- factor(fgsea_res$
