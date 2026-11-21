suppressPackageStartupMessages(library(DESeq2))

setwd("/Volumes/RES/ws/paper1_ACLF_silence")
dds <- readRDS("data/processed/dds_GSE228610_REAL.rds")

symbol_map <- read.csv("results/DE_TAA_CLP_vs_TAA_sham_REAL.csv", stringsAsFactors = FALSE)
ens2sym <- symbol_map$symbol
names(ens2sym) <- symbol_map$X

# 正确的交互对比：交互 = (TAA_CLP - TAA_sham) - (water_CLP - water_sham) = -beta3
res_int <- results(dds, contrast = c(0, 0, 0, -1), alpha = 0.05)
res_int$symbol <- ens2sym[rownames(res_int)]
write.csv(as.data.frame(res_int), "results/DE_interaction_corrected_REAL.csv")

cat("====== 关键基因的真实交互效应 ======\n")
for (g in c("Alb","Ripk3","Angpt1","Angpt2","Hgf","Hnf4a","Ttr")) {
  hit <- which(res_int$symbol == g)
  if (length(hit)) {
    r <- res_int[hit[1], ]
    cat(sprintf("%-8s log2FC=%6.2f  padj=%.4g\n", g,
                r$log2FoldChange, r$padj))
  } else cat(sprintf("%-8s 未检出\n", g))
}

sig <- sum(res_int$padj < 0.05, na.rm = TRUE)
cat("\n交互对比 padj < 0.05 的基因总数:", sig, "\n\n")

cat("====== Top 15 交互基因 ======\n")
top <- as.data.frame(res_int[order(res_int$padj), ])
print(head(top[, c("symbol", "log2FoldChange", "padj")], 15))
