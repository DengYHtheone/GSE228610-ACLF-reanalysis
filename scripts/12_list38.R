suppressPackageStartupMessages(library(DESeq2))
setwd("/Volumes/RES/ws/paper1_ACLF_silence")
res <- read.csv("results/DE_interaction_final_REAL.csv", stringsAsFactors = FALSE)
sig <- subset(res, padj < 0.05)
sig <- sig[order(-abs(sig$log2FoldChange)), ]
cat("显著交互基因数:", nrow(sig), "\n\n")
print(sig[, c("symbol", "log2FoldChange", "padj")], row.names = FALSE)
