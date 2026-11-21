suppressPackageStartupMessages(library(DESeq2))
setwd("/Volumes/RES/ws/paper1_ACLF_silence")
dds <- readRDS("data/processed/dds_GSE228610_REAL.rds")

rn <- resultsNames(dds)
cat("=== 系数真实顺序 ===\n"); print(rn)

coef_names <- setdiff(rn, "(Intercept)")
m <- regmatches(coef_names, regexec("^group_(.+)_vs_(.+)$", coef_names))
pairs <- do.call(rbind, lapply(m, function(x) x[2:3]))
ref <- names(sort(table(pairs[, 2]), decreasing = TRUE))[1]
cat("\n参照组:", ref, "\n")

v <- setNames(rep(0, length(rn)), rn)
v[paste0("group_TAA_CLP_vs_", ref)]    <- v[paste0("group_TAA_CLP_vs_", ref)] + 1
v[paste0("group_water_sham_vs_", ref)] <- v[paste0("group_water_sham_vs_", ref)] + 1
v[paste0("group_water_CLP_vs_", ref)]  <- v[paste0("group_water_CLP_vs_", ref)] - 1
cat("\n=== 交互对比向量 ===\n"); print(v)

res_int <- results(dds, contrast = as.numeric(v), alpha = 0.05)

sm <- read.csv("results/DE_TAA_CLP_vs_TAA_sham_REAL.csv", stringsAsFactors = FALSE)
ens2sym <- sm$symbol; names(ens2sym) <- sm$X
res_int$symbol <- ens2sym[rownames(res_int)]
write.csv(as.data.frame(res_int), "results/DE_interaction_final_REAL.csv")

res_pri <- results(dds, contrast = c("group", "TAA_CLP", "TAA_sham"))
res_sec <- results(dds, contrast = c("group", "water_CLP", "water_sham"))

cat("\n=== 终检表（最后一列必须等于第三列）===\n")
cat(sprintf("%-8s %10s %10s %10s %12s\n", "gene", "primary", "secondary", "interaction", "pri-sec"))
for (g in c("Alb", "Ripk3", "Angpt1", "Angpt2", "Hgf")) {
  e <- names(ens2sym)[ens2sym == g][1]
  a <- res_pri[e, "log2FoldChange"]; b <- res_sec[e, "log2FoldChange"]
  cc <- res_int[e, "log2FoldChange"]; p <- res_int[e, "padj"]
  cat(sprintf("%-8s %10.3f %10.3f %10.3f %12.3f  padj=%.4g\n", g, a, b, cc, a - b, p))
}

cat("\n交互 padj<0.05 基因数:", sum(res_int$padj < 0.05, na.rm = TRUE), "\n")
