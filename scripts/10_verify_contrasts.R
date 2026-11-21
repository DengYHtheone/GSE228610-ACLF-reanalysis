suppressPackageStartupMessages(library(DESeq2))
setwd("/Volumes/RES/ws/paper1_ACLF_silence")
dds <- readRDS("data/processed/dds_GSE228610_REAL.rds")

sm <- read.csv("results/DE_TAA_CLP_vs_TAA_sham_REAL.csv", stringsAsFactors = FALSE)
ens2sym <- sm$symbol; names(ens2sym) <- sm$X

res_pri <- results(dds, contrast = c("group", "TAA_CLP", "TAA_sham"))
res_sec <- results(dds, contrast = c("group", "water_CLP", "water_sham"))
res_int <- results(dds, contrast = c(0, 0, 0, -1))

cat(sprintf("%-8s %10s %10s %10s %12s\n", "gene", "primary", "secondary", "interaction", "pri-sec=?int"))
for (g in c("Alb", "Ripk3", "Angpt1")) {
  e <- names(ens2sym)[ens2sym == g][1]
  a <- res_pri[e, "log2FoldChange"]
  b <- res_sec[e, "log2FoldChange"]
  c <- res_int[e, "log2FoldChange"]
  cat(sprintf("%-8s %10.3f %10.3f %10.3f %12.3f\n", g, a, b, c, a - b))
}
