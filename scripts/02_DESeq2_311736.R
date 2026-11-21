setwd("/Volumes/RES/ws/paper1_ACLF_silence")
library(DESeq2)

cat("正在读取数据并准备分析（跳过注释包，直接使用 Ensembl ID）...\n")
counts <- read.delim("data/raw/GSE311736_counts.txt.gz", row.names = 1, check.names = FALSE)

# 设置我们关心的基因对应的 Ensembl ID
target_ensembl <- list(
  "Alb" = "ENSMUSG00000029368", 
  "Ttr" = "ENSMUSG00000037022", 
  "Hnf4a" = "ENSMUSG00000017950", 
  "Cebpb" = "ENSMUSG00000056501", 
  "Casp1" = "ENSMUSG00000025888", 
  "Gsdmd" = "ENSMUSG00000022534",
  "Angpt1" = "ENSMUSG00000029911", 
  "Angpt2" = "ENSMUSG00000031465", 
  "Hgf" = "ENSMUSG00000028853", 
  "Pecam1" = "ENSMUSG00000020717", 
  "Cdh5" = "ENSMUSG00000031821"
)

# ---- 分析模块 1：肝细胞 (HEP) ----
cat("\n====== 正在分析肝细胞 (HEP) ======\n")
hep_cols <- c(paste0("HEP_CLP-", 1:4), paste0("HEP_Sham-", 1:4))
hep_counts <- counts[, hep_cols]
hep_group <- factor(c(rep("CLP", 4), rep("Sham", 4)), levels = c("Sham", "CLP"))

dds_hep <- DESeqDataSetFromMatrix(countData = hep_counts, colData = data.frame(hep_group), design = ~ hep_group)
dds_hep <- DESeq(dds_hep)
res_hep <- results(dds_hep)
write.csv(as.data.frame(res_hep), "results/DESeq2_HEP.csv")

cat("肝细胞关键基因 (Log2FC, P-value):\n")
hep_genes <- c("Alb", "Ttr", "Hnf4a", "Cebpb", "Casp1", "Gsdmd")
for (g in hep_genes) {
  ens_id <- target_ensembl[[g]]
  if (ens_id %in% rownames(res_hep)) {
    cat(sprintf("%-8s Log2FC: %6.2f, P: %f\n", g, res_hep[ens_id, "log2FoldChange"], res_hep[ens_id, "pvalue"]))
  } else { cat(sprintf("%-8s 未检出\n", g)) }
}

# ---- 分析模块 2：内皮细胞 (LSEC) ----
cat("\n====== 正在分析内皮细胞 (LSEC) ======\n")
lsec_cols <- c(paste0("LSEC_CLP-", 1:4), paste0("LSEC_Sham-", 1:4))
lsec_counts <- counts[, lsec_cols]
lsec_group <- factor(c(rep("CLP", 4), rep("Sham", 4)), levels = c("Sham", "CLP"))

dds_lsec <- DESeqDataSetFromMatrix(countData = lsec_counts, colData = data.frame(lsec_group), design = ~ lsec_group)
dds_lsec <- DESeq(dds_lsec)
res_lsec <- results(dds_lsec)
write.csv(as.data.frame(res_lsec), "results/DESeq2_LSEC.csv")

cat("内皮细胞关键基因 (Log2FC, P-value):\n")
lsec_genes <- c("Angpt1", "Angpt2", "Hgf", "Pecam1", "Cdh5")
for (g in lsec_genes) {
  ens_id <- target_ensembl[[g]]
  if (ens_id %in% rownames(res_lsec)) {
    cat(sprintf("%-8s Log2FC: %6.2f, P: %f\n", g, res_lsec[ens_id, "log2FoldChange"], res_lsec[ens_id, "pvalue"]))
  } else { cat(sprintf("%-8s 未检出\n", g)) }
}
cat("\n====== 分析完成，结果已保存 ======\n")
