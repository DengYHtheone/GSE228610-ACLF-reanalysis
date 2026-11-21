setwd("/Volumes/RES/ws/paper1_ACLF_silence")
library(DESeq2)

counts <- read.delim("data/raw/GSE311736_counts.txt.gz", row.names = 1, check.names = FALSE)

# 1. 看看 Hgf 到底在哪
cat("\n====== Hgf 的原始 Counts 数据 ======\n")
hgF_id <- "ENSMUSG00000028853" # Hgf 的 ID
if(hgF_id %in% rownames(counts)) {
  print(counts[hgF_id, c(paste0("LSEC_CLP-", 1:4), paste0("LSEC_Sham-", 1:4))])
} else {
  cat("矩阵中完全找不到这个 ID\n")
}

# 2. 简单的死亡基因集评分 (基于肝细胞)
cat("\n====== 肝细胞死亡/炎症基因集得分 (使用简单 Z-score) ======\n")
hep_cols <- c(paste0("HEP_CLP-", 1:4), paste0("HEP_Sham-", 1:4))
hep_counts <- counts[, hep_cols]

# 选一组经典的凋亡和坏死性凋亡基因 (以及炎症基因做对照)
# 这里用真实存在的 Ensembl ID 替代
death_genes <- c(
  "ENSMUSG00000028163" = "Casp3", 
  "ENSMUSG00000025888" = "Casp1", 
  "ENSMUSG00000022534" = "Gsdmd", 
  "ENSMUSG00000027500" = "Stat1", 
  "ENSMUSG00000025907" = "Tnf", 
  "ENSMUSG00000030149" = "Fas" 
)
# 将 counts 转换为 logCPM
cpm_mat <- log2(t(t(hep_counts) / colSums(hep_counts) * 1e6) + 1)

cat("基因名\tSham均值\tCLP均值\t变化方向\n")
for (id in names(death_genes)) {
  if (id %in% rownames(cpm_mat)) {
    sham_mean <- mean(as.numeric(cpm_mat[id, 1:4]))
    clp_mean <- mean(as.numeric(cpm_mat[id, 5:8]))
    dir <- ifelse(clp_mean > sham_mean, "上调", "下调")
    cat(sprintf("%-8s %.2f\t%.2f\t%s\n", death_genes[id], sham_mean, clp_mean, dir))
  }
}
cat("提示：如果是单纯CLP，大量坏死基因如果不显著上调，则支持‘不是死亡’\n")

