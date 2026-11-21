# 设置工作目录
setwd("/Volumes/RES/ws/paper1_ACLF_silence")

# 1. 读取数据 (注意：这里去掉 row.names=1，让数据作为普通 dataframe 读入)
data <- read.delim("data/raw/GSE228610_Raw_gene_counts_matrix.txt.gz", check.names = FALSE)

# 2. 把第 4 列 (Gene Symbol) 设为行名
# make.unique 是为了防止有重复的基因名报错
rownames(data) <- make.unique(data[, 4])

# 3. 提取目标列 (TAA_CLP: 2个, TAA_sham: 3个)
target_cols <- c("CLP_03_02 FPKM", "CLP_03_04 FPKM", "CLP_03_07 FPKM", "CLP_03_09 FPKM", "CLP_03_11 FPKM")
expr <- data[, target_cols]

# 4. 定义分组
group <- factor(c(rep("TAA_CLP", 2), rep("TAA_sham", 3)), levels = c("TAA_sham", "TAA_CLP"))

# 5. 数据转换 (FPKM 转 log2)
expr_log2 <- log2(expr + 1)

# 6. 过滤极低表达的基因（保证在至少一半样本中表达量 > 0）
keep <- rowSums(expr_log2 > 0) >= 3
expr_log2 <- expr_log2[keep, ]

# 7. 小样本必要步骤：分位数标准化
library(limma)
expr_norm <- normalizeBetweenArrays(expr_log2, method = "quantile")

# 8. 线性模型拟合
design <- model.matrix(~ group)
fit <- lmFit(expr_norm, design)
fit <- eBayes(fit)

# 9. 提取所有结果
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
write.csv(res, "results/DE_TAA_CLP_vs_TAA_sham.csv")

# 10. 检查关键靶点
cat("\n====== 探索性差异分析完成 ======\n")
cat("结果已保存到 results/DE_TAA_CLP_vs_TAA_sham.csv\n")

cat("\n====== 关键基因效应方向 (LogFC, t, B) ======\n")
check_genes <- c("Angpt1", "Angpt2", "Hgf", "Alb", "Hnf4a", "Ttr", "Cebpb")
for (g in check_genes) {
  if (g %in% rownames(res)) {
    cat(sprintf("%-8s LogFC: %6.2f, t: %5.2f, B: %5.2f\n", g, res[g, "logFC"], res[g, "t"], res[g, "B"]))
  } else {
    cat(sprintf("%-8s 未在表达矩阵中检出\n", g))
  }
}
