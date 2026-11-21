suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggpubr)
})

setwd("/Volumes/RES/ws/paper1_ACLF_silence")
dir.create("figures", showWarnings = FALSE)

dds <- readRDS("data/processed/dds_GSE228610_REAL.rds")
samples <- read.csv("data/processed/samplesheet_GSE228610_REAL.csv", stringsAsFactors = FALSE)

vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)
samples <- samples[match(colnames(expr_matrix), samples$srr), ]

# 关键修复：直接从 group 列生成标签，不依赖可能不存在的 background 列
samples$group_label <- gsub("_", "+", samples$group)
samples$group_label <- factor(samples$group_label,
                              levels = c("water+sham", "water+CLP", "TAA+sham", "TAA+CLP"))

cat("=== 分组检查（应为 3-3-3-3）===\n")
print(table(samples$group_label, useNA = "ifany"))

# 基因查找：Ensembl -> symbol，多个转录本取表达最高的
symbol_map <- read.csv("results/DE_TAA_CLP_vs_TAA_sham_REAL.csv", stringsAsFactors = FALSE)
ens2sym <- symbol_map$symbol
names(ens2sym) <- symbol_map$X

get_expr <- function(symbol) {
  ens <- names(ens2sym)[ens2sym == symbol]
  ens <- intersect(ens, rownames(expr_matrix))
  if (length(ens) == 0) return(NULL)
  if (length(ens) > 1) {
    means <- rowMeans(expr_matrix[ens, , drop = FALSE])
    ens <- names(which.max(means))
  }
  as.numeric(expr_matrix[ens, ])
}

group_colors <- c("water+sham" = "#1f77b4", "water+CLP" = "#ff7f0e",
                  "TAA+sham" = "#2ca02c", "TAA+CLP" = "#d62728")

plot_gene <- function(symbol) {
  e <- get_expr(symbol)
  if (is.null(e)) { cat(symbol, "未找到\n"); return(NULL) }
  df <- data.frame(expression = e, group = samples$group_label)
  ggplot(df, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 3) +
    stat_compare_means(comparisons = list(c("water+sham", "water+CLP"),
                                          c("TAA+sham", "TAA+CLP"),
                                          c("water+CLP", "TAA+CLP")),
                       method = "t.test", label = "p.format") +
    scale_fill_manual(values = group_colors) +
    theme_classic() +
    labs(title = paste(symbol, "Expression"), x = "", y = "VST expression") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

p_alb <- plot_gene("Alb")
p_ripk3 <- plot_gene("Ripk3")

ggsave("figures/Fig2_Alb_expression.pdf", p_alb, width = 6, height = 4)
ggsave("figures/Fig3_Ripk3_expression.pdf", p_ripk3, width = 6, height = 4)
ggsave("figures/Fig2_Alb_expression.png", p_alb, width = 6, height = 4, dpi = 200)
ggsave("figures/Fig3_Ripk3_expression.png", p_ripk3, width = 6, height = 4, dpi = 200)

cat("✅ 图2、图3 已生成（PDF + PNG）\n")
