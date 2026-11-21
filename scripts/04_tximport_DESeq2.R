setwd("/Volumes/RES/ws/paper1_ACLF_silence")
library(tximport)
library(DESeq2)

samples <- read.csv("data/processed/samplesheet_GSE228610_REAL.csv", stringsAsFactors = FALSE)

h <- readLines("ref/cdna_headers.txt")
tx <- sub(" .*$", "", sub("^>", "", h))
toks <- strsplit(h, " ")
get_field <- function(x, prefix) {
  i <- grep(paste0("^", prefix), x)
  if (length(i)) sub(paste0("^", prefix), "", x[i[1]]) else NA_character_
}
gene <- as.character(sapply(toks, get_field, prefix = "gene:"))
sym  <- as.character(sapply(toks, get_field, prefix = "gene_symbol:"))
tx2gene <- data.frame(tx = tx, gene = gene, symbol = sym, stringsAsFactors = FALSE)

files <- file.path("data/processed/salmon", samples$srr, "quant.sf")
names(files) <- samples$srr
stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("tx", "gene")])
saveRDS(txi, "data/processed/txi_GSE228610_REAL.rds")

samples$group <- factor(samples$group, levels = c("TAA_sham", "TAA_CLP", "water_sham", "water_CLP"))
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ group)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("过滤后基因数:", nrow(dds), "\n")

dds <- DESeq(dds)
symbol_map <- tx2gene[!duplicated(tx2gene$gene), c("gene", "symbol")]
rownames(symbol_map) <- symbol_map$gene
annotate <- function(res) { res$symbol <- symbol_map[rownames(res), "symbol"]; res }

res_primary     <- annotate(results(dds, contrast = c("group", "TAA_CLP", "TAA_sham"), alpha = 0.05))
res_secondary   <- annotate(results(dds, contrast = c("group", "water_CLP", "water_sham"), alpha = 0.05))
res_interaction <- annotate(results(dds, contrast = c(0, 1, 0, -1), alpha = 0.05))

write.csv(as.data.frame(res_primary),     "results/DE_TAA_CLP_vs_TAA_sham_REAL.csv")
write.csv(as.data.frame(res_secondary),   "results/DE_water_CLP_vs_water_sham_REAL.csv")
write.csv(as.data.frame(res_interaction), "results/DE_interaction_background_REAL.csv")
saveRDS(dds, "data/processed/dds_GSE228610_REAL.rds")

cat("\n====== 关键基因 · 主对比 TAA+CLP vs TAA+sham ======\n")
key_genes <- c("Angpt1","Angpt2","Hgf","Cebpb","Alb","Hnf4a","Ttr",
               "Casp3","Casp1","Gsdmd","Ripk3","Mlkl","Tnf","Fas",
               "Cdh5","Pecam1","Col1a1","Ptprc")
for (g in key_genes) {
  hit <- which(res_primary$symbol == g)
  if (length(hit)) {
    r <- res_primary[hit[1], ]
    cat(sprintf("%-8s baseMean=%8.1f  log2FC=%6.2f  padj=%.4g\n",
                g, r[["baseMean"]], r[["log2FoldChange"]], r[["padj"]]))
  } else cat(sprintf("%-8s 未检出\n", g))
}
cat("\n结果已写入 results/ 目录\n")
