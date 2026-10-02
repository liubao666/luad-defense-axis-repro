#!/usr/bin/env Rscript
# =====================================================================
# WP6b  MCP-counter + 免疫检查点 + CYT（TCGA-LUAD 502 例）
# MCP-counter 口径: 每细胞群体 = 其 marker 基因 log2 表达的算术平均
#                   (= 线性空间几何平均, 即发表方法的几何平均口径)
# CYT 口径: Rooney et al. 细胞毒性活性 = mean(GZMA, PRF1) log2
# =====================================================================
.libPaths(c("/Users/liubao/.workbuddy/Rlibs", .libPaths()))
suppressMessages({library(data.table)})

DIR  <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP6_immune"
EXPR <- "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/TCGA-LUAD_HiSeqV2.gz"
MCPG <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP6_immune/MCPcounter-master/Signatures/genes.txt"

sam <- read.delim("/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP2_subgroup_survival/WP2_分析用样本表.tsv",
                  stringsAsFactors = FALSE)
stopifnot(nrow(sam) == 502)

E <- fread(EXPR)
genes <- E[[1]]
allc <- colnames(E)[-1]
keep <- which(substr(allc, 1, 12) %in% sam$patient & grepl("-01$", allc)) + 1L
M <- as.matrix(E[, ..keep]); rownames(M) <- genes
rm(E); gc()
dup <- duplicated(rownames(M))
if(any(dup)){
  gg <- unique(rownames(M)[dup])
  M2 <- sapply(gg, function(g) colMeans(M[rownames(M) == g, , drop = FALSE]))
  M <- rbind(M[!dup, , drop = FALSE], t(M2))
}
cat("表达矩阵:", nrow(M), "基因 x", ncol(M), "样本\n")

## ---- MCP-counter ----
g <- read.delim(MCPG, stringsAsFactors = FALSE)
colnames(g) <- c("gene", "cell")
pops <- unique(g$cell)
mcp <- sapply(pops, function(cp){
  mk <- g$gene[g$cell == cp]
  mk <- intersect(mk, rownames(M))
  if(length(mk) < 1) return(rep(NA_real_, ncol(M)))
  colMeans(M[mk, , drop = FALSE])
})
mcp <- as.data.frame(mcp)
mcp$patient <- substr(colnames(M), 1, 12)
mcp$grp <- sam$grp[match(mcp$patient, sam$patient)]
cat("MCP-counter 群体数:", length(pops), "| 每群体 marker 数:",
    paste(sapply(pops, function(cp) sum(g$gene[g$cell==cp] %in% rownames(M))), collapse = ","), "\n")
write.csv(mcp, file.path(DIR, "WP6_mcpcounter_scores.csv"), row.names = FALSE)

## ---- 检查点 + CYT ----
cp_genes <- c("CD274","PDCD1","CTLA4","LAG3","TIGIT")
missing <- setdiff(c(cp_genes, "GZMA","PRF1"), rownames(M))
if(length(missing)) cat("缺失基因:", missing, "\n")
ck <- data.frame(patient = substr(colnames(M), 1, 12))
for(gene in cp_genes) ck[[gene]] <- M[gene, ]
ck$CYT <- colMeans(M[c("GZMA","PRF1"), , drop = FALSE])
ck$grp <- sam$grp[match(ck$patient, sam$patient)]
write.csv(ck, file.path(DIR, "WP6_checkpoint_cyt.csv"), row.names = FALSE)
cat("saved WP6_mcpcounter_scores.csv / WP6_checkpoint_cyt.csv\n")
