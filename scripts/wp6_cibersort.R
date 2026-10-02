#!/usr/bin/env Rscript
# =====================================================================
# WP6a  CIBERSORT 去卷积（LM22 官方脚本 + 签名，TCGA-LUAD 502 例）
# 方案: 最终执行方案_二区落地_2026-10-01.md WP6（修订 v2）
# 输入: pathB/TCGA-LUAD_HiSeqV2.gz (Xena, log2) 中取 WP2 的 502 例
# 输出: WP6_cibersort_fractions.csv
# =====================================================================
.libPaths(c("/Users/liubao/.workbuddy/Rlibs", .libPaths()))
suppressMessages({library(data.table); library(tibble)})

DIR  <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP6_immune"
EXPR <- "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/TCGA-LUAD_HiSeqV2.gz"
LM22 <- "/Users/liubao/WorkBuddy/2026-09-23-10-23-05/生信学习计划/第8课数据/LM22.txt"
CSRF <- "/Users/liubao/WorkBuddy/2026-09-23-10-23-05/生信学习计划/第8课数据/CIBERSORT.R"

## ---- WP2 患者（12位ID + 分组）----
sam <- read.delim("/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP2_subgroup_survival/WP2_分析用样本表.tsv",
                  stringsAsFactors = FALSE)
stopifnot(nrow(sam) == 502)

## ---- 读表达（data.table 按列选 502 个肿瘤样本）----
E <- fread(EXPR)
genes <- E[[1]]
col12 <- substr(colnames(E)[-1], 1, 12)
keep <- which(col12 %in% sam$patient & grepl("-01$", colnames(E)[-1])) + 1L
cat("匹配肿瘤样本列:", length(keep), "\n")
M <- as.matrix(E[, ..keep])   # as.matrix 保留列名，勿再重赋值（偏移坑）
rownames(M) <- genes
rm(E); gc()
# 去重基因（取均值）
dup <- duplicated(rownames(M))
if(any(dup)){
  gg <- unique(rownames(M)[dup])
  M2 <- sapply(gg, function(g) colMeans(M[rownames(M) == g, , drop = FALSE]))
  M <- rbind(M[!dup, , drop = FALSE], t(M2))
}
cat("表达矩阵:", nrow(M), "基因 x", ncol(M), "样本\n")
stopifnot(all(substr(colnames(M), 1, 12) %in% sam$patient))

## ---- 运行 CIBERSORT ----
source(CSRF)
sig <- read.delim(LM22, row.names = 1, check.names = FALSE)
sig <- as.matrix(sig)
cat("LM22 签名:", nrow(sig), "基因 x", ncol(sig), "细胞类型\n")
cat("交集基因:", length(intersect(rownames(M), rownames(sig))), "\n")

t0 <- Sys.time()
res <- CIBERSORT(sig_matrix = sig, mixture_file = as.data.frame(M),
                 perm = 0, QN = TRUE, parallel = FALSE)
cat("CIBERSORT 耗时:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

out <- as.data.frame(res)
out$sample <- rownames(out)
out$patient <- substr(out$sample, 1, 12)
out$grp <- sam$grp[match(out$patient, sam$patient)]
write.csv(out, file.path(DIR, "WP6_cibersort_fractions.csv"), row.names = FALSE)
cat("saved WP6_cibersort_fractions.csv  rows:", nrow(out), "\n")
