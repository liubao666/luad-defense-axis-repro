# 模块二（干净版）：剔线粒体基因后，全肿瘤细胞虚拟敲除 IGF2BP3
.libPaths(c("~/Rlibs", .libPaths()))
library(scTenifoldNet)
library(scTenifoldKnk)
library(Matrix)

pcNet_orig <- getFromNamespace("pcNet", "scTenifoldNet")
pcNet_fast <- function(X, nComp = 3, scaleScores = TRUE, symmetric = FALSE,
                       q = 0, priorNetwork = NULL, verbose = FALSE,
                       nCores = min(4, max(1, parallel::detectCores() - 2)),
                       useRcpp = FALSE) {
  pcNet_orig(X, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric,
             q = q, priorNetwork = priorNetwork, verbose = verbose,
             nCores = nCores, useRcpp = FALSE)
}
assignInNamespace("pcNet", pcNet_fast, "scTenifoldNet")

BASE <- "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
counts <- readMM(file.path(BASE, "knk_input_noMT.mtx"))
genes <- read.table(file.path(BASE, "knk_genes_noMT.tsv"), stringsAsFactors = FALSE)$V1
rownames(counts) <- genes
cat("矩阵(剔MT):", nrow(counts), "基因 x", ncol(counts), "细胞\n")

GKO <- "ENSG00000136231"
stopifnot(GKO %in% rownames(counts))

t0 <- Sys.time()
result <- scTenifoldKnk(
  countMatrix = counts, gKO = GKO, qc = FALSE,
  nc_nNet = 10, nc_nCells = 500, nc_nComp = 3,
  nCores = min(4, max(1, parallel::detectCores() - 2))
)
cat("总耗时:", round(difftime(Sys.time(), t0, units = "mins"), 1), "分钟\n")
DR <- as.data.frame(result$diffRegulation)
write.csv(DR, file.path(BASE, "B2_knockout_IGF2BP3_noMT_full.csv"), row.names = FALSE)
cat("diffRegulation 行数:", nrow(DR), "\n完成\n")
