# 模块二（细胞子集版）：仅在 IGF2BP3+ 肿瘤细胞内虚拟敲除
# 理由：IGF2BP3 亚克隆表达（全局仅1.95%），全局网络中其调控效应被稀释；
#       在其表达的细胞子集内做 KO 才是检验「IGF2BP3 是否接线防御轴」的正确问题。
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
counts <- readMM(file.path(BASE, "knk_input_ko_cells.mtx"))
genes <- read.table(file.path(BASE, "knk_genes.tsv"), stringsAsFactors = FALSE)$V1
rownames(counts) <- genes
cat("矩阵:", nrow(counts), "基因 x", ncol(counts), "细胞 (IGF2BP3+)\n")

GKO <- "ENSG00000136231"
stopifnot(GKO %in% rownames(counts))

t0 <- Sys.time()
result <- scTenifoldKnk(
  countMatrix = counts, gKO = GKO, qc = FALSE,
  nc_nNet = 10, nc_nCells = 200, nc_nComp = 3,
  nCores = min(4, max(1, parallel::detectCores() - 2))
)
cat("总耗时:", round(difftime(Sys.time(), t0, units = "mins"), 1), "分钟\n")
DR <- as.data.frame(result$diffRegulation)
write.csv(DR, file.path(BASE, "B2_knockout_IGF2BP3_KOcells_full.csv"), row.names = FALSE)
cat("diffRegulation 行数:", nrow(DR), "\n完成\n")
