# 模块二：scTenifoldKnk 虚拟敲除 IGF2BP3（聚焦 824 基因面板）
# 工程要点：pcNet 默认走 Rcpp 后端（实测比纯R慢8倍），故在命名空间内打补丁强制 useRcpp=FALSE。
.libPaths(c("~/Rlibs", .libPaths()))
library(scTenifoldNet)
library(scTenifoldKnk)
library(Matrix)
library(future)

# --- 强制纯 R 后端 + 限核，避免 Rcpp 慢路径与多进程内存爆炸 ---
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
counts <- readMM(file.path(BASE, "knk_input.mtx"))
genes <- read.table(file.path(BASE, "knk_genes.tsv"), stringsAsFactors = FALSE)$V1
rownames(counts) <- genes
cat("矩阵:", nrow(counts), "基因 x", ncol(counts), "细胞\n")

GKO <- "ENSG00000136231"  # IGF2BP3
stopifnot(GKO %in% rownames(counts))

t0 <- Sys.time()
result <- scTenifoldKnk(
  countMatrix = counts,
  gKO       = GKO,
  qc        = FALSE,   # Python 端已过滤；默认 minPCT=0.05 会剔除非克隆表达的 IGF2BP3
  nc_nNet   = 10,
  nc_nCells = 500,
  nc_nComp  = 3,
  nCores    = min(4, max(1, parallel::detectCores() - 2))
)
cat("总耗时:", round(difftime(Sys.time(), t0, units = "mins"), 1), "分钟\n")

DR <- as.data.frame(result$diffRegulation)
print(head(DR, 20))
write.csv(DR, file.path(BASE, "B2_knockout_IGF2BP3_full.csv"), row.names = FALSE)
cat("diffRegulation 行数:", nrow(DR), "\n")
cat("完成\n")
