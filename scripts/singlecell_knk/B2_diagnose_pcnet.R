# 诊断：pcNet 单次耗时（Rcpp vs 纯R回退）
.libPaths(c("~/Rlibs", .libPaths()))
library(scTenifoldNet)
library(Matrix)

set.seed(1)
for (g in c(500, 1000, 2000)) {
  X <- matrix(rpois(500 * g, 3), nrow = g)
  rownames(X) <- paste0("G", seq_len(g))
  t1 <- system.time({
    n1 <- pcNet(X, nComp = 3, useRcpp = TRUE)
  })[["elapsed"]]
  cat(sprintf("genes=%4d Rcpp=TRUE  : %.2f s\n", g, t1))
  if (g <= 1000) {
    t2 <- system.time({
      n2 <- pcNet(X, nComp = 3, useRcpp = FALSE, nCores = 1)
    })[["elapsed"]]
    cat(sprintf("genes=%4d Rcpp=FALSE : %.2f s\n", g, t2))
  }
}
cat("Rcpp 函数存在性: ", exists("pcNetCoreRcpp"), "\n")
