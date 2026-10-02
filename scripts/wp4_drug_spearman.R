# WP4：LUAD 细胞系 9 基因防御轴评分 × 5 药 IC50/AUC（Spearman）
# 方案：最终执行方案_二区落地_2026-10-01.md WP4 行
# 口径：
#   评分 = DEFENSE9 在 LUAD 模型(DepmapModelType=='LUAD')内 z 分数均值（与全文队列内 z 一致）
#   GDSC：TCGA_DESC=='LUAD'，药物 Torin 2 / Dactolisib(=BEZ235) / MK-2206 / ZSTK474
#         细胞系名规范化匹配 Model.csv StrippedCellLineName
#   PRISM：5 药 torin-2 / NVP-BEZ235 / MK-2206 / ZSTK-474 / MLN0128(=sapanisertib)
#         按 depmap_id 直连；IC50/AUC 双口径（MTS010 批次 IC50 缺失多，WP0 已实锤）
# 预期方向：高分 → 低 IC50（更敏感）→ ρ(IC50) < 0
.libPaths(c("/Users/liubao/.workbuddy/Rlibs", .libPaths()))
TMP  <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP0_coverage/tmp"
DIR  <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP4_drug_sensitivity"
set.seed(42)

GENES <- c("SLC7A11","SLC3A2","GCLC","GCLM","GOT2","TFRC","VDAC2","IGF2BP2","IGF2BP3")

## ---- 1. 读表达矩阵（行为模型 ACH-xxxxx，列为基因），取 9 基因 ----
suppressMessages(library(data.table))
E <- fread(file.path(TMP, "OmicsExpressionProteinCodingGenesTPMLogp1.csv"))
genes_in <- sub(" \\(.*$", "", colnames(E))
gcol <- vapply(GENES, function(g) which(genes_in == g)[1], integer(1))
cat("基因列命中:", sum(!is.na(gcol)), "/9\n")
mid <- E[[1]]
idx <- gcol[!is.na(gcol)]
sub <- as.data.frame(E[, ..idx])
Em <- as.matrix(sub); colnames(Em) <- GENES[!is.na(gcol)]; rownames(Em) <- mid
cat("表达矩阵:", nrow(Em), "模型 x", ncol(Em), "基因\n")

## ---- 2. LUAD 模型 + 队列内 z 评分 ----
model <- read.csv(file.path(TMP, "Model_23Q4.csv"), stringsAsFactors = FALSE)
luad <- model$ModelID[model$DepmapModelType == "LUAD"]
El <- Em[intersect(luad, rownames(Em)), , drop = FALSE]
Z <- scale(El)                       # 每基因（列）在 LUAD 队列内 z
score <- rowMeans(Z)
sc <- data.frame(model_id = names(score), score = as.numeric(score))
cat("LUAD 模型评分:", nrow(sc), "（9 基因全部命中:", ncol(Z) == 9, "）\n")
write.csv(sc, file.path(DIR, "WP4_cellline_scores.csv"), row.names = FALSE)

norm_name <- function(x) gsub("[^A-Z0-9]", "", toupper(x))
s2m <- setNames(model$ModelID, norm_name(model$StrippedCellLineName))

## ---- 3. GDSC Spearman ----
library_gdsc <- c("Torin 2","Dactolisib","MK-2206","ZSTK474")
gdsc_res <- list()
for(fn in c("GDSC1_fitted_dose_response.csv","GDSC2_fitted_dose_response.csv")){
  ds <- if(grepl("GDSC1", fn)) "GDSC1" else "GDSC2"
  d <- read.csv(file.path(TMP, fn), stringsAsFactors = FALSE)
  d <- d[d$TCGA_DESC == "LUAD" & d$DRUG_NAME %in% library_gdsc, ]
  d$model_id <- s2m[norm_name(d$CELL_LINE_NAME)]
  d <- d[!is.na(d$model_id), ]
  d$score <- sc$score[match(d$model_id, sc$model_id)]
  d <- d[!is.na(d$score), ]
  for(dr in library_gdsc){
    dd <- d[d$DRUG_NAME == dr & !is.na(d$LN_IC50) & d$LN_IC50 != "", ]
    if(nrow(dd) < 6) next
    ic50 <- suppressWarnings(as.numeric(dd$LN_IC50))
    auc  <- suppressWarnings(as.numeric(dd$AUC))
    r1 <- cor.test(dd$score, ic50, method = "spearman", exact = FALSE)
    okA <- !is.na(auc)
    r2 <- if(sum(okA) >= 6) cor.test(dd$score[okA], auc[okA], method = "spearman", exact = FALSE) else NULL
    gdsc_res[[paste(ds, dr)]] <- data.frame(
      dataset = ds, drug = dr, n = nrow(dd),
      rho_IC50 = unname(r1$estimate), p_IC50 = r1$p.value,
      rho_AUC = if(is.null(r2)) NA else unname(r2$estimate),
      p_AUC   = if(is.null(r2)) NA else r2$p.value)
  }
}
gdsc_df <- do.call(rbind, gdsc_res)

## ---- 4. PRISM Spearman ----
prism <- read.csv(file.path(TMP, "PRISM_secondary_dose_response.csv"), stringsAsFactors = FALSE)
targets <- c("torin-2","NVP-BEZ235","MK-2206","ZSTK-474","MLN0128")
prism <- prism[prism$name %in% targets & prism$depmap_id %in% luad, ]
prism$score <- sc$score[match(prism$depmap_id, sc$model_id)]
prism <- prism[!is.na(prism$score), ]
prism_res <- list()
for(dr in targets){
  for(scr in unique(prism$screen_id[prism$name == dr])){
    dd <- prism[prism$name == dr & prism$screen_id == scr, ]
    ic50 <- suppressWarnings(as.numeric(dd$ic50))
    auc  <- suppressWarnings(as.numeric(dd$auc))
    okI <- is.finite(ic50); okA <- is.finite(auc)
    r1 <- if(sum(okI) >= 6) cor.test(dd$score[okI], ic50[okI], method = "spearman", exact = FALSE) else NULL
    r2 <- if(sum(okA) >= 6) cor.test(dd$score[okA], auc[okA], method = "spearman", exact = FALSE) else NULL
    prism_res[[paste(dr, scr)]] <- data.frame(
      dataset = paste("PRISM", scr), drug = dr, n = nrow(dd),
      n_IC50 = sum(okI),
      rho_IC50 = if(is.null(r1)) NA else unname(r1$estimate),
      p_IC50   = if(is.null(r1)) NA else r1$p.value,
      rho_AUC = if(is.null(r2)) NA else unname(r2$estimate),
      p_AUC   = if(is.null(r2)) NA else r2$p.value)
  }
}
prism_df <- do.call(rbind, prism_res)

## ---- 5. 合并落盘 ----
gdsc_df$n_IC50 <- gdsc_df$n
res <- rbind(gdsc_df[, c("dataset","drug","n","n_IC50","rho_IC50","p_IC50","rho_AUC","p_AUC")],
             prism_df[, c("dataset","drug","n","n_IC50","rho_IC50","p_IC50","rho_AUC","p_AUC")])
res$p_IC50_adj <- p.adjust(res$p_IC50, "BH")
res$p_AUC_adj  <- p.adjust(res$p_AUC, "BH")
write.csv(res, file.path(DIR, "WP4_spearman_results.csv"), row.names = FALSE)
print(res, digits = 3)

## ---- 6. 图 1：五药散点（每药取 IC50 数据最全的库）----
CANON <- c("torin-2"="Torin-2","Torin 2"="Torin-2","NVP-BEZ235"="BEZ235",
           "Dactolisib"="BEZ235","MK-2206"="MK-2206","ZSTK-474"="ZSTK474",
           "ZSTK474"="ZSTK474","MLN0128"="sapanisertib")
res$drugc <- CANON[res$drug]
best <- do.call(rbind, lapply(split(res, res$drugc), function(x){
  x <- x[order(-x$n_IC50), ]; x[1, ]
}))
best <- best[c("Torin-2","BEZ235","MK-2206","ZSTK474","sapanisertib"), ]
get_pts <- function(b){
  if(b$dataset %in% c("GDSC1","GDSC2")){
    fn <- if(b$dataset=="GDSC1") "GDSC1_fitted_dose_response.csv" else "GDSC2_fitted_dose_response.csv"
    d <- read.csv(file.path(TMP, fn), stringsAsFactors = FALSE)
    d <- d[d$TCGA_DESC=="LUAD" & d$DRUG_NAME==b$drug, ]
    d$model_id <- s2m[norm_name(d$CELL_LINE_NAME)]
    d <- d[!is.na(d$model_id), ]
    d$score <- sc$score[match(d$model_id, sc$model_id)]
    d <- d[!is.na(d$score) & !is.na(d$LN_IC50) & d$LN_IC50 != "", ]
    data.frame(score = d$score, y = suppressWarnings(as.numeric(d$LN_IC50)))
  } else {
    scr <- sub("PRISM ", "", b$dataset)
    p <- prism[prism$name == b$drug & prism$screen_id == scr, ]
    ic50 <- suppressWarnings(as.numeric(p$ic50))
    ok <- is.finite(ic50)
    data.frame(score = p$score[ok], y = log10(ic50[ok]))
  }
}
png(file.path(DIR, "WP4_score_vs_IC50.png"), width = 2400, height = 1100, res = 200)
par(mfrow = c(1, 5), mar = c(4.2, 4.2, 3, 1), cex.lab = 1.1, cex.axis = 1.0, cex.main = 1.15)
for(i in seq_len(nrow(best))){
  b <- best[i, ]
  pts <- get_pts(b)
  plot(pts$score, pts$y, pch = 19, col = rgb(0.17,0.42,0.69,0.75), cex = 1.1,
       xlab = "9-gene defence-axis score", ylab = "ln(IC50), µM",
       main = sprintf("%s (%s)", b$drugc, b$dataset))
  abline(lm(y ~ score, pts), col = "#c53030", lwd = 2)
  rho <- b$rho_IC50
  txt <- sprintf("n=%d\nSpearman ρ = %.2f\nP = %.3g", nrow(pts), rho, b$p_IC50)
  legend("topleft", bty = "n", legend = txt, cex = 0.95, adj = c(0, 1))
}
dev.off()
cat("saved WP4_score_vs_IC50.png\n")

## ---- 7. 图 2：ρ 森林汇总（IC50 与 AUC 双行）----
rr <- res[order(res$drugc), ]
rr <- rr[!is.na(rr$rho_IC50) | !is.na(rr$rho_AUC), ]
rr$lab <- paste0(rr$drugc, " | ", rr$dataset)
n_row <- nrow(rr)
png(file.path(DIR, "WP4_forest_rho.png"), width = 1800, height = 200 + 60*n_row, res = 200)
par(mar = c(4, 13, 3, 1))
plot(rr$rho_IC50, seq_len(n_row), pch = 19, cex = 1.2, xlim = c(-0.6, 0.6),
     xlab = "Spearman ρ (score vs sensitivity)", ylab = "", yaxt = "n",
     main = "Defence-axis score vs drug sensitivity in LUAD cell lines")
abline(v = 0, lty = 2, col = "grey50")
points(rr$rho_AUC, seq_len(n_row) + 0.25, pch = 17, cex = 1.0, col = "#c53030")
axis(2, at = seq_len(n_row) + 0.12, labels = rr$lab, las = 2, cex.axis = 0.85)
legend("bottomleft", bty = "n", pch = c(19, 17), col = c("black", "#c53030"),
       legend = c("IC50-based", "AUC-based"), cex = 0.95)
dev.off()
cat("saved WP4_forest_rho.png\n")
cat("DONE\n")
