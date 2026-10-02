# WP3：9 基因防御轴评分 1/3/5 年时间依赖 ROC（TCGA-LUAD + GSE68465）
# 口径与稿件一致：ext_common.R::surv_report / ext_tcga.R / ext_geo.R
#   - 评分 = 9 基因队列内 z 分数均值（DEFENSE9 固定清单）
#   - timeROC(cause=1, times=12/36/60 月)；此处另开 iid=TRUE 求 AUC 的 95%CI
# 数据：复用论文管线已落盘的逐样本评分表（不重算表达矩阵，保证逐位一致）
.libPaths(c("/Users/liubao/.workbuddy/Rlibs", .libPaths()))
suppressMessages({library(survival); library(timeROC)})

DIR  <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP3_timeROC"
EXT  <- "/Users/liubao/WorkBuddy/2026-09-23-10-23-05/生信学习计划/第16课数据/extval"
set.seed(42)

# ---- 读取两队列逐样本 (time, event, defense_score) ----
tcga <- read.csv(file.path(EXT, "TCGA_LUAD_OS_评分表.csv"), stringsAsFactors = FALSE)
gse  <- readRDS(file.path(EXT, "GSE68465_验证对象.rds"))
gse_d <- gse$out$def_OS$data[[1]]          # time(月), event, score
d_tcga <- data.frame(time = tcga$time, event = tcga$event, score = tcga$defense_score)
d_gse  <- data.frame(time = gse_d$time,  event = gse_d$event, score = gse_d$score)
cat("TCGA-LUAD OS: n=", nrow(d_tcga), " 事件=", sum(d_tcga$event), "\n", sep = "")
cat("GSE68465  OS: n=", nrow(d_gse),  " 事件=", sum(d_gse$event),  "\n", sep = "")

TIMES <- c(12, 36, 60)
COLS  <- c("12" = "#2b6cb0", "36" = "#c53030", "60" = "#2f855a")

run_roc <- function(d, label){
  d <- d[is.finite(d$time) & d$time > 0 & !is.na(d$event) & is.finite(d$score), ]
  tr <- timeROC(T = d$time, delta = d$event, marker = d$score,
                cause = 1, times = TIMES, iid = TRUE)
  tri <- timeROC(T = d$time, delta = d$event, marker = d$score,
                 cause = 1, times = TIMES, iid = FALSE)
  cat(sprintf("%s  AUC 1/3/5y = %.3f / %.3f / %.3f\n", label, tr$AUC[1], tr$AUC[2], tr$AUC[3]))
  list(d = d, tr = tr, tri = tri, label = label)
}

# AUC 95%CI（Blanche et al. 2013 iid 分解）
auc_ci <- function(tr){
  se <- tr$inference$vect_sd_1
  data.frame(year = c("1-year","3-year","5-year"),
             AUC  = round(tr$AUC, 3),
             lo   = round(pmax(0, tr$AUC - 1.96*se), 3),
             hi   = round(pmin(1, tr$AUC + 1.96*se), 3),
             row.names = NULL)
}

r_tcga <- run_roc(d_tcga, "TCGA-LUAD")
r_gse  <- run_roc(d_gse,  "GSE68465")

# ---- 一致性核对：与论文管线 iid=FALSE 点估计必须逐位一致 ----
stopifnot(all.equal(unname(round(r_tcga$tri$AUC, 3)), c(0.614, 0.618, 0.567)))
stopifnot(all.equal(unname(round(r_gse$tri$AUC, 3)),  c(0.654, 0.659, 0.601)))
cat("一致性核对通过：与 TCGA_LUAD_验证结果.txt / GEO_验证汇总.rds 逐位一致\n")

# ---- 汇总表 ----
summ <- rbind(data.frame(cohort = "TCGA-LUAD (OS)", n = nrow(r_tcga$d), events = sum(r_tcga$d$event), auc_ci(r_tcga$tr)),
              data.frame(cohort = "GSE68465 (OS)",  n = nrow(r_gse$d),  events = sum(r_gse$d$event),  auc_ci(r_gse$tr)))
write.csv(summ, file.path(DIR, "WP3_timeROC_AUC汇总.csv"), row.names = FALSE)
print(summ)

# ---- 图 1/2：每个队列一张 ROC 三联曲线 ----
draw <- function(r, fname, ttl){
  png(file.path(DIR, fname), width = 1600, height = 1500, res = 220)
  par(mar = c(4.5, 4.5, 3, 1), cex.lab = 1.15, cex.axis = 1.05)
  plot(c(0,1), c(0,1), type = "n", xlab = "1 - Specificity (FPF)",
       ylab = "Sensitivity (TPF)", main = ttl, cex.main = 1.3, asp = 1)
  abline(0, 1, lty = 2, col = "grey55")
  for(i in seq_along(TIMES)){
    x <- r$tr$FP[, i]; y <- r$tr$TP[, i]
    lines(x, y, lwd = 2.2, col = COLS[i])
  }
  legend("bottomright", bty = "n", lwd = 2.2, col = COLS,
         legend = sprintf("%s  AUC = %.3f (%.3f-%.3f)", c("1-year","3-year","5-year"),
                          r$tr$AUC, auc_ci(r$tr)$lo, auc_ci(r$tr)$hi), cex = 1.05)
  legend("topleft", bty = "n", legend = sprintf("n = %d, deaths = %d",
         nrow(r$d), sum(r$d$event)), cex = 1.05)
  dev.off()
  cat("saved", fname, "\n")
}
draw(r_tcga, "WP3_timeROC_TCGA.png", "TCGA-LUAD (OS)\ntime-dependent ROC of 9-gene defence-axis score")
draw(r_gse,  "WP3_timeROC_GSE68465.png", "GSE68465 (OS)\ntime-dependent ROC of 9-gene defence-axis score")
cat("DONE\n")
