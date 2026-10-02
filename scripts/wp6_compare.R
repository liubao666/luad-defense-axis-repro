#!/usr/bin/env Rscript
# =====================================================================
# WP6c  免疫浸润高/低防御组比较 + 图
# 家族内 BH 校正: A=CIBERSORT 22 细胞类型; B=MCP-counter 10 群体;
#                C=5 检查点 + CYT
# =====================================================================
.libPaths(c("/Users/liubao/.workbuddy/Rlibs", .libPaths()))
suppressMessages({library(ggplot2); library(gridExtra)})

DIR <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP6_immune"

compare <- function(df, value_cols, family){
  do.call(rbind, lapply(value_cols, function(v){
    hi <- df[df$grp == "high", v]; lo <- df[df$grp == "low", v]
    hi <- hi[is.finite(hi)]; lo <- lo[is.finite(lo)]
    if(length(hi) < 10 || length(lo) < 10)
      return(data.frame(feature = v, family = family, n_high = length(hi),
                        n_low = length(lo), median_high = NA, median_low = NA,
                        diff_med = NA, p_wilcox = NA, stringsAsFactors = FALSE))
    wt <- wilcox.test(hi, lo)
    data.frame(feature = v, family = family,
               n_high = length(hi), n_low = length(lo),
               median_high = median(hi), median_low = median(lo),
               diff_med = median(hi) - median(lo),
               p_wilcox = wt$p.value, stringsAsFactors = FALSE)
  }))
}

## ---- CIBERSORT ----
cs <- read.csv(file.path(DIR, "WP6_cibersort_fractions.csv"), stringsAsFactors = FALSE)
frac_cols <- setdiff(colnames(cs), c("sample","patient","grp","P-value","P.value","RMSE","Correlation"))
resA <- compare(cs, frac_cols, "CIBERSORT")
resA$p_bh <- p.adjust(resA$p_wilcox, "BH")

## ---- MCP-counter ----
mcp <- read.csv(file.path(DIR, "WP6_mcpcounter_scores.csv"), stringsAsFactors = FALSE)
mcp_cols <- setdiff(colnames(mcp), c("patient","grp"))
resB <- compare(mcp, mcp_cols, "MCPcounter")
resB$p_bh <- p.adjust(resB$p_wilcox, "BH")

## ---- 检查点 + CYT ----
ck <- read.csv(file.path(DIR, "WP6_checkpoint_cyt.csv"), stringsAsFactors = FALSE)
ck_cols <- setdiff(colnames(ck), c("patient","grp"))
resC <- compare(ck, ck_cols, "Checkpoint_CYT")
resC$p_bh <- p.adjust(resC$p_wilcox, "BH")

res <- rbind(resA, resB, resC)
res <- res[order(res$p_wilcox), ]
write.csv(res, file.path(DIR, "WP6_immune_comparison.csv"), row.names = FALSE)
print(res[, c("feature","family","median_high","median_low","diff_med","p_wilcox","p_bh")],
      row.names = FALSE, digits = 3)

## ---- 分组一致性核对 ----
stopifnot(identical(sort(cs$patient[cs$grp=="high"]), sort(mcp$patient[mcp$grp=="high"])))
cat("\n组间样本一致性 OK | 高组:", sum(cs$grp=="high"), "低组:", sum(cs$grp=="low"), "\n")

## ---- 图 1（正文候选）: 两算法 top 差异细胞面板 ----
topA <- resA[order(resA$p_wilcox), ][1:4, "feature"]
topB <- resB[order(resB$p_wilcox), ][1:4, "feature"]
mkdf <- function(df, feats, dflab){
  do.call(rbind, lapply(feats, function(v){
    data.frame(feature = v, method = dflab, group = df$grp, value = df[[v]])
  }))
}
d1 <- rbind(mkdf(cs, topA, "CIBERSORT (LM22)"), mkdf(mcp, topB, "MCP-counter"))
d1$lab <- paste0(d1$method, "\n", d1$feature)
d1$group <- factor(d1$group, levels = c("low","high"))
ann <- res[res$feature %in% c(topA, topB), c("feature","p_wilcox","p_bh")]
ann$method <- ifelse(ann$feature %in% topA, "CIBERSORT (LM22)", "MCP-counter")
ann$lab <- paste0(ann$method, "\n", ann$feature)
ann$group <- factor("low", levels = c("low","high"))
p1 <- ggplot(d1, aes(group, value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.13, size = 0.7, alpha = 0.45, color = "grey25") +
  geom_text(data = ann, aes(x = 1.5, y = Inf,
        label = sprintf("P = %.3g\nBH = %.3g", p_wilcox, p_bh)),
        inherit.aes = FALSE, vjust = 1.4, size = 2.9) +
  facet_wrap(~lab, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(low = "#2b6cb0", high = "#c53030")) +
  labs(x = NULL, y = "Estimated fraction / score",
       title = "Immune infiltration by defense-axis score group (TCGA-LUAD, n = 502)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11.5, face = "bold"),
        strip.text = element_text(size = 8.6))
ggsave(file.path(DIR, "WP6_immune_cells.png"), p1, width = 10.5, height = 5.4, dpi = 300)
cat("saved WP6_immune_cells.png\n")

## ---- 图 2（补图）: 检查点 + CYT 六面板 ----
d2 <- mkdf(ck, ck_cols, "checkpoint")
d2$lab <- d2$feature
d2$group <- factor(d2$group, levels = c("low","high"))
ann2 <- resC; ann2$lab <- ann2$feature
ann2$group <- factor("low", levels = c("low","high"))
p2 <- ggplot(d2, aes(group, value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.13, size = 0.7, alpha = 0.45, color = "grey25") +
  geom_text(data = ann2, aes(x = 1.5, y = Inf,
        label = sprintf("P = %.3g\nBH = %.3g", p_wilcox, p_bh)),
        inherit.aes = FALSE, vjust = 1.4, size = 3.1) +
  facet_wrap(~lab, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c(low = "#2b6cb0", high = "#c53030")) +
  labs(x = NULL, y = "log2 expression",
       title = "Immune checkpoint genes and cytolytic activity (CYT) by score group") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11.5, face = "bold"),
        strip.text = element_text(size = 9))
ggsave(file.path(DIR, "WP6_checkpoints_cyt.png"), p2, width = 12, height = 4.2, dpi = 300)
cat("saved WP6_checkpoints_cyt.png\n")
cat("\nDONE\n")
