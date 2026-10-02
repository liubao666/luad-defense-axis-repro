#!/usr/bin/env Rscript
# =====================================================================
# WP5  DepMap 依赖性分析（修订 v2 分组比较设计）
# 方案: 最终执行方案_二区落地_2026-10-01.md WP5
# 设计: 按 9 基因防御轴评分将 DepMap 23Q4 LUAD 细胞系分高/低组，
#       比较两组 Chronos 敲除依赖性（重点 IGF2BP3、SLC7A11）；
#       不看单基因绝对必需性（无区分度，WP0 已验证）。
# 假说方向: 防御轴若是依赖，高防御组对轴基因更必需 -> Chronos 更低。
# 口径: DepMap Public 23Q4; LUAD = DepmapModelType=="LUAD";
#       评分沿用 WP4（队列内 z 均值, 同源全文）。
# =====================================================================
.libPaths(c("/Users/liubao/.workbuddy/Rlibs", .libPaths()))
suppressMessages({library(data.table); library(ggplot2)})

DIR  <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP5_dependency"
TMP  <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP0_coverage/tmp"
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

GENES <- c("SLC7A11","SLC3A2","GCLC","GCLM","GOT2","TFRC","VDAC2","IGF2BP2","IGF2BP3")

## ---- 1. 评分（沿用 WP4 产物）----
sc <- read.csv("/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP4_drug_sensitivity/WP4_cellline_scores.csv",
               stringsAsFactors = FALSE)
sc <- sc[is.finite(sc$score), ]
cat("WP4 评分模型:", nrow(sc), "\n")

## ---- 2. 读 Chronos（data.table，整列子集）----
E <- fread(file.path(TMP, "CRISPRGeneEffect_23Q4.csv"))
genes_in <- sub(" \\(.*$", "", colnames(E))
gcol <- vapply(GENES, function(g) which(genes_in == g)[1], integer(1))
cat("基因列命中:", sum(!is.na(gcol)), "/9\n")
stopifnot(all(!is.na(gcol)))
mid <- E[[1]]
Em <- as.matrix(E[, ..gcol]); colnames(Em) <- GENES; rownames(Em) <- mid
rm(E); gc()

## ---- 3. 交集：有评分且有 Chronos 的 LUAD 模型 ----
common <- intersect(sc$model_id, rownames(Em))
score <- sc$score[match(common, sc$model_id)]; names(score) <- common
dep <- Em[common, , drop = FALSE]
cat("交集模型数:", length(common), "\n")

grp <- ifelse(score > median(score, na.rm = TRUE), "high", "low")
cat("高防御组:", sum(grp=="high"), "| 低防御组:", sum(grp=="low"), "\n")
write.csv(data.frame(model_id = common, score = score, group = grp),
          file.path(DIR, "WP5_group_assignment.csv"), row.names = FALSE)

## ---- 4. 分组比较：每基因 Wilcoxon + 中位数差 + BH 校正 ----
res <- do.call(rbind, lapply(GENES, function(g){
  v <- dep[, g]
  hi <- v[grp == "high"]; lo <- v[grp == "low"]
  wt <- wilcox.test(hi, lo)
  data.frame(gene = g,
             n_high = length(hi), n_low = length(lo),
             median_high = median(hi), median_low = median(lo),
             diff_med = median(hi) - median(lo),
             pct_high_le_minus1 = mean(hi <= -1) * 100,
             pct_low_le_minus1  = mean(lo <= -1) * 100,
             p_wilcox = wt$p.value,
             stringsAsFactors = FALSE)
}))
res$p_bh <- p.adjust(res$p_wilcox, method = "BH")
res <- res[order(res$p_wilcox), ]
print(res, row.names = FALSE)
write.csv(res, file.path(DIR, "WP5_group_comparison.csv"), row.names = FALSE)

## 内部一致性核对（防口径跑偏）
stopifnot(all(c("IGF2BP3","SLC7A11") %in% res$gene))

## ---- 5. 图 1（正文候选）: IGF2BP3 + SLC7A11 双面板 ----
key <- c("IGF2BP3","SLC7A11")
dd <- do.call(rbind, lapply(key, function(g){
  data.frame(gene = g, group = grp, chronos = dep[, g])
}))
dd$gene <- factor(dd$gene, levels = key)
dd$group <- factor(dd$group, levels = c("low","high"))
pdat <- merge(dd, res, by = "gene")
lab_df <- data.frame(
  gene = factor(key, levels = key),
  group = factor("low", levels = c("low","high")),
  lab = vapply(key, function(g){
    r <- res[res$gene == g, ]
    sprintf("Δ median = %.2f\nWilcoxon P = %.3g\nBH-adjusted P = %.3g",
            r$diff_med, r$p_wilcox, r$p_bh)
  }, character(1)))
p1 <- ggplot(dd, aes(group, chronos, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.13, size = 1.6, alpha = 0.6, color = "grey20") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_hline(yintercept = -1, linetype = 3, color = "grey70") +
  geom_text(data = lab_df, aes(x = 1.5, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.5, size = 3.4) +
  facet_wrap(~gene, scales = "free_y") +
  scale_fill_manual(values = c(low = "#2b6cb0", high = "#c53030")) +
  labs(x = NULL, y = "Chronos knockout dependency score",
       title = "CRISPR dependency of defense-axis genes in LUAD cell lines") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"))
ggsave(file.path(DIR, "WP5_dependency_IGF2BP3_SLC7A11.png"), p1,
       width = 7.5, height = 4.2, dpi = 300)
cat("saved WP5_dependency_IGF2BP3_SLC7A11.png\n")

## ---- 6. 图 2（补图）: 9 基因全面板 ----
dd9 <- do.call(rbind, lapply(GENES, function(g){
  data.frame(gene = g, group = grp, chronos = dep[, g])
}))
dd9$gene <- factor(dd9$gene, levels = GENES)
dd9$group <- factor(dd9$group, levels = c("low","high"))
p2 <- ggplot(dd9, aes(group, chronos, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.9, alpha = 0.5, color = "grey30") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  facet_wrap(~gene, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = c(low = "#2b6cb0", high = "#c53030")) +
  labs(x = NULL, y = "Chronos score",
       title = "All nine defense-axis genes, high vs low score groups") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 0))
ggsave(file.path(DIR, "WP5_dependency_9genes.png"), p2,
       width = 11, height = 5.4, dpi = 300)
cat("saved WP5_dependency_9genes.png\n")

cat("\nDONE\n")
