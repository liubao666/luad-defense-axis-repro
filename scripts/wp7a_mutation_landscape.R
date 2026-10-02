#!/usr/bin/env Rscript
# =====================================================================
# WP7a  突变景观：TMB + KEAP1/NFE2L2/TP53/KRAS 频率，高/低防御组比较
# 数据: GDC API (ssm_occurrences, TCGA-LUAD)
#   - tmb_counts.tsv         每例全突变计数 (177330 条, 557 例)
#   - tp53_kras_status.tsv   TP53/KRAS 阳性病例
#   - WP2 mutation_status.tsv  KEAP1/NFE2L2 阳性病例 (WP2 产物)
# 注: GDC 突变注释覆盖 557 例；"未在阳性清单"按 WT 处理（全队列口径一致）
# =====================================================================
.libPaths(c("/Users/liubao/.workbuddy/Rlibs", .libPaths()))
suppressMessages({library(ggplot2)})

DIR <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP7_gsea_mutation"
WP2 <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP2_subgroup_survival"

sam <- read.delim(file.path(WP2, "WP2_分析用样本表.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(sam) == 502)

mut2 <- read.delim(file.path(WP2, "mutation_status.tsv"), stringsAsFactors = FALSE)
tk <- read.delim(file.path(DIR, "tp53_kras_status.tsv"), stringsAsFactors = FALSE)
tmb <- read.delim(file.path(DIR, "tmb_counts.tsv"), stringsAsFactors = FALSE)

## 每基因阳性集合（12位患者ID）
pos <- list(
  KEAP1  = mut2$patient[mut2$KEAP1_mut == 1],
  NFE2L2 = mut2$patient[mut2$NFE2L2_mut == 1],
  TP53   = tk$patient[tk$TP53_mut == 1],
  KRAS   = tk$patient[tk$KRAS_mut == 1]
)

## 合并到 502 例
df <- data.frame(patient = sam$patient, grp = sam$grp, stringsAsFactors = FALSE)
for(g in names(pos)) df[[g]] <- as.integer(df$patient %in% pos[[g]])
df$tmb <- tmb$mutation_count[match(df$patient, tmb$patient)]
cat("TMB 匹配率:", sum(!is.na(df$tmb)), "/502\n")

## Fisher 精确检验（高/低 × mut/WT）
res <- do.call(rbind, lapply(names(pos), function(g){
  tab <- table(df$grp, df[[g]])
  ft <- fisher.test(tab)
  hi_n <- sum(df$grp == "high"); lo_n <- sum(df$grp == "low")
  data.frame(gene = g,
             freq_high = sprintf("%d/%d (%.1f%%)", tab["high","1"], hi_n, tab["high","1"]/hi_n*100),
             freq_low  = sprintf("%d/%d (%.1f%%)", tab["low","1"], lo_n, tab["low","1"]/lo_n*100),
             OR = 1/unname(ft$estimate), or_lo = 1/ft$conf.int[2], or_hi = 1/ft$conf.int[1],
             p_fisher = ft$p.value, stringsAsFactors = FALSE)
}))
res$p_bh <- p.adjust(res$p_fisher, "BH")
print(res, row.names = FALSE, digits = 3)
write.csv(res, file.path(DIR, "WP7_mutation_fisher.csv"), row.names = FALSE)

## TMB 比较
wt <- wilcox.test(tmb ~ grp, df)
cat(sprintf("\nTMB: 高组中位 %.0f vs 低组中位 %.0f | Wilcoxon P = %.3g\n",
            median(df$tmb[df$grp=="high"], na.rm=TRUE),
            median(df$tmb[df$grp=="low"], na.rm=TRUE), wt$p.value))
tmb_res <- data.frame(feature = "TMB (mutation count)",
                      median_high = median(df$tmb[df$grp=="high"], na.rm=TRUE),
                      median_low = median(df$tmb[df$grp=="low"], na.rm=TRUE),
                      p_wilcox = wt$p.value)
write.csv(tmb_res, file.path(DIR, "WP7_tmb_comparison.csv"), row.names = FALSE)

## ---- 图：突变频率 + TMB ----
fr <- do.call(rbind, lapply(names(pos), function(g){
  data.frame(gene = g, grp = c("low","high"),
             freq = c(mean(df[[g]][df$grp=="low"])*100, mean(df[[g]][df$grp=="high"])*100))
}))
fr$grp <- factor(fr$grp, levels = c("low","high"))
lab_or <- sprintf("OR=%.2g, P=%.3g", res$OR, res$p_fisher)
names(lab_or) <- res$gene
p1 <- ggplot(fr, aes(gene, freq, fill = grp)) +
  geom_col(position = position_dodge(0.7), width = 0.62, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", freq)), position = position_dodge(0.7),
            vjust = -0.4, size = 3) +
  annotate("text", x = 1:4, y = 70, size = 3.2,
           label = lab_or[levels(factor(fr$gene))]) +
  scale_fill_manual(values = c(low = "#2b6cb0", high = "#c53030")) +
  coord_cartesian(ylim = c(0, 74), clip = "off") +
  labs(x = NULL, y = "Mutation frequency (%)",
       title = "Mutation landscape by defense-axis score group (TCGA-LUAD, n = 502)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"), legend.title = element_blank(),
        plot.margin = margin(18, 8, 8, 8))
ggsave(file.path(DIR, "WP7_mutation_freq.png"), p1, width = 7.5, height = 4.6, dpi = 300)

d2 <- df[!is.na(df$tmb), ]
p2 <- ggplot(d2, aes(grp, tmb, fill = grp)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.4, color = "grey25") +
  annotate("text", x = 1.5, y = Inf, vjust = 1.4, size = 3.6,
           label = sprintf("Wilcoxon P = %.3g", wt$p.value)) +
  scale_fill_manual(values = c(low = "#2b6cb0", high = "#c53030")) +
  scale_y_log10() +
  labs(x = NULL, y = "Mutation count per case (log scale)",
       title = "Tumour mutation burden proxy by score group") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(DIR, "WP7_tmb_boxplot.png"), p2, width = 5.2, height = 4.6, dpi = 300)

write.csv(df, file.path(DIR, "WP7_mutation_table.csv"), row.names = FALSE)
cat("saved WP7_mutation_freq.png / WP7_tmb_boxplot.png\nDONE\n")
