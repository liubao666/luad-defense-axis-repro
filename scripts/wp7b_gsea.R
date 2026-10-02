#!/usr/bin/env Rscript
# =====================================================================
# WP7b  GSEA：高/低防御组，Hallmark 全 50 集 + 铁死亡/KEAP1-NRF2 目标集
# 方案: 最终执行方案_二区落地_2026-10-01.md WP7
# 基因集口径（已核查，如实命名）:
#   - HALLMARK_PI3K_AKT_MTOR_SIGNALING / EMT          (MSigDB Hallmark)
#   - GOBP_FERROPTOSIS                                (GO:BP; Hallmark 无铁死亡集)
#   - REACTOME_KEAP1_NFE2L2_PATHWAY + IBRAHIM_NRF2_UP (C2; Hallmark 无 NRF2 集)
# 排序统计量: 高/低组逐基因 Welch t（向量化解算）
# =====================================================================
.libPaths(c("/Users/liubao/.workbuddy/Rlibs", .libPaths()))
suppressMessages({library(data.table); library(msigdbr); library(fgsea); library(ggplot2)})

DIR  <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP7_gsea_mutation"
EXPR <- "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/TCGA-LUAD_HiSeqV2.gz"

sam <- read.delim("/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP2_subgroup_survival/WP2_分析用样本表.tsv",
                  stringsAsFactors = FALSE)
stopifnot(nrow(sam) == 502)

E <- fread(EXPR)
genes <- E[[1]]
allc <- colnames(E)[-1]
keep <- which(substr(allc, 1, 12) %in% sam$patient & grepl("-01$", allc)) + 1L
M <- as.matrix(E[, ..keep]); rownames(M) <- genes
rm(E); gc()
dup <- duplicated(rownames(M))
if(any(dup)){
  gg <- unique(rownames(M)[dup])
  M2 <- sapply(gg, function(g) colMeans(M[rownames(M) == g, , drop = FALSE]))
  M <- rbind(M[!dup, , drop = FALSE], t(M2))
}
grp <- sam$grp[match(substr(colnames(M), 1, 12), sam$patient)]
stopifnot(!any(is.na(grp)))
cat("表达矩阵:", nrow(M), "x", ncol(M), " | 高组:", sum(grp=="high"), "低组:", sum(grp=="low"), "\n")

## ---- 逐基因 Welch t（向量化）----
X <- M[, grp == "high"]; Y <- M[, grp == "low"]
nx <- ncol(X); ny <- ncol(Y)
mx <- rowMeans(X); my <- rowMeans(Y)
vx <- rowSums((X - mx)^2) / (nx - 1); vy <- rowSums((Y - my)^2) / (ny - 1)
se <- sqrt(vx / nx + vy / ny)
tstat <- (mx - my) / se
tstat <- tstat[is.finite(tstat)]
cat("有效基因:", length(tstat), "\n")

## ---- 基因集 ----
hall <- msigdbr(species = "Homo sapiens", collection = "H")
go   <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
c2   <- msigdbr(species = "Homo sapiens", collection = "C2")
target_sets <- c("GOBP_FERROPTOSIS",
                 "REACTOME_KEAP1_NFE2L2_PATHWAY", "IBRAHIM_NRF2_UP",
                 "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                 "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
sets <- rbind(hall, go[go$gs_name %in% target_sets, ],
              c2[c2$gs_name %in% target_sets, ])
gl <- split(sets$gene_symbol, sets$gs_name)
cat("基因集总数:", length(gl), "\n")
stopifnot(all(target_sets %in% names(gl)))

## ---- fgsea ----
set.seed(42)
fg <- fgsea(pathways = gl, stats = tstat, minSize = 10, maxSize = 500,
            nperm = 10000)
fg <- fg[order(fg$pval), ]
fg_out <- data.frame(pathway = fg$pathway, NES = fg$NES,
                     p = fg$pval, padj = fg$padj,
                     size = sapply(fg$leadingEdge, length))
write.csv(fg_out, file.path(DIR, "WP7_gsea_all_sets.csv"), row.names = FALSE)
print(head(fg_out, 15), row.names = FALSE, digits = 3)
cat("\n== 四目标集 ==\n")
print(fg_out[fg_out$pathway %in% target_sets, ], row.names = FALSE, digits = 3)

## ---- 图 1：四目标集 NES 条形 ----
tt <- fg_out[fg_out$pathway %in% target_sets, ]
tt$pathway <- factor(tt$pathway, levels = tt$pathway[order(tt$NES)])
p1 <- ggplot(tt, aes(pathway, NES, fill = NES > 0)) +
  geom_col(width = 0.65) +
  geom_text(aes(y = NES / 2, label = sprintf("padj=%.3g", padj)),
            color = "white", hjust = 0.5, size = 3.2) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#c53030", `FALSE` = "#2b6cb0")) +
  labs(x = NULL, y = "Normalized enrichment score (NES), high vs low score group",
       title = "GSEA of defense-axis-associated pathways (TCGA-LUAD)") +
  theme_bw(base_size = 11.5) +
  theme(legend.position = "none", plot.title = element_text(size = 11.5, face = "bold"))
ggsave(file.path(DIR, "WP7_gsea_target_sets.png"), p1, width = 8.5, height = 4, dpi = 300)
cat("saved WP7_gsea_target_sets.png\n")

## ---- 图 2：Hallmark top10 NES 条形（全 50 集背景）----
th <- fg_out[grepl("HALLMARK", fg_out$pathway), ]
th <- th[order(th$NES), ][c(1:5, (nrow(th)-4):nrow(th)), ]
th$pathway <- factor(th$pathway, levels = th$pathway)
p2 <- ggplot(th, aes(pathway, NES, fill = NES > 0)) +
  geom_col(width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#c53030", `FALSE` = "#2b6cb0")) +
  labs(x = NULL, y = "NES",
       title = "Top/bottom Hallmark sets (of 50) by NES") +
  theme_bw(base_size = 10.5) +
  theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 8.5))
ggsave(file.path(DIR, "WP7_gsea_hallmark_top.png"), p2, width = 8, height = 4.2, dpi = 300)
cat("saved WP7_gsea_hallmark_top.png\nDONE\n")
