# A1 复核：TCGA-LUAD 9 基因防御轴 单变量 Cox（独立重算）
# 数据：UCSC Xena TCGA hub HiSeqV2（RSEM）+ Xena survival 表
# 口径：log2(x+1) -> 队列内 z 分数 -> Cox 连续（每标准差 HR）
library(survival)

genes <- c("SLC7A11","SLC3A2","GCLC","GCLM","GOT2","TFRC","VDAC2","IGF2BP2","IGF2BP3")

expr <- read.delim(gzfile("TCGA-LUAD_HiSeqV2.gz"), check.names = FALSE, row.names = 1)
surv <- read.delim(gzfile("TCGA-LUAD_survival.tsv.gz"), check.names = FALSE)

# 行名形如 SYMBOL|entrez，取 SYMBOL；列名形如 TCGA-XX-XXXX-01A
sym <- sub("\\|.*$", "", rownames(expr))
keep <- sym %in% genes
expr <- expr[keep, ]; sym <- sym[keep]
rownames(expr) <- sym

samples <- colnames(expr)
sid15 <- samples
sid <- substr(samples, 1, 12)                 # TCGA-69-7978-01（12 位，无 aliquot 字母）
tum <- substr(samples, 14, 15) == "01"        # 原发肿瘤
ss  <- substr(surv$sample, 1, 12)             # TCGA-NJ-A4YI-01A -> 12 位
surv <- surv[!duplicated(ss), ]
rownames(surv) <- substr(surv$sample, 1, 12)
common <- intersect(sid[tum], rownames(surv))
m <- match(common, sid)
expr <- expr[, m]
surv <- surv[common, ]

cat("样本量:", nrow(surv), " 死亡事件:", sum(surv$OS), "\n\n")

X <- log2(as.matrix(expr) + 1)
X <- t(scale(t(X)))   # 每基因队列内 z

cat(sprintf("%-10s %6s %6s %6s %10s\n", "gene", "HR", "lo95", "hi95", "P"))
res <- data.frame(gene=genes, HR=NA, lo=NA, hi=NA, P=NA)
for (g in genes) {
  x <- X[g, ]
  ok <- complete.cases(x, surv$OS.time, surv$OS)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = surv[ok, ])
  s <- summary(fit)$conf.int
  p <- summary(fit)$coefficients[, "Pr(>|z|)"]
  res[res$gene==g, ] <- data.frame(gene=g, HR=round(s[1,1],3), lo=round(s[1,3],3),
                                   hi=round(s[1,4],3), P=signif(p, 3))
  cat(sprintf("%-10s %6.3f %6.3f %6.3f %10.3g\n", g, s[1,1], s[1,3], s[1,4], p))
}
write.csv(res, "A1_九基因单变量Cox_复核.csv", row.names = FALSE)
nsig <- sum(res$P < 0.05)
cat("\nP<0.05 的基因数:", nsig, "（", paste(res$gene[res$P<0.05], collapse=", "), "）\n")
