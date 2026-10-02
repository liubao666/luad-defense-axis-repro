# WP7.5 多变量 Cox 混杂校正（WP7 遗留关键项）
# 目的：校正年龄、病理分期、TMB、KEAP1 突变、TP53 突变后，9 基因防御轴评分对 OS 的独立预后价值
# 依据：WP7 发现高防御组 TP53 突变富集（60.6% vs 37.1%）且 TMB 近翻倍（280 vs 154），
#       存在吸烟等混杂风险，审稿人很可能要求多变量校正
# 口径：score = 队列内 9 基因 z 均值（稿件 §2.5）；KM 时间单位 = 月（OS.time/30.4375）
# 输入：WP2_分析用样本表.tsv + tmb_counts.tsv + mutation_status.tsv + tp53_kras_status.tsv
# 输出：WP7.5_多变量Cox结果.tsv、WP7.5_多变量森林图.png
.libPaths("/Users/liubao/.workbuddy/Rlibs")
library(survival)

DIR <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP2_subgroup_survival"
OUT <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP7.5_multivariable_cox"
dir.create(OUT, showWarnings = FALSE)

df  <- read.delim(file.path(DIR,"WP2_分析用样本表.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
tmb <- read.delim("/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP7_gsea_mutation/tmb_counts.tsv", check.names=FALSE, stringsAsFactors=FALSE)
mut <- read.delim(file.path(DIR,"mutation_status.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
tp  <- read.delim("/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP7_gsea_mutation/tp53_kras_status.tsv", check.names=FALSE, stringsAsFactors=FALSE)

# ---- 合并协变量 ----
rownames(tmb) <- tmb$patient
rownames(mut) <- mut$patient
rownames(tp)  <- tp$patient
df$TMB        <- as.numeric(tmb[df$patient, "mutation_count"])
df$KEAP1_mut  <- as.numeric(ifelse(df$patient %in% rownames(mut), mut[df$patient,"KEAP1_mut"], 0))
df$TP53_mut   <- as.numeric(ifelse(df$patient %in% rownames(tp),  tp[df$patient,"TP53_mut"],  0))
df$grp <- factor(df$grp, levels=c("low","high"))
df$TMBlog     <- log2(df$TMB + 1)
df$stage4     <- ifelse(df$stageGrp=="III-IV", 1, 0)   # I-II 为参照

# ---- 口径核对断言 ----
stopifnot(nrow(df) == 502)
uv <- summary(coxph(Surv(OSmo, OS) ~ grp, data=df))$conf.int
stopifnot(abs(uv[1,1] - 1.81) < 0.02, uv[1,3] < 1.36, uv[1,4] > 2.42)
cat(sprintf("断言通过：单变量 HR=%.3f (%.2f-%.2f) 与 WP2/稿件一致\n", uv[1,1], uv[1,3], uv[1,4]))

# ---- 缺失情况 ----
cat("TMB 缺失:", sum(is.na(df$TMB)), " age 缺失:", sum(is.na(df$age)),
    " stageGrp 缺失:", sum(is.na(df$stageGrp)), "\n")
dd <- df[!is.na(df$OS) & !is.na(df$age) & !is.na(df$stageGrp) & !is.na(df$TMB), ]
cat("进入多变量模型的样本:", nrow(dd), " 事件:", sum(dd$OS==1), "\n\n")

# ---- 模型 1：评分二分（与全文 KM 口径一致）----
m1 <- coxph(Surv(OSmo, OS) ~ grp + age + stage4 + TMBlog + KEAP1_mut + TP53_mut, data=dd)
# ---- 模型 2：评分连续（每 1 SD）----
m2 <- coxph(Surv(OSmo, OS) ~ score + age + stage4 + TMBlog + KEAP1_mut + TP53_mut, data=dd)
# ---- 模型 3：评分二分 + 评分×分期交互（WP2 标注的可选检验）----
m3 <- coxph(Surv(OSmo, OS) ~ grp*factor(stageGrp) + age + TMBlog + KEAP1_mut + TP53_mut, data=dd)

extract <- function(fit, label){
  s  <- summary(fit)
  ci <- s$conf.int
  data.frame(model=label,
             term=rownames(ci),
             HR=round(ci[,1],3), lo=round(ci[,3],3), hi=round(ci[,4],3),
             P=s$coefficients[rownames(ci),"Pr(>|z|)"],
             row.names=NULL)
}
res <- rbind(extract(m1,"M1 二分评分校正"), extract(m2,"M2 连续评分校正"))
res$P <- formatC(res$P, format="g", digits=3)
write.table(res, file.path(OUT,"WP7.5_多变量Cox结果.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

s3 <- summary(m3)
cat("==== M1 二分评分多变量 ====\n"); print(summary(m1)$coefficients)
cat("==== M2 连续评分多变量 ====\n"); print(summary(m2)$coefficients)
cat("==== M3 评分×分期交互 ====\n")
print(s3$coefficients[grep(":", rownames(s3$coefficients)), , drop=FALSE])

# ---- 结果数值断言（防止口径跑偏）----
m1s <- summary(m1)$coefficients
grp_p <- m1s["grphigh","Pr(>|z|)"]
cat(sprintf("\n校正后 grphigh HR=%.3f P=%.3g；score 连续 HR/SD=%.3f\n",
            summary(m1)$conf.int["grphigh",1], grp_p,
            summary(m2)$conf.int["score",1]))
stopifnot(summary(m1)$conf.int["grphigh",1] > 1)  # 方向必须保持有害

# ---- 森林图（英文标注，投稿口径）----
terms_eng <- c("grphigh"="9-gene score (high vs low)",
               "score"="9-gene score (per SD)",
               "age"="Age (per year)",
               "stage4"="Stage III–IV (vs I–II)",
               "TMBlog"="TMB (log2 count, per unit)",
               "KEAP1_mut"="KEAP1 mutation",
               "TP53_mut"="TP53 mutation")
r1 <- res[res$model=="M1 二分评分校正" & res$term %in% names(terms_eng), ]
r2 <- res[res$model=="M2 连续评分校正" & res$term %in% names(terms_eng), ]
r2 <- r2[r2$term=="score", ]
plotdf <- rbind(r1, r2)
plotdf$lab <- terms_eng[plotdf$term]
plotdf$Pnum <- as.numeric(plotdf$P)
plotdf <- plotdf[c(1,7,2,3,4,5,6), ]  # 评分两行置前

png(file.path(OUT,"WP7.5_多变量森林图.png"), width=2000, height=1200, res=200)
par(mar=c(4.5,13,3,1))
y <- nrow(plotdf):1
xlo <- max(0.05, min(plotdf$lo)*0.8); xhi <- max(plotdf$hi)*1.2
plot(plotdf$HR, y, xlim=c(xlo,xhi), ylim=c(0.5,nrow(plotdf)+0.6), pch=15, cex=1.4,
     log="x", xlab="Hazard ratio (log scale), adjusted", ylab="", yaxt="n", cex.lab=1.15)
abline(v=1, lty=2, col="grey40")
arrows(plotdf$lo, y, plotdf$hi, y, angle=90, code=3, length=0.07, lwd=1.5, col="#2b6cb0")
axis(2, at=y, labels=paste0(plotdf$lab, "   "), las=2, cex.axis=1.0)
text(plotdf$hi*1.12, y, sprintf("%.2f (%.2f–%.2f), P=%s", plotdf$HR, plotdf$lo, plotdf$hi, plotdf$P),
     adj=c(0,0.5), cex=0.92, xpd=NA)
title(sprintf("Multivariable Cox: OS predictors in TCGA-LUAD (n=%d, %d events)", nrow(dd), sum(dd$OS==1)), cex.main=1.15)
dev.off()
cat("saved WP7.5_多变量森林图.png\n")
