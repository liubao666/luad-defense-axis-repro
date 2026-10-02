# WP2 亚组预后分层：TCGA-LUAD 9 基因防御轴评分 × 5 类亚组 KM
# 口径（稿件 §2.5）：log2(x+1) -> 队列内各基因 z -> 9 基因算术平均；队列内中位数二分
# 输入：Xena HiSeqV2 + Xena survival + GDC 临床 + GDC MC3 突变（WP2_subgroup_survival/ 下）
# 输出：KM 图 5 张、森林图 1 张、结果表 WP2_亚组KM结果.tsv
library(survival)

genes <- c("SLC7A11","SLC3A2","GCLC","GCLM","GOT2","TFRC","VDAC2","IGF2BP2","IGF2BP3")
DIR <- "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP2_subgroup_survival"
PB  <- "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB"

expr <- read.delim(gzfile(file.path(PB,"TCGA-LUAD_HiSeqV2.gz")), check.names=FALSE, row.names=1)
surv <- read.delim(gzfile(file.path(PB,"TCGA-LUAD_survival.tsv.gz")), check.names=FALSE)
clin <- read.delim(file.path(DIR,"clinical_gdc.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
mut  <- read.delim(file.path(DIR,"mutation_status.tsv"),  check.names=FALSE, stringsAsFactors=FALSE)

# ---- 表达矩阵：取 9 基因、原发瘤(01)、12位患者 ----
sym <- sub("\\|.*$","",rownames(expr))
expr <- expr[sym %in% genes, ]; sym <- sym[sym %in% genes]
rownames(expr) <- sym
samples <- colnames(expr)
sid <- substr(samples,1,12)
tum <- substr(samples,14,15)=="01"
sid <- sid[tum]; expr <- expr[,tum]

# ---- 生存：12位去重 ----
ss <- substr(surv$sample,1,12)
surv <- surv[!duplicated(ss), ]; rownames(surv) <- ss[!duplicated(ss)]

common <- intersect(sid, rownames(surv))
m <- match(common, sid)
E <- expr[,m]; S <- surv[common,]

# ---- 9 基因评分 ----
X <- log2(as.matrix(E)+1)
X <- t(scale(t(X)))
score <- colMeans(X[genes,])
grp <- ifelse(score >= median(score, na.rm=TRUE), "high","low")

# ---- 合并临床与突变（突变表只含检出突变的患者，未记录者按 WT 处理）----
rownames(clin) <- clin$patient
all_pat <- common
keap1 <- setNames(mut$KEAP1_mut, mut$patient)
nfe2l2 <- setNames(mut$NFE2L2_mut, mut$patient)
egfrv  <- setNames(mut$EGFR_mut, mut$patient)
kmut <- ifelse(all_pat %in% names(keap1), keap1[all_pat], 0)
nmut <- ifelse(all_pat %in% names(nfe2l2), nfe2l2[all_pat], 0)
emut <- ifelse(all_pat %in% names(egfrv), egfrv[all_pat], 0)

df <- data.frame(
  patient = common,
  OS.time = as.numeric(S$OS.time),
  OS      = as.numeric(S$OS),
  score   = as.numeric(score),
  grp     = factor(grp, levels=c("low","high")),
  age     = as.numeric(clin[common,"age"]),
  ageGrp  = ifelse(is.na(as.numeric(clin[common,"age"])), NA,
             ifelse(as.numeric(clin[common,"age"]) > 66, ">66","<=66")),
  gender  = clin[common,"gender"],
  stage   = clin[common,"ajcc_stage"],
  stageGrp= ifelse(clin[common,"ajcc_stage"] %in% c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"), "I-II",
             ifelse(clin[common,"ajcc_stage"] %in% c("Stage IIIA","Stage IIIB","Stage IV"), "III-IV", NA)),
  knrf2   = ifelse(pmax(kmut, nmut)==1, "mut","WT"),
  egfr    = ifelse(emut==1, "mut","WT")
)

df$OSmo <- df$OS.time / 30.4375
write.table(df, file.path(DIR,"WP2_分析用样本表.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat("分析样本量:", nrow(df), " 死亡事件:", sum(df$OS==1,na.rm=TRUE), "\n")
cat("评分中位数:", round(median(df$score),3), " 高/低:", sum(df$grp=="high"), "/", sum(df$grp=="low"), "\n")
cat("ageGrp NA:", sum(is.na(df$ageGrp)), " gender空:", sum(df$gender==""|is.na(df$gender)),
    " stageGrp NA:", sum(is.na(df$stageGrp)), " knrf2 NA:", sum(is.na(df$knrf2)), " egfr NA:", sum(is.na(df$egfr)), "\n\n")

# ---- 亚组定义 ----
subgroups <- list(
  Age    = list(var="ageGrp", label="年龄",      levels=c("<=66",">66")),
  Gender = list(var="gender", label="性别",      levels=c("female","male")),
  Stage  = list(var="stageGrp", label="病理分期", levels=c("I-II","III-IV")),
  KEAP1_NRF2 = list(var="knrf2", label="KEAP1/NFE2L2 突变", levels=c("WT","mut")),
  EGFR   = list(var="egfr", label="EGFR 突变", levels=c("WT","mut"))
)

km_one <- function(dd, title){
  fit <- survfit(Surv(OSmo, OS) ~ grp, data=dd)
  lr  <- survdiff(Surv(OSmo, OS) ~ grp, data=dd)
  p   <- 1 - pchisq(lr$chisq, 1)
  hr  <- summary(coxph(Surv(OSmo, OS) ~ grp, data=dd))$conf.int
  list(fit=fit, p=p, HR=hr[1,1], lo=hr[1,3], hi=hr[1,4],
       n=nrow(dd), ev=sum(dd$OS==1),
       n_high=sum(dd$grp=="high"), ev_high=sum(dd$OS[dd$grp=="high"]==1),
       n_low=sum(dd$grp=="low"),  ev_low=sum(dd$OS[dd$grp=="low"]==1))
}

res <- list()
# 全队列参考行（用于与稿件 §3.6 口径核对方向与量级）
dd_all <- df[!is.na(df$OS), ]
r_all <- km_one(dd_all, "All")
r_all$subgroup <- "全队列(参考)"; r_all$level <- "全部"
res[["ALL"]] <- r_all
for(sg in names(subgroups)){
  info <- subgroups[[sg]]
  for(lv in info$levels){
    dd <- df[df[[info$var]]==lv & !is.na(df[[info$var]]) & !is.na(df$OS), ]
    if(nrow(dd) < 20 || sum(dd$OS==1) < 5){ res[[paste(sg,lv)]] <- NULL; next }
    r <- km_one(dd, paste(info$label, lv))
    r$subgroup <- info$label; r$level <- lv
    res[[paste(sg,lv)]] <- r
  }
}

cat(sprintf("%-22s %5s %5s %6s %6s %6s %10s\n","亚组","n","事件","HR","lo95","hi95","logrankP"))
out <- data.frame()
for(k in names(res)){
  r <- res[[k]]
  cat(sprintf("%-22s %5d %5d %6.2f %6.2f %6.2f %10.3g\n",
              k, r$n, r$ev, r$HR, r$lo, r$hi, r$p))
  out <- rbind(out, data.frame(subgroup=r$subgroup, level=r$level, n=r$n, events=r$ev,
    n_high=r$n_high, events_high=r$ev_high, n_low=r$n_low, events_low=r$ev_low,
    HR=r$HR, lo=r$lo, hi=r$hi, logrank_P=r$p))
}
write.table(out, file.path(DIR,"WP2_亚组KM结果.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ---- 作图：每亚组一张 2 面板 KM（英文标注，投稿口径）----
cols <- c("low"="#2b6cb0","high"="#c53030")
eng_label <- c("Age"="Age", "Gender"="Sex", "Stage"="Pathological stage",
               "KEAP1_NRF2"="KEAP1/NFE2L2", "EGFR"="EGFR")
eng_level <- c("<=66"="≤66", ">66"=">66", "female"="female", "male"="male",
               "I-II"="Stage I–II", "III-IV"="Stage III–IV", "WT"="WT", "mut"="mutant")
for(sg in names(subgroups)){
  info <- subgroups[[sg]]
  png(file.path(DIR, sprintf("WP2_KM_%s.png", sg)), width=2200, height=1100, res=200)
  par(mfrow=c(1,2), mar=c(4.5,4.5,3,1), cex.lab=1.15, cex.axis=1.05, cex.main=1.2)
  for(i in seq_along(info$levels)){
    lv <- info$levels[i]
    r <- res[[paste(sg,lv)]]
    if(is.null(r)){ plot.new(); text(0.5,0.5,"insufficient data"); next }
    plot(r$fit, col=cols, lwd=2.2, mark.time=TRUE,
         xlab="Overall survival (months)", ylab="Survival probability",
         main=sprintf("%s: %s", eng_label[sg], eng_level[lv]),
         cex.main=1.15)
    legend("topright", bty="n", lwd=2.2, col=cols,
           legend=c(sprintf("Low score (n=%d, events=%d)", r$n_low, r$ev_low),
                    sprintf("High score (n=%d, events=%d)", r$n_high, r$ev_high)),
           cex=0.95)
    txt <- sprintf("HR = %.2f (95%% CI %.2f–%.2f)\nlog-rank P = %.3g", r$HR, r$lo, r$hi, r$p)
    text(par("usr")[1] + 0.04*diff(par("usr")[1:2]),
         grconvertY(0.20,"ndc","user"), txt, cex=1.0, adj=c(0,0))
  }
  dev.off()
  cat("saved WP2_KM_", sg, ".png\n", sep="")
}

# ---- 森林图：亚组单元 HR ----
png(file.path(DIR,"WP2_forest.png"), width=1800, height=1200, res=200)
par(mar=c(4,12,3,2))
o <- out[out$subgroup!="全队列(参考)", ]
o <- o[order(o$HR), ]
y <- nrow(o):1
xlo <- max(0.1, min(o$lo)*0.85); xhi <- max(o$hi)*1.15
plot(o$HR, y, xlim=c(xlo, xhi), ylim=c(0.5,nrow(o)+0.5), pch=15, cex=1.3,
     xlab="HR (high vs low score), log scale", ylab="", yaxt="n", log="x", cex.lab=1.1)
abline(v=1, lty=2, col="grey40")
arrows(o$lo, y, o$hi, y, angle=90, code=3, length=0.06, lwd=1.4, col="#2b6cb0")
sg_eng <- c("年龄"="Age", "性别"="Sex", "病理分期"="Pathological stage",
            "KEAP1/NFE2L2 突变"="KEAP1/NFE2L2", "EGFR 突变"="EGFR")
axis(2, at=y, labels=paste(sg_eng[o$subgroup], eng_level[o$level], sep=" | "), las=2, cex.axis=0.9)
title("OS hazard of high 9-gene defense-axis score across TCGA-LUAD subgroups", cex.main=1.1)
dev.off()
cat("saved WP2_forest.png\n")
