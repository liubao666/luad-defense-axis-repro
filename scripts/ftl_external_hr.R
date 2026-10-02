## 复审补充：FTL 在 GSE68465 / GSE50081 的外部队列单变量 HR（OS，每标准差）
.libPaths(c("/Users/liubao/WorkBuddy/2026-09-23-10-23-05/Rlib", .libPaths()))
suppressMessages(library(survival))

BASE <- "/Users/liubao/WorkBuddy/2026-09-23-10-23-05/生信学习计划"
OUT  <- file.path(BASE, "第16课数据/extval")
D8   <- file.path(BASE, "第8课数据")
source(file.path(OUT, "ext_geo_common.R"))

ann <- read.csv(file.path(D8, "GPL570_探针对应基因.csv"), stringsAsFactors = FALSE)
P2G <- setNames(ann$gene, ann$probe)

num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "", as.character(x))))

run_ftl <- function(gse, ad_col, ad_val, os_t, os_e) {
  gm <- parse_geo_meta(file.path(OUT, paste0(gse, "_series_matrix.txt.gz")))
  m <- gm$meta
  E <- collapse_probe(parse_geo_expr(gm$lines, gm$table_line), P2G)
  is_ad <- m[[ad_col]] %in% ad_val
  t_os <- num(m[[os_t]])
  e_os <- as.integer(tolower(as.character(m[[os_e]])) %in% c("dead","deceased","1"))
  scale_t <- if (max(t_os, na.rm = TRUE) < 40) 12 else 1
  t_os <- t_os * scale_t
  keep <- is_ad & is.finite(t_os) & t_os > 0
  s <- m$gsm[keep]
  common <- intersect(s, colnames(E))
  stopifnot("FTL" %in% rownames(E))
  x <- as.numeric(E["FTL", common])
  z <- as.numeric(scale(x))
  tt <- t_os[match(common, m$gsm)]
  ee <- e_os[match(common, m$gsm)]
  cx <- coxph(Surv(tt, ee) ~ z)
  s <- summary(cx)
  data.frame(gse = gse, n = sum(keep), events = sum(ee),
             HR_perSD = round(s$conf.int[1, "exp(coef)"], 3),
             lo = round(s$conf.int[1, "lower .95"], 3),
             hi = round(s$conf.int[1, "upper .95"], 3),
             p = formatC(s$coefficients[1, "Pr(>|z|)"], format = "g", digits = 3))
}

r1 <- run_ftl("GSE50081", "histology", "adenocarcinoma",
              "survival time", "status")
r2 <- run_ftl("GSE68465", "disease_state", "Lung Adenocarcinoma",
              "months_to_last_contact_or_death", "vital_status")
res <- rbind(r1, r2)
print(res)
write.csv(res, file.path(OUT, "FTL_外部队列单变量HR.csv"), row.names = FALSE)
cat("DONE\n")
