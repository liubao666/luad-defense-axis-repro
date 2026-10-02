# -*- coding: utf-8 -*-
"""对症检验：转录后调控框架下，IGF2BP3 表达与防御轴基因的单细胞共变（Spearman）。
在 IGF2BP3+ 肿瘤细胞内计算（该调控子实际表达的细胞）。"""
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
import pandas as pd

BASE = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
genes = [l.strip() for l in open(f"{BASE}/knk_genes.tsv")]
gidx = {g: i for i, g in enumerate(genes)}
AXIS = {'SLC7A11':'ENSG00000151012','SLC3A2':'ENSG00000168003','GCLC':'ENSG00000001084',
        'GCLM':'ENSG00000023909','GOT2':'ENSG00000125166','TFRC':'ENSG00000072274',
        'VDAC2':'ENSG00000165637','IGF2BP2':'ENSG00000073792'}

mat = csr_matrix(mmread(f"{BASE}/knk_input.mtx")).toarray()  # 824 x 15000
ko = mat[gidx['ENSG00000136231'], :]
pos = ko > 0
print("IGF2BP3+ 肿瘤细胞 n=%d / %d" % (pos.sum(), len(ko)))
sub = mat[:, pos]  # 824 x n_pos

rows = []
ko_expr = sub[gidx['ENSG00000136231'], :]
for s, e in AXIS.items():
    x = sub[gidx[e], :].astype(float)
    if np.std(x) == 0:
        rows.append((s, np.nan, np.nan, int((x > 0).sum())))
        continue
    rho, p = spearmanr(ko_expr, x)
    rows.append((s, rho, p, int((x > 0).sum())))
df = pd.DataFrame(rows, columns=['gene', 'spearman_rho', 'p_value', 'n_expressing'])
df['p_adj'] = np.minimum(df['p_value'] * len(df), 1.0)
print(df.to_string(index=False))
sig = df[df['p_adj'] < 0.05]
print("\nFDR<0.05 的防御轴成员: %d/8" % len(sig))
print(sig[['gene','spearman_rho','p_adj']].to_string(index=False))
df.to_csv(f"{BASE}/B2_IGF2BP3_防御轴_单细胞共变.tsv", sep="\t", index=False)

# 对照：随机 8 个同表达量基因的平均相关（判断 8 个基因整体是否正相关偏多）
rng = np.random.default_rng(1)
det = (mat > 0).mean(axis=1)
ko_det = (ko > 0).mean()
cand = [i for i in range(len(genes)) if abs(det[i] - ko_det) < 0.02 and genes[i] != 'ENSG00000136231']
rand_rho = []
for _ in range(200):
    pick = rng.choice(cand, 8, replace=False)
    rs = [spearmanr(ko_expr, mat[i, pos].astype(float))[0] for i in pick]
    rand_rho.append(np.nanmean(rs))
print("\n随机8基因(表达量匹配)平均rho: 中位=%.3f, IGF2BP3 vs 防御轴平均rho=%.3f" %
      (np.nanmedian(rand_rho), df['spearman_rho'].mean()))
