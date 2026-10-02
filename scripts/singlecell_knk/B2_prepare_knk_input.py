# -*- coding: utf-8 -*-
"""模块二数据准备：E-MTAB-6653 (Lambrechts 2018 人 LUAD scRNA)
从全矩阵筛选肿瘤细胞子集 -> 基因过滤 -> 输出 scTenifoldKnk 输入矩阵 (genes x cells, mtx)。
"""
import gzip
import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix, csc_matrix
import pandas as pd

BASE = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
RNG = np.random.default_rng(20260930)

DEFENSE_ENSG = {
    "SLC7A11": "ENSG00000151012",
    "SLC3A2": "ENSG00000168003",
    "GCLC": "ENSG00000001084",
    "GCLM": "ENSG00000023909",
    "GOT2": "ENSG00000125166",
    "TFRC": "ENSG00000072274",
    "VDAC2": "ENSG00000165637",
    "IGF2BP2": "ENSG00000073792",
    "IGF2BP3": "ENSG00000136231",
}
READERS_ENSG = {
    "IGF2BP1": "ENSG00000153875",
    "YTHDF1": "ENSG00000149658",
    "YTHDF2": "ENSG00000198492",
    "YTHDF3": "ENSG00000185728",
    "YTHDC1": "ENSG00000163439",
    "YTHDC2": "ENSG00000054054",
    "HNRNPA2B1": "ENSG00000122566",
    "HNRNPC": "ENSG00000092199",
    "FMR1": "ENSG00000102081",
    "RBMX": "ENSG00000138169",
}

print("[1/5] 读取基因/细胞注释 ...")
genes = [l.split("\t")[0] for l in gzip.open(f"{BASE}/LUAD6653_rows.gz", "rt")]
cells = [l.strip() for l in gzip.open(f"{BASE}/LUAD6653_cols.gz", "rt")]
meta = pd.read_csv(f"{BASE}/LUAD6653_metadata.tsv", sep="\t")
meta = meta.loc[:, ~meta.columns.duplicated()].set_index("id")
print(f"  基因 {len(genes)} | 细胞 {len(cells)} | metadata {meta.shape[0]}")

# 肿瘤细胞：sampling_site 含 tumour
tum_mask = meta.loc[cells, "sampling_site"].str.match(r"^tumour").values
tum_cells = [c for c, m in zip(cells, tum_mask) if m]
print(f"[2/5] 肿瘤细胞 {len(tum_cells)} / {len(cells)}")

# 随机抽取 15000 个肿瘤细胞（scTenifoldKnk 内部子采样 500x10，输入过多无额外收益）
if len(tum_cells) > 15000:
    tum_cells = list(RNG.choice(np.array(tum_cells), 15000, replace=False))
print(f"  抽样后 {len(tum_cells)}")

print("[3/5] 读取表达矩阵 (63M nnz, 约 30-60s) ...")
mat = mmread(gzip.open(f"{BASE}/LUAD6653_counts.mtx.gz", "rt"))  # genes x cells, coo
mat = csc_matrix(mat)
cell_idx = [cells.index(c) for c in tum_cells]
sub = mat[:, cell_idx].tocsr()  # genes x selected cells
print(f"  子矩阵 {sub.shape[0]} x {sub.shape[1]}, nnz={sub.nnz}")

print("[4/5] 基因过滤 ...")
detected = (sub > 0)
det_rate = np.asarray(detected.mean(axis=1)).ravel()  # 每个基因被检出的细胞比例
keep = det_rate >= 0.02
total = np.asarray(sub.sum(axis=1)).ravel()
# 保证防御轴 + 阅读器 + KO 基因不被过滤掉
# 铁死亡核心调控子（教科书级公认分子，用于聚焦面板生物学相关性）
FerroptosisCore = {
    "GPX4": "ENSG00000167468", "ACSL4": "ENSG00000123983", "NFE2L2": "ENSG00000116044",
    "KEAP1": "ENSG00000079999", "HMOX1": "ENSG00000100292", "NCOA4": "ENSG00000126368",
    "ALOX15": "ENSG00000161905", "AIFM2": "ENSG00000042286", "GCH1": "ENSG00000180176",
    "HSPB1": "ENSG00000106211", "LPCAT3": "ENSG00000111684", "CBS": "ENSG00000160200",
    "PRNP": "ENSG00000171867", "FANCD2": "ENSG00000144554", "HSPA5": "ENSG00000044574",
}
force = set(DEFENSE_ENSG.values()) | set(READERS_ENSG.values()) | set(FerroptosisCore.values())
gi = {g: i for i, g in enumerate(genes)}
force_idx = {gi[e] for e in force if e in gi}
keep[list(force_idx)] = True
n_pass = int(keep.sum())
print(f"  检出率>=2%: {n_pass} 基因")

# 按总表达量补足到 TOPN（强制保留基因除外）。
# 关键工程约束：scTenifoldKnk 的 pcNet 是「每基因一次 SVD」，实测 500 基因=12.6s/网络，
# 3000 基因×10 网络 >2 小时且内存溢出。聚焦面板 ~800 基因可在数分钟跑完，
# 且对「IGF2BP3 敲除是否扰动防御轴」这一具体问题，通路的聚焦面板统计功效反而更高。
TOPN = 800
if n_pass > TOPN:
    score = total.copy()
    score[~keep] = -1
    thr = np.sort(score[keep])[::-1][TOPN - 1]
    keep2 = keep & ((score >= thr) | np.isin(np.arange(len(genes)), list(force_idx)))
else:
    keep2 = keep
idx = np.where(keep2)[0]
sub2 = sub[idx, :]
print(f"  最终 {sub2.shape[0]} 基因 x {sub2.shape[1]} 细胞")

# 防御轴基因存在性核验
present = {s: (e in gi and keep2[gi[e]]) for s, e in DEFENSE_ENSG.items()}
print(f"  防御轴基因保留情况: {present}")

print("[5/5] 写出 knk_input.mtx (genes x cells) ...")
mmwrite(f"{BASE}/knk_input.mtx", csr_matrix(sub2).astype(float))
with open(f"{BASE}/knk_genes.tsv", "w") as f:
    f.write("\n".join(genes[i] for i in idx) + "\n")
with open(f"{BASE}/knk_cells.tsv", "w") as f:
    f.write("\n".join(tum_cells) + "\n")
# 防御轴基因 symbol->ENSG 对照（供 R 与作图使用）
with open(f"{BASE}/knk_axis_map.tsv", "w") as f:
    f.write("symbol\tensg\tretained\n")
    for s, e in DEFENSE_ENSG.items():
        f.write(f"{s}\t{e}\t{present[s]}\n")
print("完成。")
