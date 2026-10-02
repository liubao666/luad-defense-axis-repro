# -*- coding: utf-8 -*-
"""模块二补充：只在 IGF2BP3+ 肿瘤细胞内做虚拟敲除（生物学正确的问题）
全局网络里 IGF2BP3 仅 ~10/500 细胞表达，效应被稀释；在其表达的细胞子集内重跑。
"""
import gzip
import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix, csc_matrix
import pandas as pd

BASE = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
RNG = np.random.default_rng(20260930)

# 复用主脚本的面板基因（knk_genes.tsv）与细胞（knk_cells.tsv）
genes = [l.strip() for l in open(f"{BASE}/knk_genes.tsv")]
cells = [l.strip() for l in open(f"{BASE}/knk_cells.tsv")]
meta = pd.read_csv(f"{BASE}/LUAD6653_metadata.tsv", sep="\t")
meta = meta.loc[:, ~meta.columns.duplicated()].set_index("id")

# 读取完整子矩阵（824 基因 x 15000 肿瘤细胞）
mat = mmread(f"{BASE}/knk_input.mtx")
mat = csr_matrix(mat).tocsr()  # 824 x 15000, 行=基因
gidx = {g: i for i, g in enumerate(genes)}
ko_row = mat[gidx["ENSG00000136231"], :].toarray().ravel()
pos_mask = ko_row > 0
pos_cells = [c for c, m in zip(cells, pos_mask) if m]
print("IGF2BP3+ 肿瘤细胞:", len(pos_cells), "/", len(cells))
# 该子集内防御轴检出率
AXIS = {"SLC7A11":"ENSG00000151012","SLC3A2":"ENSG00000168003","GCLC":"ENSG00000001084",
        "GCLM":"ENSG00000023909","GOT2":"ENSG00000125166","TFRC":"ENSG00000072274",
        "VDAC2":"ENSG00000165637","IGF2BP2":"ENSG00000073792","IGF2BP3":"ENSG00000136231"}
for s, e in AXIS.items():
    r = mat[gidx[e], pos_mask].toarray().ravel()
    print("  %s 检出率 %.1f%%" % (s, 100 * (r > 0).mean()))

sub = mat[:, pos_mask]
mmwrite(f"{BASE}/knk_input_ko_cells.mtx", sub)
with open(f"{BASE}/knk_cells_ko.tsv", "w") as f:
    f.write("\n".join(pos_cells) + "\n")
print("已写出 knk_input_ko_cells.mtx:", sub.shape)
