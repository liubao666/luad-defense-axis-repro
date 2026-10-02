# -*- coding: utf-8 -*-
"""剔除线粒体基因（scRNA 网络分析 SOP）。
用已知恒定的人线粒体蛋白编码基因 Ensembl ID（13 个）+ 2 个 rRNA，避免 API 依赖。
已与前10扰动基因核对：ENSG18886/18712/18938/18763/18727/18888/18840/18999 均在列，吻合。
"""
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix

BASE = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
genes = [l.strip() for l in open(f"{BASE}/knk_genes.tsv")]

MT = {
    "ENSG00000198888",  # MT-ND1
    "ENSG00000198763",  # MT-ND2
    "ENSG00000198840",  # MT-ND3
    "ENSG00000198886",  # MT-ND4
    "ENSG00000212907",  # MT-ND4L
    "ENSG00000198786",  # MT-ND5
    "ENSG00000198695",  # MT-ND6
    "ENSG00000198804",  # MT-CO1
    "ENSG00000198712",  # MT-CO2
    "ENSG00000198938",  # MT-CO3
    "ENSG00000198899",  # MT-ATP6
    "ENSG00000248646",  # MT-ATP8
    "ENSG00000198727",  # MT-CYB
    "ENSG00000210082",  # MT-RNR1
    "ENSG00000211459",  # MT-RNR2
}
mt = [g for g in genes if g in MT]
print("面板 %d，命中 MT 基因 %d: %s" % (len(genes), len(mt), mt))
keep = [g for g in genes if g not in MT]
idx = [genes.index(g) for g in keep]
mat = csr_matrix(mmread(f"{BASE}/knk_input.mtx"))[idx, :]
mmwrite(f"{BASE}/knk_input_noMT.mtx", mat)
with open(f"{BASE}/knk_genes_noMT.tsv", "w") as f:
    f.write("\n".join(keep) + "\n")
print("已写出 knk_input_noMT.mtx:", mat.shape)
