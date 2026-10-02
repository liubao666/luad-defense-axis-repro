# -*- coding: utf-8 -*-
"""Fig_B5：IGF2BP3 与防御轴的单细胞表达共变（转录后框架的正确检验）
左：8 个防御轴成员的 Spearman rho 条形图（按显著性着色）
右：代表性散点 IGF2BP3 vs GCLC（rho 最高之一）
"""
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.sans-serif'] = ['Arial Unicode MS', 'PingFang SC', 'Heiti SC', 'sans-serif']
rcParams['axes.unicode_minus'] = False

BASE = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
OUT = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/spatial_data"

genes = [l.strip() for l in open(f"{BASE}/knk_genes.tsv")]
gidx = {g: i for i, g in enumerate(genes)}
AXIS = [('SLC7A11','ENSG00000151012'),('SLC3A2','ENSG00000168003'),('GCLC','ENSG00000001084'),
        ('GCLM','ENSG00000023909'),('GOT2','ENSG00000125166'),('TFRC','ENSG00000072274'),
        ('VDAC2','ENSG00000165637'),('IGF2BP2','ENSG00000073792')]

mat = csr_matrix(mmread(f"{BASE}/knk_input.mtx")).toarray()
ko = mat[gidx['ENSG00000136231'], :]
pos = ko > 0
sub = mat[:, pos]
ko_expr = sub[gidx['ENSG00000136231'], :].astype(float)

res = []
for s, e in AXIS:
    x = sub[gidx[e], :].astype(float)
    rho, p = spearmanr(ko_expr, x)
    padj = min(p * len(AXIS), 1.0)
    res.append((s, rho, padj, x))
df = pd.DataFrame([(r[0], r[1], r[2]) for r in res], columns=['gene','rho','padj'])

plt.rcParams.update({'axes.labelsize': 13, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 11})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.6), gridspec_kw={'width_ratios':[1.5,1]})

# 左：条形图
dd = df.sort_values('rho')
colors = ['#c0392b' if pa < 0.05 else '#95a5a6' for pa in dd['padj']]
bars = ax1.barh(dd['gene'], dd['rho'], color=colors, edgecolor='#7f8c8d')
for b, (rho, padj) in zip(bars, zip(dd['rho'], dd['padj'])):
    star = '***' if padj < 1e-10 else ('**' if padj < 1e-5 else ('*' if padj < 0.05 else 'ns'))
    ax1.text(b.get_width()+0.008, b.get_y()+b.get_height()/2, 'ρ=%.2f %s' % (rho, star),
             va='center', fontsize=13)
ax1.set_xlabel('Spearman correlation with IGF2BP3 (in IGF2BP3+ tumour cells, n=%d)' % pos.sum())
ax1.set_title('Defence-axis members co-vary with IGF2BP3\nred = FDR<0.05', fontsize=13)
ax1.axvline(0, color='#555', lw=0.8)
ax1.set_xlim(0, max(dd['rho'])*1.28)
from matplotlib.patches import Patch
ax1.legend(handles=[Patch(color='#c0392b', label='FDR<0.05 (7/8)'), Patch(color='#95a5a6', label='n.s. (IGF2BP2)')],
           loc='lower right', fontsize=11)

# 右：代表性散点 IGF2BP3 vs GCLC
gi = [i for i,(s,e) in enumerate(AXIS) if s=='GCLC'][0]
x = res[gi][3]
# log1p 变换改善可视化
lx, ly = np.log1p(x), np.log1p(ko_expr)
rho_g, p_g = spearmanr(ly, lx)
ax2.scatter(lx, ly, s=14, alpha=0.6, color='#2980b9', edgecolor='none')
# 拟合线
m, b0 = np.polyfit(lx, ly, 1)
xs = np.linspace(lx.min(), lx.max(), 50)
ax2.plot(xs, m*xs+b0, color='#c0392b', lw=1.8)
ax2.set_xlabel('GCLC expression (log1p CPM)')
ax2.set_ylabel('IGF2BP3 expression (log1p CPM)')
ax2.set_title('Representative: IGF2BP3 vs GCLC\nρ=%.2f, FDR=%.2e' % (rho_g, min(p_g*len(AXIS),1)), fontsize=13)

plt.tight_layout()
plt.savefig(f"{OUT}/Fig_B5_IGF2BP3_防御轴单细胞共变.png", dpi=300, bbox_inches='tight')
print("Fig_B5 已输出")
print(df.to_string(index=False))
