# -*- coding: utf-8 -*-
"""路径 B 模块一：NSCLC Visium 空间转录组 —— 防御轴/m6A 阅读器的空间定位
数据：E-MTAB-13530（ArrayExpress，4 例 NSCLC Visium 样本，10x Space Ranger 输出）
产出：
  Fig_B1_防御轴评分空间图.png   （4 例，防御轴 9 基因 spot 得分）
  Fig_B2_阅读器评分空间图.png   （4 例，m6A 阅读器 12 基因 spot 得分）
  Fig_B3_肿瘤区vs基质区箱线图.png
  Fig_B4_spot级防御-阅读器相关.png
  B1_spot_scores.tsv            （每 spot 得分表，供正文补表）
"""
import json, os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sc.settings.verbosity = 0
SAMPLES = ['P24_T1', 'P17_T1', 'P10_T1', 'P11_T1']
DATA = os.path.dirname(os.path.abspath(__file__))

DEFENSE = ['SLC7A11', 'SLC3A2', 'GCLC', 'GCLM', 'GOT2', 'TFRC', 'VDAC2', 'IGF2BP2', 'IGF2BP3']
READERS = ['IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'YTHDF1', 'YTHDF2', 'YTHDF3',
           'YTHDC1', 'YTHDC2', 'HNRNPA2B1', 'HNRNPC', 'FMR1', 'RBMX']
EPITHELIAL = ['EPCAM', 'KRT18', 'KRT19']

adatas = {}
for s in SAMPLES:
    a = sc.read_10x_h5(f'{DATA}/{s}.h5', genome=None)
    a.var_names_make_unique()
    pos = pd.read_csv(f'{DATA}/{s}_spatial/tissue_positions_list.csv', header=None,
                      names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row', 'pxl_col'])
    pos = pos.set_index('barcode')
    # 统一 barcode：两端都去掉 -1 后缀再匹配
    pos.index = pos.index.str.replace('-1$', '', regex=True)
    keep = [b for b in a.obs_names if b.replace('-1', '') in pos.index]
    a = a[keep].copy()
    pos = pos.loc[[b.replace('-1', '') for b in keep]]
    a.obs['in_tissue'] = pos['in_tissue'].values
    a.obsm['spatial'] = pos[['pxl_col', 'pxl_row']].values
    a.obs['sample'] = s
    # QC
    sc.pp.filter_cells(a, min_genes=200)
    a = a[a.obs['in_tissue'] == 1].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    adatas[s] = a
    print(f'{s}: {a.n_obs} spots')

all_scores = []
for s, a in adatas.items():
    present_def = [g for g in DEFENSE if g in a.var_names]
    present_read = [g for g in READERS if g in a.var_names]
    sc.tl.score_genes(a, present_def, score_name='defense_score')
    sc.tl.score_genes(a, present_read, score_name='reader_score')
    sc.tl.score_genes(a, [g for g in EPITHELIAL if g in a.var_names], score_name='epi_score')
    df = a.obs[['sample', 'defense_score', 'reader_score', 'epi_score']].copy()
    df['x'] = a.obsm['spatial'][:, 0]
    df['y'] = a.obsm['spatial'][:, 1]
    all_scores.append(df)
scores = pd.concat(all_scores)
scores.to_csv(f'{DATA}/B1_spot_scores.tsv', sep='\t')

# 肿瘤区定义：epi_score > 0（EPCAM/KRT18/KRT19 表达高于均值）
scores['region'] = np.where(scores['epi_score'] > 0, 'tumor-rich', 'stroma-rich')

def spatial_panel(ax, df, col, title, cmap='RdYlBu_r'):
    v = df[col].values
    sc_map = ax.scatter(df['x'], -df['y'], c=v, cmap=cmap, s=6, edgecolors='none')
    ax.set_title(title, fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect('equal')
    plt.colorbar(sc_map, ax=ax, fraction=0.04, pad=0.02)

for col, name, ttl in [('defense_score', 'Fig_B1_防御轴评分空间图', 'Defence-axis score'),
                       ('reader_score', 'Fig_B2_阅读器评分空间图', 'm6A-reader score')]:
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    for ax, s in zip(axes, SAMPLES):
        sub = scores[scores['sample'] == s]
        spatial_panel(ax, sub, col, f'{s}  {ttl}')
    fig.suptitle(ttl + ' — NSCLC Visium (E-MTAB-13530)', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f'{DATA}/{name}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

# 箱线图：肿瘤区 vs 基质区
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, col, ttl in zip(axes, ['defense_score', 'reader_score'],
                        ['Defence-axis score', 'm6A-reader score']):
    data = [scores.loc[scores['region'] == r, col].values for r in ['tumor-rich', 'stroma-rich']]
    bp = ax.boxplot(data, tick_labels=['tumor-rich\n(EPCAM/KRT18/KRT19 high)', 'stroma-rich'], showfliers=False)
    from scipy import stats
    p = stats.mannwhitneyu(data[0], data[1]).pvalue
    ax.set_title(f'{ttl}\nMann–Whitney P = {p:.2e}')
    ax.set_ylabel('spot score')
fig.suptitle('Defense program is tumor-region-enriched in LUAD space', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(f'{DATA}/Fig_B3_肿瘤区vs基质区箱线图.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# spot 级相关
fig, ax = plt.subplots(figsize=(5.5, 5))
from scipy import stats
rho, p = stats.spearmanr(scores['defense_score'], scores['reader_score'])
ax.scatter(scores['defense_score'], scores['reader_score'], s=4, alpha=0.3, edgecolors='none')
ax.set_xlabel('Defence-axis score'); ax.set_ylabel('m6A-reader score')
ax.set_title(f'Per-spot Spearman ρ = {rho:.3f}, P = {p:.2e}\n(n = {len(scores)} spots, 4 tumors)')
fig.tight_layout()
fig.savefig(f'{DATA}/Fig_B4_spot级防御-阅读器相关.png', dpi=300, bbox_inches='tight')
plt.close(fig)

print('\n=== 汇总 ===')
print('spots 总数:', len(scores))
print(scores.groupby('region')[['defense_score', 'reader_score']].median())
print('spot 级 Spearman ρ =', round(rho, 3), 'P =', p)
