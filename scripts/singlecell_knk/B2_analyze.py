# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
BASE = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
df = pd.read_csv(f"{BASE}/B2_knockout_IGF2BP3_full.csv")
df = df.sort_values('distance', ascending=False).reset_index(drop=True)
df['rnk'] = np.arange(1, len(df) + 1)
AXIS = {'SLC7A11':'ENSG00000151012','SLC3A2':'ENSG00000168003','GCLC':'ENSG00000001084',
        'GCLM':'ENSG00000023909','GOT2':'ENSG00000125166','TFRC':'ENSG00000072274',
        'VDAC2':'ENSG00000165637','IGF2BP2':'ENSG00000073792','IGF2BP3':'ENSG00000136231'}
print('=== 防御轴成员扰动排名（总 n=%d，按 distance 降序）===' % len(df))
for s, e in AXIS.items():
    r = df[df.gene == e]
    if len(r):
        print("%-10s rank=%3d  Z=%6.2f  FC=%8.2f  p.adj=%.3g" %
              (s, int(r['rnk'].iloc[0]), r['Z'].iloc[0], r['FC'].iloc[0], r['p.adj'].iloc[0]))
axis_rank = df[(df.gene.isin(AXIS.values())) & (df.gene != 'ENSG00000136231')]['rnk'].values
others = df[~df.gene.isin(AXIS.values())]['rnk'].values
u, p = mannwhitneyu(axis_rank, others, alternative='less')
print('\n防御轴(8,剔KO)平均排名=%.1f vs 其他=%.1f | Mann-Whitney p=%.3g' % (axis_rank.mean(), others.mean(), p))
mt = {'ENSG00000198899':'MT-ATP6','ENSG00000198886':'MT-ATP6','ENSG00000198712':'MT-CO3',
      'ENSG00000198938':'MT-ND4','ENSG00000198763':'MT-ND4L','ENSG00000198727':'MT-ND5',
      'ENSG00000198888':'MT-ND3','ENSG00000198840':'MT-ND2','ENSG00000198804':'MT-CO1'}
top20 = df.head(20)
print('前20中线粒体基因: %d/20' % sum(1 for g in top20.gene if g in mt))
print('前8基因:', [mt.get(g, g) for g in df.head(8).gene])
df.to_csv(f"{BASE}/B2_knockout_IGF2BP3_ranked.csv", index=False)
print('\n已保存 ranked 版')
