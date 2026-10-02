# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
BASE = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
AXIS = {'SLC7A11':'ENSG00000151012','SLC3A2':'ENSG00000168003','GCLC':'ENSG00000001084',
        'GCLM':'ENSG00000023909','GOT2':'ENSG00000125166','TFRC':'ENSG00000072274',
        'VDAC2':'ENSG00000165637','IGF2BP2':'ENSG00000073792','IGF2BP3':'ENSG00000136231'}
SYM = {v:k for k,v in AXIS.items()}

def report(path, tag):
    df = pd.read_csv(path).sort_values('distance', ascending=False).reset_index(drop=True)
    df['rnk'] = np.arange(1, len(df) + 1)
    df['sym'] = df['gene'].map(SYM).fillna(df['gene'])
    print('==== %s (n=%d) ====' % (tag, len(df)))
    for s, e in AXIS.items():
        r = df[df.gene == e]
        if len(r):
            print("%-10s rank=%3d/%d  Z=%6.2f  FC=%9.4f  p.adj=%.3g" %
                  (s, int(r['rnk'].iloc[0]), len(df), r['Z'].iloc[0], r['FC'].iloc[0], r['p.adj'].iloc[0]))
    ar = df[(df.gene.isin(AXIS.values())) & (df.gene != 'ENSG00000136231')]['rnk'].values
    ot = df[~df.gene.isin(AXIS.values())]['rnk'].values
    u, p = mannwhitneyu(ar, ot, alternative='less')
    print('防御轴(8,剔KO)平均排名=%.1f vs 其他=%.1f | Mann-Whitney p=%.3g' % (ar.mean(), ot.mean(), p))
    print('前10:', df.head(10)['sym'].tolist())
    print('FC 非零的防御轴成员数: %d/8' % sum(
        1 for e in list(AXIS.values()) if e != 'ENSG00000136231'
        and abs(df[df.gene == e]['FC'].iloc[0]) > 0))
    print()
    return df

df = report(f"{BASE}/B2_knockout_IGF2BP3_noMT_full.csv", "剔MT·全肿瘤细胞·虚拟敲除")
df.to_csv(f"{BASE}/B2_knockout_IGF2BP3_noMT_ranked.csv", index=False)
