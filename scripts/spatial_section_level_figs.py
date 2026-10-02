# -*- coding: utf-8 -*-
"""Major5 图内修正：删除 spot 级推断 P 值，改标切片级统计（n=4 sections）。
重出 Fig_B3（肿瘤区vs基质区）与 Fig_B4（防御-阅读器相关），spot 级仅作描述性展示。"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

DATA = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/spatial_data"
scores = pd.read_csv(f"{DATA}/B1_spot_scores.tsv", sep="\t").rename(columns={"\xef\xbb\xbfspot": "spot"})
if "spot" not in scores.columns:
    scores = scores.rename(columns={scores.columns[0]: "spot"})
scores["region"] = np.where(scores["epi_score"] > 0, "tumor-rich", "stroma-rich")
SAMPLES = ["P24_T1", "P17_T1", "P10_T1", "P11_T1"]

# 切片级统计
per = []
for s, g in scores.groupby("sample"):
    t = g.loc[g.region == "tumor-rich"]; m = g.loc[g.region == "stroma-rich"]
    per.append((s, t.defense_score.median() - m.defense_score.median(),
                t.reader_score.median() - m.reader_score.median(),
                stats.spearmanr(g.defense_score, g.reader_score).statistic))
per = pd.DataFrame(per, columns=["section", "def_diff", "read_diff", "rho"])
w_def, p_def = stats.wilcoxon(per.def_diff, alternative="greater")
w_read, p_read = stats.wilcoxon(per.read_diff, alternative="greater")
z = np.arctanh(per.rho.values)
t_z, p_z = stats.ttest_1samp(z, 0)

# ---- Fig_B3：箱线图保留（描述性），标题只放切片级统计 ----
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
for ax, col, ttl, pval in zip(axes, ["defense_score", "reader_score"],
                              ["Defence-axis score", "m6A-reader score"], [p_def, p_read]):
    data = [scores.loc[scores.region == r, col].values for r in ["tumor-rich", "stroma-rich"]]
    ax.boxplot(data, tick_labels=["tumor-rich\n(EPCAM/KRT18/KRT19 high)", "stroma-rich"], showfliers=False)
    ax.set_title(f"{ttl}\nsection-level Wilcoxon P = {pval:.3f} (n = 4 sections;\ntumour-minus-stroma median positive in 4/4)", fontsize=9)
    ax.set_ylabel("spot score")
fig.suptitle("Defense program is directionally tumor-region-enriched in NSCLC sections\n"
             "(spot-level boxplots are descriptive only: spots within a section are not independent)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(f"{DATA}/Fig_B3_肿瘤区vs基质区箱线图.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---- Fig_B4：散点保留（描述性），只标每例 ρ 区间与切片级检验 ----
fig, ax = plt.subplots(figsize=(5.5, 5.2))
ax.scatter(scores.defense_score, scores.reader_score, s=4, alpha=0.25, edgecolors="none")
ax.set_xlabel("Defence-axis score"); ax.set_ylabel("m6A-reader score")
ax.set_title(f"per-section Spearman ρ = {per.rho.min():.3f}–{per.rho.max():.3f} (mean {per.rho.mean():.2f})\n"
             f"section-level test P = {p_z:.3f} (n = 4 sections)\n"
             f"spot display is descriptive only — spots are not independent", fontsize=9)
fig.tight_layout()
fig.savefig(f"{DATA}/Fig_B4_spot级防御-阅读器相关.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("重出完成:", f"P_def={p_def:.4f}", f"P_read={p_read:.4f}", f"P_z={p_z:.3f}")
