# -*- coding: utf-8 -*-
"""WP1 空间第二批遗留分析
A. 4 例 E-MTAB-13530（De Zuani et al. 2024, Nat Commun）按切片分层：各自 评分/相关/箱线
B. 第 5 例 lepidic Visium 样本（lepidic.rar，无坐标，仅签名级推断）：
   防御轴 9 基因 / m6A 阅读器 / CTM2024 Table S1 lepidic-acinar 签名打分
C. 输出 per-section 统计表（数值断言与 pathB B1_section_level_stats.tsv 核对）
红线：spot 级仅描述；切片级统计 n=4；空间结果定性为组织层面探索性线索。
"""
import json, os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sc.settings.verbosity = 0
BASE = "/Users/liubao/WorkBuddy/2026-10-01-13-01-59/WP1_spatial"
PATHB_SP = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/spatial_data"
SAMPLES = ["P10_T1", "P17_T1", "P24_T1", "P11_T1"]  # 3 LUAD/NSCLC-NOS + 1 LUSC

DEFENSE = ["SLC7A11","SLC3A2","GCLC","GCLM","GOT2","TFRC","VDAC2","IGF2BP2","IGF2BP3"]
READERS = ["IGF2BP1","IGF2BP2","IGF2BP3","YTHDF1","YTHDF2","YTHDF3",
           "YTHDC1","YTHDC2","HNRNPA2B1","HNRNPC","FMR1","RBMX"]
EPITHELIAL = ["EPCAM","KRT18","KRT19"]
SIG = json.load(open(f"{BASE}/CTM_lepidic_acinar_signatures.json"))
LEP, ACI = SIG["lepidic"], SIG["acinar"]

# ============ A. 4 例按切片分层 ============
scores = pd.read_csv(f"{PATHB_SP}/B1_spot_scores.tsv", sep="\t")
scores = scores.rename(columns={scores.columns[0]: "spot"})
scores["region"] = np.where(scores["epi_score"] > 0, "tumor-rich", "stroma-rich")

# 数值断言：与 pathB 已产出的切片级统计逐位核对
old = pd.read_csv(f"{PATHB_SP}/B1_section_level_stats.tsv", sep="\t").set_index("section")
per_rows = []
for s in SAMPLES:
    g = scores[scores["sample"] == s]
    t = g[g.region == "tumor-rich"]; m = g[g.region == "stroma-rich"]
    rho = stats.spearmanr(g.defense_score, g.reader_score).statistic
    per_rows.append(dict(section=s, n_spots=len(g), n_tumor=len(t), n_stroma=len(m),
        defense_median_tumor=t.defense_score.median(),
        defense_median_stroma=m.defense_score.median(),
        defense_diff=t.defense_score.median() - m.defense_score.median(),
        reader_diff=t.reader_score.median() - m.reader_score.median(),
        spot_spearman_rho=rho))
per = pd.DataFrame(per_rows)
for _, r in per.iterrows():
    o = old.loc[r["section"]]
    assert abs(r["defense_diff"] - o["defense_diff"]) < 2e-4, (r["section"], "def_diff mismatch")
    assert abs(r["spot_spearman_rho"] - o["spot_spearman_rho"]) < 2e-4, (r["section"], "rho mismatch")
print("断言通过：4 例切片级统计与 pathB B1_section_level_stats.tsv 一致")
per.to_csv(f"{BASE}/WP1_per_section_stats.tsv", sep="\t", index=False, float_format="%.4f")

# 切片级检验（n=4，与 B1_regen 同口径）
w_def, p_def = stats.wilcoxon(per.defense_diff, alternative="greater")
w_read, p_read = stats.wilcoxon(per.reader_diff, alternative="greater")
z = np.arctanh(per.spot_spearman_rho.values)
t_z, p_z = stats.ttest_1samp(z, 0)
print(f"section-level: defense diff P={p_def:.4f} (4/4 positive), reader diff P={p_read:.4f}, rho Fisher-z P={p_z:.4f}")

# Fig WP1-A：每例 肿瘤区 vs 基质区 箱线（defense / reader 双面板 × 4 例）
fig, axes = plt.subplots(2, 4, figsize=(17, 7.5))
for j, s in enumerate(SAMPLES):
    sub = scores[scores["sample"] == s]
    for i, (col, ttl) in enumerate([("defense_score", "Defence-axis score"),
                                    ("reader_score", "m6A-reader score")]):
        data = [sub.loc[sub.region == r, col].values for r in ["tumor-rich", "stroma-rich"]]
        axes[i, j].boxplot(data, tick_labels=["tumor\nrich", "stroma\nrich"], showfliers=False)
        axes[i, j].set_title(f"{s}\n{ttl}", fontsize=9)
        if j == 0: axes[i, j].set_ylabel("spot score")
fig.suptitle("Per-section tumour-rich vs stroma-rich spot scores (E-MTAB-13530; "
             "spot display descriptive only — section-level tests in text, n = 4)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(f"{BASE}/WP1_FigA_每例肿瘤区vs基质区箱线.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# Fig WP1-B：每例 defense-reader 相关散点
fig, axes = plt.subplots(1, 4, figsize=(18, 4.6))
for ax, s in zip(axes, SAMPLES):
    g = scores[scores["sample"] == s]
    rho = per.loc[per.section == s, "spot_spearman_rho"].iloc[0]
    ax.scatter(g.defense_score, g.reader_score, s=4, alpha=0.3, edgecolors="none")
    ax.set_title(f"{s}\nSpearman ρ = {rho:.3f}", fontsize=10)
    ax.set_xlabel("Defence-axis score"); ax.set_ylabel("m6A-reader score")
fig.suptitle("Per-section defence–reader co-variation (spot display descriptive only)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(f"{BASE}/WP1_FigB_每例防御阅读器相关.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ============ B. lepidic 第 5 例签名打分 ============
a = sc.read_10x_mtx(f"{BASE}/lepidic/filtered_feature_bc_matrix", var_names="gene_symbols", cache=False)
a.var_names_make_unique()
print(f"lepidic 样本: {a.n_obs} spots × {a.n_vars} genes")
sc.pp.filter_cells(a, min_genes=200)
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)

def cov(df_genes):
    return [g for g in df_genes if g in a.var_names]
print("基因覆盖: defense", len(cov(DEFENSE)), "/9; readers", len(cov(READERS)), "/12;",
      "epi", len(cov(EPITHELIAL)), "/3; lepidic_sig", len(cov(LEP)), f"/{len(LEP)}; acinar_sig", len(cov(ACI)), f"/{len(ACI)}")
missing_epi = [g for g in EPITHELIAL if g not in a.var_names]
print("epi 缺失基因:", missing_epi)
assert len(cov(EPITHELIAL)) >= 2 and len(cov(LEP)) >= 3

sc.tl.score_genes(a, cov(DEFENSE), score_name="defense_score")
sc.tl.score_genes(a, cov(READERS), score_name="reader_score")
sc.tl.score_genes(a, cov(EPITHELIAL), score_name="epi_score")
sc.tl.score_genes(a, cov(LEP), score_name="lepidic_sig")
sc.tl.score_genes(a, cov(ACI), score_name="acinar_sig")

lep_df = a.obs[["defense_score","reader_score","epi_score","lepidic_sig","acinar_sig"]].copy()
lep_df["sample"] = "LEP"
lep_df["region"] = np.where(lep_df["epi_score"] > 0, "tumor-rich", "stroma-rich")
lep_df.to_csv(f"{BASE}/WP1_lepidic_spot_scores.tsv", sep="\t", float_format="%.5f")

# 与 4 例合并的描述性分布（评分各自样本内计算，跨样本比较仅描述性）
allk = scores[["sample","region","defense_score","reader_score"]].copy()
allk["lepidic_sig"] = np.nan; allk["acinar_sig"] = np.nan
lep_sub = lep_df[["sample","region","defense_score","reader_score","lepidic_sig","acinar_sig"]]
combo = pd.concat([allk, lep_sub], ignore_index=True)

# lepidic 样本内部：region 比较 + defense-reader 相关
t = lep_df[lep_df.region=="tumor-rich"]; m = lep_df[lep_df.region=="stroma-rich"]
res_lep = {}
for col in ["defense_score","reader_score","lepidic_sig","acinar_sig"]:
    res_lep[col] = dict(t_med=t[col].median(), s_med=m[col].median(),
                        diff=t[col].median()-m[col].median(),
                        p=stats.mannwhitneyu(t[col], m[col]).pvalue)
rho_lep = stats.spearmanr(lep_df.defense_score, lep_df.reader_score).statistic
p_rho_lep = stats.spearmanr(lep_df.defense_score, lep_df.reader_score).pvalue
print("\nlepidic 样本内部 (tumor-rich vs stroma-rich):")
for k, v in res_lep.items():
    print(f"  {k}: diff={v['diff']:+.4f} (MW P={v['p']:.3g})")
print(f"  defense-reader Spearman ρ = {rho_lep:.3f} (P={p_rho_lep:.3g})")

# Fig WP1-C：左-跨样本 defense score 分布；右-lepidic 内部四评分 region 比较
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
order = ["P24_T1","P17_T1","P10_T1","P11_T1","LEP"]
data = [combo.loc[combo["sample"]==s, "defense_score"].values for s in order]
axes[0].boxplot(data, tick_labels=[s.replace("_T1","")+("\n(lepidic)" if s=="LEP" else "") for s in order], showfliers=False)
axes[0].set_ylabel("Defence-axis spot score"); axes[0].set_title("Defence-axis score across sections\n(scored within each section; descriptive only)", fontsize=10)
cols = ["defense_score","reader_score","lepidic_sig","acinar_sig"]
lab = ["Defence\naxis","m6A\nreaders","Lepidic\nsignature","Acinar\nsignature"]
data2 = [[lep_df.loc[lep_df.region==r, c].values for c in cols] for r in ["tumor-rich","stroma-rich"]]
pos = [1,2,3,4, 6,7,8,9]
bp1 = axes[1].boxplot(data2[0], positions=[1,2,3,4], widths=0.7, showfliers=False)
bp2 = axes[1].boxplot(data2[1], positions=[6,7,8,9], widths=0.7, showfliers=False)
for b in bp1["boxes"]: b.set_color("#c53030")
for b in bp2["boxes"]: b.set_color("#2b6cb0")
axes[1].set_xticks(pos); axes[1].set_xticklabels(lab+lab, fontsize=8)
axes[1].axvline(5, color="grey", lw=0.8, ls="--")
ytxt = axes[1].get_ylim()[0] + 0.92*(axes[1].get_ylim()[1]-axes[1].get_ylim()[0])
axes[1].text(2.5, ytxt, "tumor-rich", ha="center", color="#c53030", fontsize=9)
axes[1].text(7.5, ytxt, "stroma-rich", ha="center", color="#2b6cb0", fontsize=9)
axes[1].set_ylabel("spot score")
axes[1].set_title(f"Lepidic-pattern section (5th Visium sample): signature scores by region\n"
                  f"lepidic signature: Δmedian={res_lep['lepidic_sig']['diff']:+.3f}; "
                  f"acinar: {res_lep['acinar_sig']['diff']:+.3f} (signature-level inference only)", fontsize=10)
fig.suptitle("Histologic-pattern context: lepidic section signature scoring", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(f"{BASE}/WP1_FigC_lepidic签名打分.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# 汇总表
summ = pd.DataFrame([
    dict(section="LEP(lepidic)", n_spots=len(lep_df), n_tumor=len(t), n_stroma=len(m),
         defense_median_tumor=t.defense_score.median(), defense_median_stroma=m.defense_score.median(),
         defense_diff=res_lep["defense_score"]["diff"],
         reader_diff=res_lep["reader_score"]["diff"],
         spot_spearman_rho=rho_lep,
         lepidic_sig_diff=res_lep["lepidic_sig"]["diff"], acinar_sig_diff=res_lep["acinar_sig"]["diff"])])
pd.concat([per, summ], ignore_index=True).to_csv(f"{BASE}/WP1_all_sections_stats.tsv", sep="\t", index=False, float_format="%.4f")
print("\nsaved WP1_FigA/B/C + WP1_per_section_stats.tsv + WP1_lepidic_spot_scores.tsv + WP1_all_sections_stats.tsv")
