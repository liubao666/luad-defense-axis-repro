# -*- coding: utf-8 -*-
"""模块二出图：scTenifoldKnk 虚拟敲除 IGF2BP3 结果 -> Fig_B5 排名图 + 防御轴成员排名表"""
import json
import urllib.request
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/knk"
OUT = "/Users/liubao/WorkBuddy/2026-09-20-10-55-46/生信论文审稿/pathB/spatial_data"

AXIS = {"SLC7A11": "ENSG00000151012", "SLC3A2": "ENSG00000168003", "GCLC": "ENSG00000001084",
        "GCLM": "ENSG00000023909", "GOT2": "ENSG00000125166", "TFRC": "ENSG00000072274",
        "VDAC2": "ENSG00000165637", "IGF2BP2": "ENSG00000073792", "IGF2BP3": "ENSG00000136231"}

df = pd.read_csv(f"{BASE}/B2_knockout_IGF2BP3_full.csv")
print("列名:", list(df.columns), "行数:", len(df))
# scTenifoldKnk 输出: gene / FC / p.value / p.adj / Z（不同版本列名略有差异，自动识别）
gcol = [c for c in df.columns if c.lower() in ("gene", "genes")][0]
zcol = [c for c in df.columns if c.lower() in ("z", "zscore", "z.score", "dz")][0]
pcol = [c for c in df.columns if "p.adj" in c.lower() or "padj" in c.lower() or "fdr" in c.lower()]
pcol = pcol[0] if pcol else None
df = df.rename(columns={gcol: "ensg", zcol: "z"})
df = df.sort_values("z", ascending=False).reset_index(drop=True)
df["rank"] = np.arange(1, len(df) + 1)

# ENSG -> symbol：防御轴已知，Top 30 经 mygene 补映射
sym_map = {e: s for s, e in AXIS.items()}
need = [g for g in df.head(30)["ensg"] if g not in sym_map]
for i in range(0, len(need), 200):
    batch = need[i:i + 200]
    q = " OR ".join(f"ensembl.gene:{g}" for g in batch)
    req = urllib.request.Request(
        "https://mygene.info/v3/query",
        data=f"q={urllib.parse.quote(q)}&fields=symbol,ensembl.gene&species=human&size=500".encode(),
        headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            hits = json.load(r).get("hits", [])
        for h in hits:
            ens = h.get("ensembl", {})
            ens = ens.get("gene", ens) if isinstance(ens, dict) else None
            if isinstance(ens, list):
                ens = ens[0] if ens else None
            if ens and h.get("symbol"):
                sym_map.setdefault(ens, h["symbol"])
    except Exception as ex:
        print("mygene 批次失败:", ex)
df["symbol"] = df["ensg"].map(sym_map).fillna(df["ensg"])

# ---- 防御轴成员排名表 ----
axis_df = df[df["ensg"].isin(AXIS.values())][["symbol", "ensg", "rank", "z"] + ([pcol] if pcol else [])]
axis_df = axis_df.sort_values("rank")
axis_df.to_csv(f"{OUT}/B2_防御轴成员_扰动排名表.tsv", sep="\t", index=False)
print("\n防御轴成员扰动排名:\n", axis_df.to_string(index=False))

# ---- Fig_B5 排名图 ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6.5), gridspec_kw={"width_ratios": [2.2, 1]})
x = df["rank"].values
y = df["z"].values
colors = np.where(df["ensg"].isin(AXIS.values()), "#d62728", "#bdbdbd")
ax1.scatter(x, y, c=colors, s=8, alpha=0.7)
ax1.axhline(0, color="#666666", lw=0.8, ls="--")
ax1.set_xlabel("Perturbation rank (all genes)")
ax1.set_ylabel("Perturbation z-score (dZ)")
ax1.set_title("Virtual knockout of IGF2BP3 (scTenifoldKnk)\nred = 9-gene defence-axis members")
# 标注 Top 8 + 全部轴成员
lab = pd.concat([df.head(8), df[df["ensg"].isin(AXIS.values())]]).drop_duplicates("ensg")
for _, r in lab.iterrows():
    ax1.annotate(r["symbol"], (r["rank"], r["z"]), textcoords="offset points",
                 xytext=(4, 4), fontsize=7.5, color="#d62728" if r["ensg"] in AXIS.values() else "#333333")
red_proxy = plt.Line2D([], [], marker="o", ls="", color="#d62728", label="Defence-axis member (n=9)")
grey_proxy = plt.Line2D([], [], marker="o", ls="", color="#bdbdbd", label="Other genes")
ax1.legend(handles=[red_proxy, grey_proxy], fontsize=8, loc="lower left")

# 右图：轴成员 z 值条形图
ad = axis_df.sort_values("z")
ax2.barh(ad["symbol"], ad["z"], color="#d62728", edgecolor="#8c1a13")
ax2.axvline(0, color="#666666", lw=0.8)
ax2.set_xlabel("Perturbation z-score")
ax2.set_title("Defence-axis members\nrank by |perturbation|")
for i, (z, rk) in enumerate(zip(ad["z"], ad["rank"])):
    ax2.text(z + (0.3 if z >= 0 else -0.3), i, f"#{rk}", va="center",
             ha="left" if z >= 0 else "right", fontsize=7.5, color="#555555")
plt.tight_layout()
plt.savefig(f"{OUT}/Fig_B5_虚拟敲除IGF2BP3_扰动排名.png", dpi=300, bbox_inches="tight")
print("\nFig_B5 已输出")
