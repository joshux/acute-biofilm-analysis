#!/usr/bin/env python3
"""Corrected figure: quadrant placement on the Gifford-faithful growth axis."""
import json, os, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ABF = "/agent/workspace/abf"
os.chdir(ABF)
rows = json.load(open("quadrant_corrected.json"))
INVIVO = {"ERR2275087","ERR2275088","ERR2275089","ERR2275090","ERR2275091",
          "ERR2275092","ERR2275093","ERR2275094","ERR2275095","ERR2275096",
          "ERR2275097","ERR2275098","ERR2275099","ERR2275100","ERR2275101"}
EXP = {"ERR2275075","ERR2275076","ERR2275081","ERR2275082","ERR2591773",
       "ERR2591774","ERR2591777","ERR2591778","ERR2591781","ERR2591782",
       "ERR2591785","ERR2591786"}
STAT = {"ERR2275079","ERR2275080","ERR2275085","ERR2275086","ERR2591772",
        "ERR2591775","ERR2591776","ERR2591779","ERR2591780","ERR2591783","ERR2591784"}

def nsys(r):
    return sum(1 for s in ("alginate", "psl", "pel") if r[{"alginate":"alg","psl":"psl","pel":"pel"}[s]] >= 3)

groups = [("acute BALF (n=13)", "acute", "#c0392b", "o"),
          ("chronic CF sputum in vivo (n=15)", "chronic-invivo", "#1f6f8b", "s"),
          ("chronic lab exponential (n=12)", "chronic-exp", "#7f8c8d", "^"),
          ("chronic lab stationary (n=11)", "chronic-stat", "#6c3483", "v")]

fig, ax = plt.subplots(figsize=(8.6, 6.4))
rng = np.random.default_rng(11)
for lab, g, col, mk in groups:
    grp = [r for r in rows if r["group"] == g and r["tot_pao1"] >= 1000
           and r["pctRP_coding"] is not None]
    x = np.array([r["pctRP_coding"] for r in grp]) + rng.uniform(-0.25, 0.25, len(grp))
    y = np.array([nsys(r) for r in grp]) + rng.uniform(-0.08, 0.08, len(grp))
    ax.scatter(x, y, s=62, c=col, marker=mk, alpha=0.82, edgecolors="white",
               linewidths=0.8, label=lab, zorder=3)

ax.axvline(10, color="#444", ls="--", lw=1.2, zorder=1)
ax.axvline(5, color="#888", ls=":", lw=1.2, zorder=1)
ax.axhline(2, color="#444", ls="--", lw=1.2, zorder=1)
ax.text(10.25, 3.42, "FAST >10%RP", fontsize=9, color="#444")
ax.text(0.4, 3.42, "SLOW <5%RP", fontsize=9, color="#666")
ax.text(26.6, 2.12, "matrix-ON", fontsize=9, color="#444", ha="right")
ax.axvspan(10, 29, color="#f7d9d3", alpha=0.22, zorder=0)
ax.text(28.4, 0.28, "Kolpen acute\n(matrix-ON + FAST)", fontsize=9, ha="right",
        color="#8a4b3f", style="italic")
ax.text(0.6, 2.85, "Kolpen chronic\n(matrix-ON + SLOW)", fontsize=9, color="#1f4e5f",
        style="italic")

# group medians
for lab, g, col, mk in groups:
    grp = [r for r in rows if r["group"] == g and r["tot_pao1"] >= 1000
           and r["pctRP_coding"] is not None]
    if grp:
        ax.plot([st.median([r["pctRP_coding"] for r in grp])] * 2, [-0.32, -0.05],
                color=col, lw=2.5, zorder=2, solid_capstyle="butt")

ax.scatter([-0.7] * 9, [-0.5] * 9, s=42, c="#b7950b", marker="x",
           label="healthy BALF (n=9, no signal)", zorder=3)
ax.set_xlabel("growth axis   %RP  = RP reads / Pseudomonas coding reads", fontsize=11)
ax.set_ylabel("matrix axis   systems detected (alginate / psl / pel)", fontsize=11)
ax.set_title("Corrected two-dimension placement (Gifford-faithful %RP)\n"
             "lab anchors validate the metric: exponential 13.6% FAST, stationary 2.1% SLOW",
             fontsize=11, pad=12)
ax.set_xlim(-1.5, 29); ax.set_ylim(-0.7, 3.6)
ax.set_yticks([0, 1, 2, 3])
ax.legend(loc="lower left", fontsize=8.3, framealpha=0.93)
ax.grid(alpha=0.15, zorder=0)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
plt.tight_layout()
plt.savefig("fig_quadrant_corrected.png", dpi=170)
plt.savefig("fig_quadrant_corrected.svg")
print("wrote fig_quadrant_corrected.png / .svg")
