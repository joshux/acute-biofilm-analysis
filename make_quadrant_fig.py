#!/usr/bin/env python3
"""Figure: quadrant placement across groups (growth axis x matrix axis)."""
import json, glob, os, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ABF = "/agent/workspace/abf"
os.chdir(ABF)
INVIVO = {"ERR2275087","ERR2275088","ERR2275089","ERR2275090","ERR2275091",
          "ERR2275092","ERR2275093","ERR2275094","ERR2275095","ERR2275096",
          "ERR2275097","ERR2275098","ERR2275099","ERR2275100","ERR2275101"}
EXP = {"ERR2275075","ERR2275076","ERR2275081","ERR2275082","ERR2591773",
       "ERR2591774","ERR2591777","ERR2591778","ERR2591781","ERR2591782",
       "ERR2591785","ERR2591786"}
acute = [json.load(open(f)) for f in glob.glob("out/SRR2734*_pseu_v3.json")]
chron = [json.load(open(f)) for f in glob.glob("out/ERR*_pseu_v3.json")]
healthy = [json.load(open(f)) for f in glob.glob("out/SRR5677*_pseu_v3.json")]
iv = [d for d in chron if d["srr"] in INVIVO]
ex = [d for d in chron if d["srr"] in EXP]


def nsys(d):
    return sum(1 for s in ("alginate", "psl", "pel") if d["systems"][s]["genes_detected"] >= 3)


def jit(v, amt=0.06):
    return np.array(v) + np.random.default_rng(7).uniform(-amt, amt, len(v))


fig, ax = plt.subplots(figsize=(8.2, 6.2))

groups = [("acute BALF (n=15)", acute, "#c0392b", "o"),
          ("chronic CF sputum (n=15)", iv, "#1f6f8b", "s"),
          ("chronic lab exponential (n=12)", ex, "#7f8c8d", "^")]
for lab, grp, col, mk in groups:
    x = jit([d["f_rp"] * 100 for d in grp])
    y = jit([nsys(d) for d in grp], 0.09)
    ax.scatter(x, y, s=64, c=col, marker=mk, alpha=0.8, edgecolors="white",
               linewidths=0.8, label=lab, zorder=3)

# Gifford anchors
ax.axvline(10, color="#444", ls="--", lw=1.2, zorder=1)
ax.axvline(5, color="#888", ls=":", lw=1.2, zorder=1)
ax.axhline(2, color="#444", ls="--", lw=1.2, zorder=1)
ax.text(10.3, 3.35, "FAST >10%RP", fontsize=9, color="#444")
ax.text(0.3, 3.35, "SLOW <5%RP", fontsize=9, color="#666")
ax.text(26.5, 2.12, "matrix-ON", fontsize=9, color="#444", ha="right")

# quadrant shading
ax.axvspan(10, 30, ymin=0, ymax=1, color="#f7d9d3", alpha=0.25, zorder=0)
ax.text(27.5, 0.25, "Kolpen\n(matrix-ON + FAST)", fontsize=9, ha="right",
        color="#8a4b3f", style="italic")
ax.text(1.0, 2.9, "classical\n(matrix-ON + SLOW)", fontsize=9, color="#1f4e5f",
        style="italic")

# healthy: zero signal -> mark at axis edge
ax.scatter([-0.6]*len(healthy), [-0.28]*len(healthy), s=48, c="#b7950b", marker="x",
           label="healthy BALF (n=9, no signal)", zorder=3)

ax.set_xlabel("growth axis   %RP  (ribosomal-protein transcript fraction)", fontsize=11)
ax.set_ylabel("matrix axis   systems detected (alginate / psl / pel)", fontsize=11)
ax.set_title("Two-dimension placement: acute vs chronic vs healthy\n"
             "Kolpen 2022 predicts acute = matrix-ON + FAST;  chronic = matrix-ON + SLOW",
             fontsize=11.5, pad=12)
ax.set_xlim(-1.5, 29)
ax.set_ylim(-0.6, 3.6)
ax.set_yticks([0, 1, 2, 3])
ax.legend(loc="lower left", fontsize=8.5, framealpha=0.92)
ax.grid(alpha=0.15, zorder=0)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
plt.tight_layout()
plt.savefig("fig_quadrant.png", dpi=170)
plt.savefig("fig_quadrant.svg")
print("wrote fig_quadrant.png")

# summary print for the write-up
print("\nmedian pctRP  acute %.2f%%   chronic-in-vivo %.2f%%   chronic-exp %.2f%%"
      % (100*st.median([d["f_rp"] for d in acute]),
         100*st.median([d["f_rp"] for d in iv]),
         100*st.median([d["f_rp"] for d in ex])))
print("matrix-ON  acute %d/15   chronic-in-vivo %d/15"
      % (sum(1 for d in acute if nsys(d) >= 2), sum(1 for d in iv if nsys(d) >= 2)))
