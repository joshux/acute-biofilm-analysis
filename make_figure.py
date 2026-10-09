#!/usr/bin/env python3
"""Figure: matrix vs growth investment in acute bacterial BALF."""
import json, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

a = json.load(open("/tmp/abf/analysis_final.json"))
coh = a["cohort"]
allr = a["all"]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# --- panel 1: balance distribution -------------------------------------
ax = axes[0]
bf = [r["B_full"] for r in coh]
br = [r["B_rp"] for r in coh]
y = np.arange(len(coh))
order = np.argsort(bf)
ax.barh(y - 0.2, [bf[i] for i in order], height=0.4, color="#2b6cb0", label="vs full growth panel")
ax.barh(y + 0.2, [br[i] for i in order], height=0.4, color="#dd6b20", label="vs ribosomal protein only")
ax.set_yticks(y)
ax.set_yticklabels([coh[i]["srr"][-5:] for i in order], fontsize=8)
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("B = log(matrix) − log(growth)   [log reads/kb]")
ax.set_title("Matrix vs growth balance\nper Pseudomonas-dominant sample", fontsize=11)
ax.legend(fontsize=8, loc="lower right")

# --- panel 2: depth confound -------------------------------------------
ax = axes[1]
lt = [math.log10(r["target"]) for r in coh]
lg = [math.log10(r["R_grw"]) for r in coh]
lm = [math.log10(r["R_mat"]) for r in coh]
ax.scatter(lt, lg, s=42, color="#2b6cb0", label="growth machinery")
ax.scatter(lt, lm, s=42, color="#38a169", marker="s", label="matrix machinery")
for x, yv, c in [(lt, lg, "#2b6cb0"), (lt, lm, "#38a169")]:
    z = np.polyfit(x, yv, 1)
    xs = np.linspace(min(x), max(x), 20)
    ax.plot(xs, np.polyval(z, xs), color=c, lw=1.4, alpha=0.7)
ax.set_xlabel("log10 target reads in panel (depth)")
ax.set_ylabel("log10 reads per kb")
ax.set_title("Both arms track sequencing depth\n(r ≈ 0.91 growth, 0.77 matrix)", fontsize=11)
ax.legend(fontsize=8)

# --- panel 3: partial correlation --------------------------------------
ax = axes[2]
labels = ["raw\nPearson", "partialling\nout depth"]
vals = [0.501, -0.759]
cols = ["#38a169", "#c53030"]
ax.bar(labels, vals, color=cols, width=0.5)
ax.axhline(0, color="k", lw=1)
for i, v in enumerate(vals):
    ax.text(i, v + (0.06 if v > 0 else -0.12), f"{v:+.2f}", ha="center", fontweight="bold")
ax.set_ylabel("correlation of log matrix vs log growth")
ax.set_ylim(-1.0, 0.75)
ax.set_title("The apparent co-expression is a\ndepth artifact: it inverts when depth is removed", fontsize=11)

plt.tight_layout()
plt.savefig("/tmp/abf/fig_balance.png", dpi=160, facecolor="white")
plt.savefig("/tmp/abf/fig_balance.svg", facecolor="white")
print("wrote fig_balance.png + fig_balance.svg")
