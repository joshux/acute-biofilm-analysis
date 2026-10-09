#!/usr/bin/env python3
"""Chronic arm, split by the dataset's OWN experimental groups.

PRJEB24688 (Rossi 2018 Nat Commun) contains three expression clusters:
  in vivo  = CF sputum (INVIVO)
  in vitro exponential = EXP + 3H cultures
  in vitro stationary  = STAT + 24H cultures

The lab cultures are the internal FAST/SLOW anchors, exactly analogous to
Gifford's growth-rate anchors. This removes the cross-dataset normalisation
problem: the question becomes where in vivo sputum sits between the two lab
anchors measured in the same experiment.

Outputs chronic_groups.csv and prints the ordering.
"""
import json, os, glob, csv, statistics as st

ABF = "/agent/workspace/abf"
os.chdir(ABF)

# subgroup assignment from ENA sample titles (verified 2026-10-10)
INVIVO = {"ERR2275087","ERR2275088","ERR2275089","ERR2275090","ERR2275091",
          "ERR2275092","ERR2275093","ERR2275094","ERR2275095","ERR2275096",
          "ERR2275097","ERR2275098","ERR2275099","ERR2275100","ERR2275101"}
EXP = {"ERR2275075","ERR2275076","ERR2275081","ERR2275082",                       # EXP
       "ERR2591773","ERR2591774","ERR2591777","ERR2591778",
       "ERR2591781","ERR2591782","ERR2591785","ERR2591786"}                       # 3H
STAT = {"ERR2275079","ERR2275080","ERR2275085","ERR2275086",                      # STAT
        "ERR2591772","ERR2591775","ERR2591776","ERR2591779",
        "ERR2591780","ERR2591783","ERR2591784"}                                   # 24H

rows = []
for f in sorted(glob.glob("out/ERR*_pseu_v3.json")):
    d = json.load(open(f))
    s = d["srr"]
    grp = "in_vivo_sputum" if s in INVIVO else "in_vitro_exp" if s in EXP else \
          "in_vitro_stat" if s in STAT else "?"
    rows.append({"srr": s, "subgroup": grp, "target": d["target_reads"],
                 "f_rp": d["f_rp"], "matrix_on": d["matrix_on"],
                 "alg": d["systems"]["alginate"]["genes_detected"],
                 "psl": d["systems"]["psl"]["genes_detected"],
                 "pel": d["systems"]["pel"]["genes_detected"]})

with open("chronic_groups.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print(f"{'sample':12s}{'subgroup':16s}{'target':>9s}{'%RP':>8s}{'matrix':>7s}{'alg':>5s}{'psl':>5s}{'pel':>5s}")
print("-" * 67)
for r in sorted(rows, key=lambda r: (r["subgroup"], -r["f_rp"])):
    print(f"{r['srr']:12s}{r['subgroup']:16s}{r['target']:>9,}{r['f_rp']*100:>7.2f}%"
          f"{'ON' if r['matrix_on'] else 'off':>7s}{r['alg']:>5}{r['psl']:>5}{r['pel']:>5}")

print("\n=== %RP by subgroup (the growth axis, same panel, same experiment) ===")
for g in ("in_vitro_exp", "in_vivo_sputum", "in_vitro_stat"):
    v = [r["f_rp"] for r in rows if r["subgroup"] == g]
    if not v: continue
    print(f"  {g:16s} n={len(v):2d}  median {100*st.median(v):5.2f}%  "
          f"mean {100*st.mean(v):5.2f}%  range [{100*min(v):.2f}, {100*max(v):.2f}]")

# rank test: is in_vivo below in_vitro_exp?
iv = [r["f_rp"] for r in rows if r["subgroup"] == "in_vivo_sputum"]
ex = [r["f_rp"] for r in rows if r["subgroup"] == "in_vitro_exp"]
sx = [r["f_rp"] for r in rows if r["subgroup"] == "in_vitro_stat"]
print(f"\n  in_vivo < in_vitro_exp ?  median {100*st.median(iv):.2f}% vs {100*st.median(ex):.2f}%")
print(f"  in_vivo > in_vitro_stat?  median {100*st.median(iv):.2f}% vs {100*st.median(sx):.2f}%")

# Mann-Whitney U (exact-ish via normal approx) for in_vivo vs exp
def mwu(a, b):
    import itertools
    n1, n2 = len(a), len(b)
    U = sum(1 for x in a for y in b if x > y) + 0.5*sum(1 for x in a for y in b if x == y)
    mu = n1*n2/2
    import math
    sd = math.sqrt(n1*n2*(n1+n2+1)/12)
    z = (U-mu)/sd if sd else 0
    # two-sided p from normal
    p = math.erfc(abs(z)/math.sqrt(2))
    return U, z, p

U, z, p = mwu(iv, ex)
print(f"  Mann-Whitney in_vivo vs in_vitro_exp: U={U:.0f} z={z:+.2f} p={p:.4f}")
U2, z2, p2 = mwu(iv, sx)
print(f"  Mann-Whitney in_vivo vs in_vitro_stat: U={U2:.0f} z={z2:+.2f} p={p2:.4f}")
json.dump({"in_vivo": iv, "in_vitro_exp": ex, "in_vitro_stat": sx},
          open("chronic_growth.json", "w"), indent=1)
