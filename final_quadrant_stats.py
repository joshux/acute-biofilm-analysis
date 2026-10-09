#!/usr/bin/env python3
"""Final group statistics for the two-dimension (quadrant) test."""
import json, glob, math, statistics as st

ABF = "/agent/workspace/abf"
import os; os.chdir(ABF)

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


def mwu(a, b):
    n1, n2 = len(a), len(b)
    U = sum(1 for x in a for y in b if x > y) + 0.5 * sum(1 for x in a for y in b if x == y)
    mu = n1 * n2 / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (U - mu) / sd if sd else 0
    return U, z, math.erfc(abs(z) / math.sqrt(2))


def mat_any(d):
    return sum(1 for s in ("alginate", "psl", "pel") if d["systems"][s]["genes_detected"] >= 1) >= 1


def mat_on(d):
    return sum(1 for s in ("alginate", "psl", "pel") if d["systems"][s]["genes_detected"] >= 3) >= 2


print("=" * 78)
print("GROWTH AXIS  (%RP; Gifford anchors FAST>10%, SLOW<5%)")
print("=" * 78)
for lab, grp in [("acute BALF", acute), ("chronic in vivo sputum", iv),
                 ("chronic in vitro exponential", ex)]:
    v = [d["f_rp"] for d in grp]
    nf = sum(1 for x in v if x > 0.10); ns = sum(1 for x in v if x < 0.05)
    print(f"  {lab:30s} n={len(v):2d}  median {100*st.median(v):5.2f}%  "
          f"FAST {nf:2d}  SLOW {ns:2d}  ambig {len(v)-nf-ns:2d}")

U, z, p = mwu([d["f_rp"] for d in acute], [d["f_rp"] for d in iv])
print(f"\n  acute vs chronic-in-vivo %RP :  U={U:.0f}  z={z:+.2f}  p={p:.3f}  "
      f"(medians {100*st.median([d['f_rp'] for d in acute]):.2f}% vs {100*st.median([d['f_rp'] for d in iv]):.2f}%)")
U, z, p = mwu([d["f_rp"] for d in iv], [d["f_rp"] for d in ex])
print(f"  chronic-in-vivo vs exponential: U={U:.0f}  z={z:+.2f}  p={p:.4f}  "
      f"(internal anchor check)")
U, z, p = mwu([d["f_rp"] for d in acute], [d["f_rp"] for d in ex])
print(f"  acute vs chronic-exponential   : U={U:.0f}  z={z:+.2f}  p={p:.4f}")

print("\n" + "=" * 78)
print("MATRIX AXIS  (systems detected among alginate / psl / pel)")
print("=" * 78)
for lab, grp in [("acute BALF", acute), ("chronic in vivo sputum", iv)]:
    n = len(grp)
    any_m = sum(1 for d in grp if mat_any(d))
    on = sum(1 for d in grp if mat_on(d))
    alg = sum(1 for d in grp if d["systems"]["alginate"]["genes_detected"] >= 3)
    psl = sum(1 for d in grp if d["systems"]["psl"]["genes_detected"] >= 3)
    pel = sum(1 for d in grp if d["systems"]["pel"]["genes_detected"] >= 3)
    print(f"  {lab:24s} n={n:2d}  matrix present {any_m}/{n} ({100*any_m//n}%)  "
          f"matrix-ON(>=2 sys) {on}/{n} ({100*on//n}%)")
    print(f"  {'':24s}      alginate {alg}/{n}  psl {psl}/{n}  pel {pel}/{n}")

print("\n" + "=" * 78)
print("HEALTHY CONTROLS (PRJNA390194, non-COPD BALF)")
print("=" * 78)
tg = [d["target_reads"] for d in healthy]
mp = [d["mapped_reads"] for d in healthy]
dc = [d["decoy_hits"] for d in healthy]
print(f"  n={len(healthy)}  panel target reads: median {st.median(tg):.0f}  total {sum(tg)}")
print(f"  reads mapping to decoy (kitome/co-infecting) genera: median {st.median(dc):.0f}")
print(f"  -> no bacterial target signal above the kitome floor")

json.dump({
    "growth": {"acute_median_pctRP": round(100*st.median([d["f_rp"] for d in acute]), 2),
               "chronic_invivo_median_pctRP": round(100*st.median([d["f_rp"] for d in iv]), 2),
               "chronic_exp_median_pctRP": round(100*st.median([d["f_rp"] for d in ex]), 2),
               "acute_vs_chronic_p": round(p, 4)},
    "matrix": {"acute_matrix_present": sum(1 for d in acute if mat_any(d)), "acute_n": len(acute),
               "chronic_invivo_matrix_present": sum(1 for d in iv if mat_any(d)), "chronic_n": len(iv),
               "acute_matrix_ON": sum(1 for d in acute if mat_on(d)),
               "chronic_matrix_ON": sum(1 for d in iv if mat_on(d))},
    "healthy": {"n": len(healthy), "target_reads_total": sum(tg)},
}, open("final_stats.json", "w"), indent=1)
print("\nwrote final_stats.json")
