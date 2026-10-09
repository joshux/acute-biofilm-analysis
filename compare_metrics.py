#!/usr/bin/env python3
"""Compare three growth metrics across groups.

  %RP_panel  = RP / (matrix+rp+og) on the multi-genome panel   <- the original (flawed)
  %RP_genome = RP / all PAO1-mapping reads                     <- Gifford-faithful
  %RP_coding = RP / PAO1-mapping reads outside rRNA/tRNA       <- rRNA-robust (cross-dataset clean)
"""
import json, os, math, statistics as st

ABF = "/agent/workspace/abf"
os.chdir(ABF)

INVIVO = {"ERR2275087","ERR2275088","ERR2275089","ERR2275090","ERR2275091",
          "ERR2275092","ERR2275093","ERR2275094","ERR2275095","ERR2275096",
          "ERR2275097","ERR2275098","ERR2275099","ERR2275100","ERR2275101"}
EXP = {"ERR2275075","ERR2275076","ERR2275081","ERR2275082","ERR2591773",
       "ERR2591774","ERR2591777","ERR2591778","ERR2591781","ERR2591782",
       "ERR2591785","ERR2591786"}
STAT = {"ERR2275079","ERR2275080","ERR2275085","ERR2275086","ERR2591772",
        "ERR2591775","ERR2591776","ERR2591779","ERR2591780","ERR2591783","ERR2591784"}
acute = [l.strip() for l in open("acute_runs.txt") if l.strip()]


def mwu(a, b):
    n1, n2 = len(a), len(b)
    U = sum(1 for x in a for y in b if x > y) + 0.5 * sum(1 for x in a for y in b if x == y)
    mu = n1 * n2 / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (U - mu) / sd if sd else 0
    return math.erfc(abs(z) / math.sqrt(2))


def vals(samples, key):
    out = []
    for s in samples:
        f = f"out/{s}_{key}.json"
        if not os.path.exists(f):
            continue
        d = json.load(open(f))
        v = d.get({"rp": "pctRP", "rp2": "pctRP_coding"}[key])
        if key == "rp" and d["total_pao1"] == 0:
            continue
        if key == "rp2" and not d.get("coding_reads"):
            continue
        if v is not None:
            out.append(v)
    return out


groups = [("acute BALF", acute), ("chronic in vivo", sorted(INVIVO)),
          ("chronic lab exp", sorted(EXP)), ("chronic lab stat", sorted(STAT))]

print("=" * 88)
print("GROWTH METRIC COMPARISON  (medians; n in parens)")
print("=" * 88)
print(f"{'group':20s}{'%RP_panel':>16s}{'%RP_genome':>16s}{'%RP_coding':>16s}")
tab = {}
for g, ss in groups:
    a = vals(ss, "rp"); b = vals(ss, "rp2")
    p = [json.load(open(f"out/{s}_pseu_v3.json"))["f_rp"] * 100
         for s in ss if os.path.exists(f"out/{s}_pseu_v3.json")]
    tab[g] = {"panel": p, "genome": a, "coding": b}
    print(f"{g:20s}{st.median(p):>9.1f}%({len(p):2d}){st.median(a):>9.1f}%({len(a):2d})"
          f"{st.median(b):>9.1f}%({len(b):2d})")

print("\n--- contrasts on the rRNA-robust metric (%RP_coding) ---")
pairs = [("acute BALF", "chronic in vivo"), ("acute BALF", "chronic lab exp"),
         ("chronic lab exp", "chronic lab stat"), ("chronic in vivo", "chronic lab stat"),
         ("chronic in vivo", "chronic lab exp")]
for a, b in pairs:
    va, vb = tab[a]["coding"], tab[b]["coding"]
    if va and vb:
        print(f"  {a:20s} vs {b:20s}  {st.median(va):5.1f}% vs {st.median(vb):5.1f}%   "
              f"Mann-Whitney p={mwu(va, vb):.4f}")

print("\n--- Gifford anchor check (metric validity) ---")
print(f"  lab exponential (fast growth expected): {st.median(tab['chronic lab exp']['coding']):.1f}%")
print(f"  lab stationary  (slow growth expected): {st.median(tab['chronic lab stat']['coding']):.1f}%")

json.dump({g: {k: v for k, v in d.items()} for g, d in tab.items()},
          open("growth_metric_compare.json", "w"), indent=1)
print("\nwrote growth_metric_compare.json")
