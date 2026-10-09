#!/usr/bin/env python3
"""Compare the panel-based %RP (flawed) against the Gifford-faithful %RP
(RP reads / total Pseudomonas genome reads) across the three groups."""
import json, glob, os, math, statistics as st

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
healthy = [l.strip() for l in open("healthy_all.txt") if l.strip()]


def mwu(a, b):
    n1, n2 = len(a), len(b)
    U = sum(1 for x in a for y in b if x > y) + 0.5 * sum(1 for x in a for y in b if x == y)
    mu = n1 * n2 / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (U - mu) / sd if sd else 0
    return U, z, math.erfc(abs(z) / math.sqrt(2))


def load(s):
    f = f"out/{s}_rp.json"
    return json.load(open(f)) if os.path.exists(f) else None


groups = {
    "acute BALF": [(s, "acute") for s in acute],
    "chronic in vivo sputum": [(s, "invivo") for s in INVIVO],
    "chronic lab exponential": [(s, "exp") for s in EXP],
    "chronic lab stationary": [(s, "stat") for s in STAT],
    "healthy BALF": [(s, "healthy") for s in healthy],
}

print("=" * 92)
print("CORRECTED growth axis:  %RP = RP reads / total Pseudomonas genome reads   (Gifford 2014)")
print("=" * 92)
print(f"{'group':26s}{'n':>4s}{'median %RP':>12s}{'range':>18s}{'total reads':>14s}")
vals = {}
for g, items in groups.items():
    v = []
    tot = 0
    for s, _ in items:
        d = load(s)
        if d and d["total_pao1"] > 0:
            v.append(d["pctRP"]); tot += d["total_pao1"]
    vals[g] = v
    if v:
        print(f"{g:26s}{len(v):>4d}{st.median(v):>11.1f}%"
              f"{('['+format(min(v),'.1f')+', '+format(max(v),'.1f')+']'):>18s}{tot:>14,}")
    else:
        print(f"{g:26s}{0:>4d}{'--':>12s}{'(no signal)':>18s}{0:>14,}")

print("\n--- key contrasts (corrected) ---")
for a, b in [("acute BALF", "chronic in vivo sputum"),
             ("chronic in vivo sputum", "chronic lab exponential"),
             ("acute BALF", "chronic lab exponential"),
             ("chronic lab stationary", "chronic in vivo sputum")]:
    if vals[a] and vals[b]:
        U, z, p = mwu(vals[a], vals[b])
        print(f"  {a:24s} vs {b:24s}  medians {st.median(vals[a]):5.1f}% vs {st.median(vals[b]):5.1f}%"
              f"   p={p:.4f}")

# depth check on corrected metric
print("\n--- depth vs corrected %RP (acute) ---")
acute_rows = [(s, load(s)) for s in acute if load(s) and load(s)["total_pao1"] > 0]
def sp(x, y):
    def rank(a):
        s = sorted(range(len(a)), key=lambda i: a[i]); r = [0]*len(a)
        for p, i in enumerate(s): r[i] = p
        return r
    rx, ry = rank(x), rank(y); n = len(x)
    mx, my = sum(rx)/n, sum(ry)/n
    num = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))
    return num/den if den else float("nan")
if len(acute_rows) > 3:
    print(f"  Spearman(log total, %RP) = {sp([math.log(d['total_pao1']) for _,d in acute_rows],[d['pctRP'] for _,d in acute_rows]):+.3f}  n={len(acute_rows)}")
    print(f"  depth range: {min(d['total_pao1'] for _,d in acute_rows):,} - {max(d['total_pao1'] for _,d in acute_rows):,}")

json.dump({k: v for k, v in vals.items()}, open("rp_corrected_groups.json", "w"), indent=1)
print("\nwrote rp_corrected_groups.json")
