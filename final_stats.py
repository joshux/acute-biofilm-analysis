#!/usr/bin/env python3
"""Final statistics for the cohort, with correct sign test and Wilcoxon."""
import json, math, statistics as st
from math import comb

a = json.load(open("/tmp/abf/analysis_final.json"))
coh = a["cohort"]


def sign_test(v):
    v = [x for x in v if x != 0]
    n = len(v); k = sum(1 for x in v if x > 0)
    up = sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n
    dn = sum(comb(n, i) for i in range(0, k + 1)) / 2 ** n
    return min(1.0, 2 * min(up, dn))


def wilcoxon(v):
    v = [x for x in v if x != 0]
    n = len(v)
    idx = sorted(range(n), key=lambda i: abs(v[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(v[idx[j + 1]]) == abs(v[idx[i]]):
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[idx[k]] = avg
        i = j + 1
    Wp = sum(ranks[i] for i in range(n) if v[i] > 0)
    Wm = sum(ranks[i] for i in range(n) if v[i] < 0)
    W = min(Wp, Wm)
    # normal approximation
    mu = n * (n + 1) / 4
    sd = math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z = (W - mu) / sd if sd else 0
    return W, z


print(f"COHORT n={len(coh)}  (Pseudomonas RNA fraction >= 0.5, target reads >= 100)\n")
for key, lab in [("B_full", "B = log R_matrix - log R_growth"), ("B_rp", "B = log R_matrix - log R_rp")]:
    v = [r[key] for r in coh]
    W, z = wilcoxon(v)
    print(f"{lab}")
    print(f"   median {st.median(v):+.3f}   mean {st.mean(v):+.3f}   sd {st.stdev(v):.3f}")
    print(f"   range [{min(v):+.2f}, {max(v):+.2f}]   n>0 = {sum(1 for x in v if x>0)}/{len(v)}")
    print(f"   sign test p = {sign_test(v):.3f}   Wilcoxon W={W:.1f} z={z:+.2f}")
    print()

print("Per-sample class shares (fraction of target reads):")
print(f"{'sample':13s}{'mat':>7s}{'rp':>7s}{'og':>7s}")
for r in coh:
    print(f"{r['srr']:13s}{r['f_mat']:>7.3f}{r['f_rp']:>7.3f}{r['f_og']:>7.3f}")

mat = [r["f_mat"] for r in coh]
print(f"\nmatrix share: median {st.median(mat):.3f}  range [{min(mat):.3f},{max(mat):.3f}]")
print(f"growth share (rp+og): median {st.median([1-m for m in mat]):.3f}")

# matrix gene census
import collections
genes = collections.Counter(); reads = collections.Counter()
for r in coh:
    d = json.load(open(f"/tmp/abf/{r['srr']}_pseu_q3.json"))
    for g, cls, acc, c, L in d["top_genes"]:
        if cls == "matrix":
            genes[g] += 1; reads[g] += c
print(f"\nmatrix genes detected: {len(genes)} distinct / 101 in panel; {sum(reads.values())} reads total")
for g, c in sorted(genes.items(), key=lambda x: -x[1])[:20]:
    print(f"   {g:8s} {c:2d}/10 samples  {reads[g]:4d} reads")
json.dump({"genes": dict(genes), "reads": dict(reads)}, open("/tmp/abf/matrix_gene_census.json", "w"), indent=1)
