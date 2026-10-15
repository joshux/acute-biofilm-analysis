#!/usr/bin/env python3
"""arm_stats.py — arm-level statistics on the corrected screen + partial item12."""
import json, collections, itertools, os

S = json.load(open("screen3_ALL.json"))
MC = {"pseudomonas", "klebsiella", "escherichia", "staphylococcus", "acinetobacter"}

# deepest run per patient, matrix-capable read count
byp = collections.defaultdict(list)
for srr, r in S.items():
    if "error" in r: continue
    byp[r["patient"]].append(r)

def fisher(a, b, c, d):
    """two-sided Fisher exact p (small tables), a/b = arm1 pos/neg, c/d = arm2."""
    from math import comb
    n = a + b + c + d
    def p(a, b, c, d):
        return comb(a+b, a) * comb(c+d, c) / comb(n, a+c)
    p0 = p(a, b, c, d)
    tot = 0.0
    for aa in range(0, a+b+1):
        for cc in range(0, c+d+1):
            bb = (a+b) - aa; dd = (c+d) - cc
            if aa+bb != a+b or cc+dd != c+d: continue
            pp = p(aa, bb, cc, dd)
            if pp <= p0 + 1e-12: tot += pp
    return min(1.0, tot)

rows = {}
for p, rs in byp.items():
    # best run per patient = most matrix-capable reads
    best = max(rs, key=lambda r: sum(v for k, v in r["counts"].items() if k in MC))
    n = sum(v for k, v in best["counts"].items() if k in MC)
    dom = max(best["counts"], key=best["counts"].get) if best["counts"] else None
    rows[p] = {"status": best["status"], "srr": best["srr"], "matcap": n, "dom": dom,
               "total": best["total"]}

for thr in (50, 100, 300, 1000):
    a = sum(1 for r in rows.values() if r["status"] == "2BP" and r["matcap"] >= thr)
    b = sum(1 for r in rows.values() if r["status"] == "2BP" and r["matcap"] < thr)
    c = sum(1 for r in rows.values() if r["status"] == "No-BP" and r["matcap"] >= thr)
    d = sum(1 for r in rows.values() if r["status"] == "No-BP" and r["matcap"] < thr)
    print(f"matrix-capable >= {thr:>5}: 2BP {a}/{a+b}  No-BP {c}/{c+d}  Fisher p={fisher(a,b,c,d):.4f}")

# depth balance
import statistics
for arm in ("2BP", "No-BP"):
    v = [r["total"] for r in rows.values() if r["status"] == arm]
    print(f"{arm}: n={len(v)} median total panel reads={statistics.median(v):.0f} "
          f"mean={statistics.mean(v):.0f}")

# dominant organism by arm
print("\ndominant organism (best matrix-capable run per patient):")
for arm in ("2BP", "No-BP"):
    c = collections.Counter(r["dom"] for r in rows.values() if r["status"] == arm)
    print(f"  {arm}: {dict(c.most_common())}")

# --- partial item12 ---
print("\n=== item12 (corrected-panel two-axis) results so far ===")
D = "item12"
if os.path.isdir(D):
    res = []
    for f in sorted(os.listdir(D)):
        if not f.endswith(".json"): continue
        d = json.load(open(os.path.join(D, f)))
        if d.get("pctRP") is None: continue
        res.append(d)
    for d in sorted(res, key=lambda x: -x["target_reads"]):
        print(f"  {d['srr']} {d['dominant']:15s} reads={d['target_reads']:>7} "
              f"%RP={d['pctRP']:>6} on={d['systems_on']} eps={d.get('eps_on')}")
    print(f"  n quantified = {len(res)}")
    print(f"  matrix-ON (any 2 systems): {sum(1 for d in res if d['matrix_on'])}/{len(res)}")
    print(f"  matrix-ON (2 EPS systems): {sum(1 for d in res if d.get('matrix_on_eps'))}/{len(res)}")
