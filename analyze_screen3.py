#!/usr/bin/env python3
"""analyze_screen3.py — compare the corrected 40k screen against the committed claims."""
import json, collections

S = json.load(open("screen3_ALL.json"))

AXIS = ["pseudomonas", "klebsiella", "escherichia", "staphylococcus",
        "acinetobacter", "haemophilus", "streptococcus", "stenotrophomonas"]

# --- per-patient deepest run ---
byp = collections.defaultdict(list)
for srr, r in S.items():
    if "error" in r: continue
    byp[r["patient"]].append(r)

print("=== deepest run per patient: dominant matrix-capable organism (>=20 reads) ===")
deep = {}
for p, rs in byp.items():
    rs = sorted(rs, key=lambda x: -x["read_count"])
    top = rs[0]
    cnt = top["counts"]
    dom = max(cnt, key=cnt.get) if cnt else None
    deep[p] = (top, dom, cnt)
    st = top["status"]
    print(f"{p:7s} {st:6s} {top['srr']:13s} total={top['total']:>6} dom={str(dom):16s} "
          f"top3={sorted(cnt.items(), key=lambda x:-x[1])[:3]}")

# --- the HS243 claim ---
print("\n=== HS243 CLAIM CHECK ===")
for srr, r in sorted(S.items(), key=lambda x: -x[1].get("read_count", 0)):
    if r.get("patient") == "HS243":
        print(f"  {srr} reads={r['read_count']:>10,} total_panel={r['total']:>6} counts={r['counts']}")

# --- arm-level: how many patients have a matrix-capable dominant at depth ---
print("\n=== ARM-LEVEL dominant-organism composition (deepest run/patient) ===")
for arm in ("2BP", "No-BP"):
    c = collections.Counter()
    for p, (top, dom, cnt) in deep.items():
        if top["status"] == arm:
            c[dom] += 1
    print(f"  {arm:6s} n={sum(c.values())}: {dict(c.most_common())}")

# --- depth distribution by arm ---
print("\n=== total panel reads at 40k, by arm (deepest run/patient) ===")
for arm in ("2BP", "No-BP"):
    v = sorted([top["total"] for p, (top, dom, cnt) in deep.items() if top["status"] == arm])
    if v:
        import statistics
        print(f"  {arm:6s} n={len(v)} median={statistics.median(v):.0f} "
              f"min={v[0]} max={v[-1]}  >=100: {sum(1 for x in v if x>=100)}")

json.dump({p: {"srr": t["srr"], "status": t["status"], "dom": d, "counts": c,
               "total": t["total"]} for p, (t, d, c) in deep.items()},
          open("screen3_deepest.json", "w"), indent=1)
