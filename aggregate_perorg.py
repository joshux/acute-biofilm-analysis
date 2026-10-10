#!/usr/bin/env python3
"""Aggregate the per-organism two-axis results across the acute cohort."""
import json, os, glob, collections, statistics as st

ABF = "/agent/workspace/abf"
os.chdir(ABF)

rows = []
for f in glob.glob("final_out/*.json"):
    d = json.load(open(f))
    rows.append(d)

usable = [d for d in rows if d.get("pctRP") is not None and d.get("coding", 0) >= 1000]
print(f"samples scored: {len(rows)}   usable (coding>=1000): {len(usable)}")

print("\n=== dominant organism across scored samples ===")
for g, n in collections.Counter(d["dominant"] for d in rows).most_common():
    print(f"  {g:18s} {n}")

print("\n=== per-organism results (coding>=1000) ===")
byorg = collections.defaultdict(list)
for d in usable:
    byorg[d["org"]].append(d)
print(f"{'org':7s}{'n':>4s}{'median %RP':>12s}{'%RP range':>18s}{'matrix-ON':>11s}")
for org, v in sorted(byorg.items(), key=lambda x: -len(x[1])):
    p = [d["pctRP"] for d in v]
    mon = sum(1 for d in v if d["matrix_on"])
    print(f"{org:7s}{len(v):>4d}{st.median(p):>11.1f}%"
          f"{('['+format(min(p),'.1f')+', '+format(max(p),'.1f')+']'):>18s}{mon:>6d}/{len(v)}")

print("\n=== per-sample table ===")
print(f"{'sample':13s}{'organism':14s}{'coding':>8s}{'%RP':>8s}{'matrix':>8s}  systems")
for d in sorted(usable, key=lambda d: (d["org"], -d["pctRP"])):
    print(f"{d['srr']:13s}{d['dominant']:14s}{d['coding']:>8,}{d['pctRP']:>7.1f}%"
          f"{('ON' if d['matrix_on'] else 'off'):>8s}  {';'.join(d['systems_on'])}")

# quadrant table across all organisms
print("\n=== QUADRANT (all organisms pooled) ===")
FAST, SLOW = 10.0, 5.0
print(f"{'':12s}{'FAST':>7s}{'ambig':>7s}{'SLOW':>7s}")
for mc in ("ON", "off"):
    cells = {"FAST": 0, "ambig": 0, "SLOW": 0}
    for d in usable:
        on = "ON" if d["matrix_on"] else "off"
        if on != mc:
            continue
        g = "FAST" if d["pctRP"] > FAST else ("SLOW" if d["pctRP"] < SLOW else "ambig")
        cells[g] += 1
    print(f"{mc:12s}{cells['FAST']:>7d}{cells['ambig']:>7d}{cells['SLOW']:>7d}")

json.dump({"n_scored": len(rows), "n_usable": len(usable),
           "by_org": {k: len(v) for k, v in byorg.items()}},
          open("perorg_summary.json", "w"), indent=1)
