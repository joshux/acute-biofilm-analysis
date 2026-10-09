#!/usr/bin/env python3
"""Quadrant placement across the three groups.

Growth axis : f_rp (== %RP / 100).  Gifford anchors: FAST > 0.10, SLOW < 0.05.
Matrix axis : matrix-ON = >=2 of {alginate, psl, pel} each with >=3 genes detected.

Groups: acute (PRJNA1056765 BALF) | chronic (PRJEB24688 CF sputum) | healthy (PRJNA390194 non-COPD BALF)

Outputs: quadrant_table.csv, group_summary.json, and a printed report.
"""
import json, os, glob, csv, statistics as st, math
from math import comb

ABF = "/agent/workspace/abf"
os.chdir(ABF)

FAST, SLOW = 0.10, 0.05
MIN_TARGET = 100          # depth gate, same for every group
SYSTEMS = ["alginate", "psl", "pel"]

acute = {l.strip() for l in open("acute_runs.txt") if l.strip()}
healthy = {l.strip() for l in open("healthy_runs.txt") if l.strip()}


def group_of(srr):
    if srr in acute:
        return "acute"
    if srr in healthy:
        return "healthy"
    if srr.startswith("ERR"):
        return "chronic"
    return "?"


def growth_call(f):
    if f > FAST: return "FAST"
    if f < SLOW: return "SLOW"
    return "ambig"


def matrix_call(systems_on, systems):
    # >=2 of 3 systems with >=3 genes each
    if len(systems_on) >= 2: return "ON"
    ndet = sum(1 for s in SYSTEMS if systems[s]["genes_detected"] >= 1)
    if ndet >= 1: return "alg-only" if systems_on == ["alginate"] else "partial"
    return "OFF"


rows = []
for f in sorted(glob.glob("out/*_v3.json")):
    d = json.load(open(f))
    g = group_of(d["srr"])
    if g == "?": continue
    rows.append({
        "srr": d["srr"], "group": g,
        "target": d["target_reads"], "f_rp": d["f_rp"], "f_mat": d["f_mat"],
        "growth": growth_call(d["f_rp"]),
        "matrix": matrix_call(d["systems_on"], d["systems"]),
        "systems_on": ";".join(d["systems_on"]),
        "alg": d["systems"]["alginate"]["genes_detected"],
        "psl": d["systems"]["psl"]["genes_detected"],
        "pel": d["systems"]["pel"]["genes_detected"],
        "mapped": d["mapped_reads"],
    })

rows.sort(key=lambda r: (r["group"], -r["f_rp"]))

with open("quadrant_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# ---------------- report ------------------------------------------------
print(f"{'sample':13s}{'grp':8s}{'target':>8s}{'%RP':>8s}{'growth':>7s}{'matrix':>9s}"
      f"{'alg':>4s}{'psl':>4s}{'pel':>4s}{'systems':>18s}")
print("-" * 88)
for r in rows:
    print(f"{r['srr']:13s}{r['group']:8s}{r['target']:>8,}{r['f_rp']*100:>7.2f}%{r['growth']:>7s}"
          f"{r['matrix']:>9s}{r['alg']:>4}{r['psl']:>4}{r['pel']:>4}{r['systems_on']:>18s}")

# ---------------- group summaries --------------------------------------
summary = {}
print("\n" + "=" * 88)
for g in ("acute", "chronic", "healthy"):
    gr = [r for r in rows if r["group"] == g]
    gd = [r for r in gr if r["target"] >= MIN_TARGET]
    n = len(gd)
    frp = [r["f_rp"] for r in gd]
    nfast = sum(1 for r in gd if r["growth"] == "FAST")
    nslow = sum(1 for r in gd if r["growth"] == "SLOW")
    namb = n - nfast - nslow
    non = sum(1 for r in gd if r["matrix"] == "ON")
    summary[g] = {
        "n_total": len(gr), "n_depth_ok": n,
        "median_pctRP": round(100 * st.median(frp), 2) if frp else None,
        "pctRP_range": [round(100 * min(frp), 2), round(100 * max(frp), 2)] if frp else None,
        "FAST": nfast, "SLOW": nslow, "ambiguous": namb,
        "matrix_ON": non, "matrix_ON_rate": round(non / n, 3) if n else None,
        "matrix_alg_only": sum(1 for r in gd if r["matrix"] == "alg-only"),
        "matrix_OFF": sum(1 for r in gd if r["matrix"] == "OFF"),
    }
    print(f"\n### {g.upper()}   (depth-ok n={n} of {len(gr)}; target>= {MIN_TARGET})")
    if n:
        print(f"  growth  : median %RP {summary[g]['median_pctRP']:.2f}%  "
              f"range {summary[g]['pctRP_range']}  FAST {nfast}  ambiguous {namb}  SLOW {nslow}")
        print(f"  matrix  : ON {non}/{n} ({summary[g]['matrix_ON_rate']*100:.0f}%)  "
              f"alginate-only {summary[g]['matrix_alg_only']}  OFF {summary[g]['matrix_OFF']}")

# ---------------- quadrant tables --------------------------------------
print("\n" + "=" * 88)
print("QUADRANT TABLES  (rows = matrix, cols = growth)")
for g in ("acute", "chronic", "healthy"):
    gd = [r for r in rows if r["group"] == g and r["target"] >= MIN_TARGET]
    if not gd: continue
    print(f"\n{g.upper()} (n={len(gd)})")
    hdr = f"{'':12s}{'FAST':>8s}{'ambig':>8s}{'SLOW':>8s}"
    print(hdr)
    for mc in ("ON", "alg-only", "partial", "OFF"):
        cells = {gc: sum(1 for r in gd if r["matrix"] == mc and r["growth"] == gc)
                 for gc in ("FAST", "ambig", "SLOW")}
        if sum(cells.values()) == 0: continue
        print(f"{mc:12s}{cells['FAST']:>8d}{cells['ambig']:>8d}{cells['SLOW']:>8d}")

json.dump({"rows": rows, "summary": summary}, open("quadrant_summary.json", "w"), indent=1)
print("\nwrote quadrant_table.csv / quadrant_summary.json")
