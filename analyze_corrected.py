#!/usr/bin/env python3
"""CORRECTED quadrant analysis (growth axis fixed).

Growth = %RP_coding = RP reads / PAO1-mapping reads outside rRNA/tRNA.
Gifford anchors, recalibrated on this dataset's own lab cultures:
  lab exponential (fast) = 13.6%   lab stationary (slow) = 2.1%
  -> FAST >= 10%, SLOW <= 5% (unchanged thresholds; the metric scale is now right)

Matrix = systems detected among alginate/psl/pel (unchanged).
"""
import json, os, glob, csv, math, statistics as st

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

FAST, SLOW = 10.0, 5.0
MIN_TOT = 1000          # depth gate on total PAO1-mapping reads


def growth(p):
    if p is None: return "?"
    if p > FAST: return "FAST"
    if p < SLOW: return "SLOW"
    return "ambig"


def sys_on(d):
    return [s for s in ("alginate", "psl", "pel") if d["systems"][s]["genes_detected"] >= 3]


def mat_call(d):
    on = sys_on(d)
    if len(on) >= 2: return "ON"
    if len(on) == 1: return "alg-only" if on == ["alginate"] else "partial"
    nd = sum(1 for s in ("alginate", "psl", "pel") if d["systems"][s]["genes_detected"] >= 1)
    return "partial" if nd else "OFF"


rows = []
for s in acute + sorted(INVIVO) + sorted(EXP) + sorted(STAT) + healthy:
    fq = f"out/{s}_pseu_v3.json"; fr = f"out/{s}_rp2.json"
    if not (os.path.exists(fq) and os.path.exists(fr)):
        continue
    d = json.load(open(fq)); r = json.load(open(fr))
    grp = ("acute" if s in acute else "chronic-invivo" if s in INVIVO else
           "chronic-exp" if s in EXP else "chronic-stat" if s in STAT else "healthy")
    rows.append({"srr": s, "group": grp, "tot_pao1": r["total"], "coding": r["coding_reads"],
                 "pctRP_coding": r["pctRP_coding"], "growth": growth(r["pctRP_coding"]),
                 "matrix": mat_call(d), "systems_on": ";".join(sys_on(d)),
                 "alg": d["systems"]["alginate"]["genes_detected"],
                 "psl": d["systems"]["psl"]["genes_detected"],
                 "pel": d["systems"]["pel"]["genes_detected"]})

with open("quadrant_table_corrected.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print(f"{'sample':13s}{'group':16s}{'PAO1reads':>10s}{'%RP':>7s}{'growth':>7s}{'matrix':>9s}{'systems':>16s}")
print("-" * 80)
for r in sorted(rows, key=lambda r: (r["group"], -(r["pctRP_coding"] or 0))):
    pc = f"{r['pctRP_coding']:.1f}%" if r["pctRP_coding"] is not None else "--"
    print(f"{r['srr']:13s}{r['group']:16s}{r['tot_pao1']:>10,}{pc:>7s}{r['growth']:>7s}"
          f"{r['matrix']:>9s}{r['systems_on']:>16s}")

print("\n" + "=" * 80)
print("GROUP SUMMARY (growth gated on >= %d PAO1 reads; matrix on all)" % MIN_TOT)
for g in ("acute", "chronic-invivo", "chronic-exp", "chronic-stat", "healthy"):
    gr = [r for r in rows if r["group"] == g]
    gd = [r for r in gr if r["tot_pao1"] >= MIN_TOT and r["pctRP_coding"] is not None]
    v = [r["pctRP_coding"] for r in gd]
    print(f"  {g:16s} n={len(gr):2d} depth-ok={len(gd):2d}  median %RP "
          f"{st.median(v):5.1f}%" if v else f"  {g:16s} n={len(gr):2d} depth-ok=0", end="")
    if v:
        print(f"  FAST {sum(1 for x in v if x>FAST)}  ambig {sum(1 for x in v if SLOW<=x<=FAST)}"
              f"  SLOW {sum(1 for x in v if x<SLOW)}")
    else:
        print()

print("\nQUADRANT (matrix x growth), depth-ok samples")
for g in ("acute", "chronic-invivo"):
    gd = [r for r in rows if r["group"] == g and r["tot_pao1"] >= MIN_TOT]
    print(f"\n  {g} (n={len(gd)})")
    print(f"  {'':12s}{'FAST':>7s}{'ambig':>7s}{'SLOW':>7s}")
    for mc in ("ON", "alg-only", "partial", "OFF"):
        cells = {gc: sum(1 for r in gd if r["matrix"] == mc and r["growth"] == gc)
                 for gc in ("FAST", "ambig", "SLOW")}
        if sum(cells.values()):
            print(f"  {mc:12s}{cells['FAST']:>7d}{cells['ambig']:>7d}{cells['SLOW']:>7d}")

json.dump(rows, open("quadrant_corrected.json", "w"), indent=1)
print("\nwrote quadrant_table_corrected.csv / quadrant_corrected.json")
