#!/usr/bin/env python3
"""item5_depth_dilution.py — the experiment the literature says is missing.

Take ONE deep, Pseudomonas-dominant CF sputum run (Rossi ERR2275090, ~400k
target reads) and quantify it repeatedly at DECREASING subsampled depths.
If %RP drifts systematically as depth falls, the growth axis is load-sensitive
and a "SLOW" call at low depth may be an artifact rather than biology.

Also runs the same dilution on the Cornforth in-vivo sputum run for a second
organism/chemistry.

Output: item5/dilution_<SRR>.json  -> {depth: {target_reads, coding, rp, pctRP}}
"""
import os, json, subprocess, gzip, collections

os.chdir("/agent/workspace/abf")
os.makedirs("item5", exist_ok=True)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"

IV = json.load(open("gff_intervals3.json"))
iv = IV["pseu"]

DEPTHS = [2_000_000, 1_000_000, 500_000, 250_000, 100_000, 50_000, 20_000, 10_000]


def hit(lst, pos, width=300):
    lo, hi = 0, len(lst)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e = lst[mid][0], lst[mid][1]
        if pos < s - width: hi = mid
        elif pos > e + width: lo = mid + 1
        else: return lst[mid]
    return None


def quantify(fq, org="pseu"):
    t = rna = rp = dec = 0
    sysdet = collections.Counter()
    with os.popen(f"{MM2} -x sr -t 2 -a --secondary=no gcomp3_{org}.mmi {fq} 2>/dev/null") as fh:
        for line in fh:
            if line[0] == "@": continue
            f = line.split("\t", 6)
            if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20: continue
            if f[2].startswith("decoy|"): dec += 1; continue
            d = iv.get(f[2])
            if not d: continue
            t += 1
            pos = int(f[3])
            if hit(d["rna"], pos): rna += 1; continue
            m = hit(d["matrix"], pos)
            if m: sysdet[m[2]] += 1
            if hit(d["rp"], pos): rp += 1
    den = t - rna
    return {"target_reads": t, "rrna": rna, "coding": den, "rp": rp, "decoy": dec,
            "pctRP": round(100 * rp / den, 2) if den else None,
            "matrix_systems": dict(sysdet)}


SRCS = [("ERR2275090", "item4b/ERR2275090_1.fastq.gz"),
        ("SRR6833349", "item4/SRR6833349_1.fastq.gz")]

for srr, fq in SRCS:
    out = f"item5/dilution_{srr}.json"
    if os.path.exists(out):
        print("have", srr); continue
    if not os.path.exists(fq):
        print("nofile", fq); continue
    res = {}
    for d in DEPTHS:
        sub = f"item5/{srr}_{d}.fq"
        # take d reads (4 lines each)
        subprocess.run(f"zcat {fq} 2>/dev/null | head -n {4*d} > {sub}",
                       shell=True, timeout=900)
        if not os.path.exists(sub) or os.path.getsize(sub) == 0:
            continue
        r = quantify(sub)
        r["requested"] = d
        res[str(d)] = r
        os.remove(sub)
        print(f"{srr} depth={d:>9,} target={r['target_reads']:>7} coding={r['coding']:>7} "
              f"pctRP={r['pctRP']}", flush=True)
    json.dump(res, open(out, "w"), indent=1)
print("ITEM5 DONE", flush=True)
