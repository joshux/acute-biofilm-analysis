#!/usr/bin/env python3
"""Quantify the rRNA burden: what share of PAO1-mapping reads fall in the
highest-coverage windows (rRNA operons are the classic extreme-coverage peaks)?

If acute samples carry far more rRNA than chronic (Tang depleted human rRNA only;
Rossi depleted bacterial rRNA too), then acute's Gifford denominator is inflated
and its %RP is an UNDER-estimate -- i.e. the acute-fast finding is conservative.
"""
import sys, os, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
WIN = 5000


def run(srr):
    fastq = f"fastq/{srr}.fastq.gz"
    cov = collections.Counter()
    total = 0
    cmd = f"{MM2} -x sr -t 2 -a --secondary=no compfull_pseu.mmi {fastq} 2>/dev/null"
    with os.popen(cmd) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.split("\t", 6)
            if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20:
                continue
            if f[2].startswith("decoy|"):
                continue
            total += 1
            cov[int(f[3]) // WIN] += 1
    if not total:
        return None
    ranked = sorted(cov.items(), key=lambda x: -x[1])
    top1 = sum(c for _, c in ranked[:1]) / total
    top5 = sum(c for _, c in ranked[:5]) / total
    top20 = sum(c for _, c in ranked[:20]) / total
    return {"srr": srr, "total": total, "n_windows": len(cov),
            "top1_share": round(100 * top1, 1),
            "top5_share": round(100 * top5, 1),
            "top20_share": round(100 * top20, 1),
            "top_window_reads": ranked[0][1]}


for srr in ["SRR27343249", "SRR27343409", "ERR2275089", "ERR2275092", "ERR2591781", "ERR2275086"]:
    r = run(srr)
    print(r, flush=True)
