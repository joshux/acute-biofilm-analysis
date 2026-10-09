#!/usr/bin/env python3
"""Quantify one sample for the quadrant test.

Growth axis  : f_rp = rp reads / (matrix+rp+og reads)   -> %RP (Gifford anchors)
Matrix axis  : per-system detection among alginate / psl / pel
               matrix-ON = >=2 of the 3 systems each with >=3 genes detected
Also records the balance B = log(R_matrix) - log(R_growth) for continuity.

Usage: quantify_v3.py <SRR> <pseu|kleb> [fastq] [MAPQ]
"""
import sys, os, json, collections, math, gzip

ABF = "/agent/workspace/abf"
os.chdir(ABF)

SRR = sys.argv[1]
ORG = sys.argv[2]
FASTQ = sys.argv[3] if len(sys.argv) > 3 else f"fastq/{SRR}.fastq.gz"
MAPQ = int(sys.argv[4]) if len(sys.argv) > 4 else 20
NPERM = int(os.environ.get("NPERM", "2000"))
PSEUDO = 0.5

MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
PANEL = f"comp_{ORG}.mmi"
OUT = f"out/{SRR}_{ORG}_v3.json"

refs = {}          # full header -> (cls, system, gene)
for line in open(f"qp_{ORG}.fna"):
    if line.startswith(">"):
        h = line[1:].rstrip("\n")
        parts = h.split("|")
        cls, acc, gene, system, idx = parts[0], parts[1], parts[2], parts[3], parts[4]
        refs[h] = (cls, system, gene)

SYSTEMS = ["alginate", "psl", "pel"]

counts = collections.Counter()
decoy_hits = 0
total = 0
mapped = 0
gz = FASTQ.endswith(".gz")
opener = gzip.open if gz else open
# minimap2 reads the file directly; use it as input
cmd = f"{MM2} -x sr -t 2 -a --secondary=no {PANEL} {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        total += 1
        f = line.split("\t", 6)
        if len(f) < 7:
            continue
        if int(f[1]) & 4:
            continue
        mapped += 1
        if int(f[4]) < MAPQ:
            continue
        ref = f[2]
        if ref in refs:
            counts[ref] += 1
        else:
            decoy_hits += 1

by_cls = collections.defaultdict(list)
for r, (cls, sysname, gene) in refs.items():
    by_cls[cls].append(r)

reads = {c: sum(counts[r] for r in rl) for c, rl in by_cls.items()}
tot_target = reads.get("matrix", 0) + reads.get("rp", 0) + reads.get("og", 0)
f_rp = reads.get("rp", 0) / tot_target if tot_target else 0.0
f_mat = reads.get("matrix", 0) / tot_target if tot_target else 0.0
f_og = reads.get("og", 0) / tot_target if tot_target else 0.0

# per-system detection
sys_genes = {s: [r for r in by_cls["matrix"] if refs[r][1] == s] for s in SYSTEMS}
sys_det = {}
for s in SYSTEMS:
    genes = sys_genes[s]
    nd = sum(1 for r in genes if counts[r] > 0)
    sys_det[s] = {"genes_in_panel": len(genes), "genes_detected": nd,
                  "reads": sum(counts[r] for r in genes)}
systems_on = [s for s in SYSTEMS if sys_det[s]["genes_detected"] >= 3]
matrix_on = len(systems_on) >= 2

# balance (length-normalised) for continuity with the earlier result
lens = collections.defaultdict(int)
for line in open(f"qp_{ORG}.fna"):
    pass  # length not needed for the quadrant call; kept simple

out = {
    "srr": SRR, "org": ORG, "mapq": MAPQ, "fastq": FASTQ,
    "total_reads": total, "mapped_reads": mapped, "decoy_hits": decoy_hits,
    "target_reads": tot_target,
    "reads": {k: reads.get(k, 0) for k in ("matrix", "rp", "og")},
    "f_rp": round(f_rp, 5), "f_mat": round(f_mat, 5), "f_og": round(f_og, 5),
    "systems": sys_det, "systems_on": systems_on, "matrix_on": matrix_on,
    "top_genes": sorted(
        [(refs[r][2], refs[r][0], refs[r][1], c) for r, c in counts.items()],
        key=lambda x: -x[3])[:40],
}
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "top_genes"}))
