#!/usr/bin/env python3
"""Gifford-faithful %RP: ribosomal-protein reads / TOTAL target-organism reads.

The panel-based metric (RP / (rp+og+matrix)) has two flaws:
  1. multi-genome competition drops conserved RP reads at MAPQ>=20 (median MAPQ 0)
  2. the denominator is a small curated panel, not the whole transcriptome
Both bias %RP. This computes the real thing:

  index   = full PAO1 genome + 34 non-Pseudomonas decoy genomes
  total   = reads mapping to PAO1 (MAPQ>=20)
  rp      = those overlapping an annotated ribosomal-protein gene interval
  %RP     = rp / total          <- Gifford 2014 definition

The single Pseudomonas genome removes intra-species competition, so conserved
Pseudomonas RP reads are retained; the decoys keep non-Pseudomonas conserved
reads out.

usage: quantify_rp.py <SRR> [fastq]
"""
import sys, os, json, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
MAPQ = 20

SRR = sys.argv[1]
FASTQ = sys.argv[2] if len(sys.argv) > 2 else f"fastq/{SRR}.fastq.gz"

# --- RP intervals on PAO1 --------------------------------------------------
import re
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")
panels = json.load(open("gff/panels.json"))
pao1 = "GCF_000006765.1"
rp_iv = collections.defaultdict(list)   # seqid -> [(start,end)]
n_rp = 0
for g in panels[pao1]:
    if g["cat"].startswith("growth") and RP_RE.match(g["gene"]):
        rp_iv[g["seqid"]].append((g["start"], g["end"]))
        n_rp += 1
for k in rp_iv:
    rp_iv[k].sort()
print(f"PAO1 RP gene intervals: {n_rp}", file=sys.stderr)

# --- map ------------------------------------------------------------------
cmd = f"{MM2} -x sr -t 2 -a --secondary=no compfull_pseu.mmi {FASTQ} 2>/dev/null"
total = 0
rp = 0
per_contig = collections.Counter()
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < MAPQ:
            continue
        ref = f[2]
        if ref.startswith("decoy|"):
            continue                      # non-Pseudomonas -> not target
        total += 1
        per_contig[ref] += 1
        # alignment position
        pos = int(f[3])
        ivs = rp_iv.get(ref)
        if not ivs:
            continue
        # binary search for overlap
        lo, hi = 0, len(ivs)
        hit = False
        while lo < hi:
            mid = (lo + hi) // 2
            s, e = ivs[mid]
            if pos < s - 300:
                hi = mid
            elif pos > e + 300:
                lo = mid + 1
            else:
                hit = True
                break
        if hit:
            rp += 1

out = {"srr": SRR, "total_pao1": total, "rp_reads": rp,
       "pctRP": round(100 * rp / total, 2) if total else None,
       "contigs": dict(per_contig)}
json.dump(out, open(f"out/{SRR}_rp.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "contigs"}))
