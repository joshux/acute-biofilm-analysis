#!/usr/bin/env python3
"""rRNA-excluded growth axis.

%RP_rna = RP reads / (all PAO1-mapping reads OUTSIDE rRNA/tRNA features)

Removing rRNA/tRNA from the denominator makes the metric independent of the
rRNA-depletion difference between the two datasets (Tang: human-only depletion;
Rossi: bacterial rRNA depleted). This is the cleanest available cross-dataset
comparison.

usage: quantify_rp2.py <SRR>
"""
import sys, os, json, collections, gzip

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
MAPQ = 20
SRR = sys.argv[1]
FASTQ = f"fastq/{SRR}.fastq.gz"

# --- load PAO1 features from the GFF --------------------------------------
import re
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")
rna_iv = collections.defaultdict(list)     # seqid -> [(s,e)]
rp_iv = collections.defaultdict(list)
n_rrna = n_rp = 0
with gzip.open("pao1.gff.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        seqid, typ, s, e = c[0], c[2], int(c[3]), int(c[4])
        if typ in ("rRNA", "tRNA", "tmRNA"):
            rna_iv[seqid].append((s, e)); n_rrna += 1
        elif typ == "CDS":
            m = re.search(r"gene=([^;]+)", c[8]) or re.search(r"Name=([^;]+)", c[8])
            g = m.group(1) if m else ""
            if RP_RE.match(g):
                rp_iv[seqid].append((s, e)); n_rp += 1


def overlaps(ivs, pos):
    if not ivs:
        return False
    lo, hi = 0, len(ivs)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid]
        if pos < s - 300:
            hi = mid
        elif pos > e + 300:
            lo = mid + 1
        else:
            return True
    return False


cmd = f"{MM2} -x sr -t 2 -a --secondary=no compfull_pseu.mmi {FASTQ} 2>/dev/null"
tot = rna = rp = 0
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < MAPQ:
            continue
        if f[2].startswith("decoy|"):
            continue
        tot += 1
        pos = int(f[3]); ref = f[2]
        if overlaps(rna_iv.get(ref, []), pos):
            rna += 1
            continue
        if overlaps(rp_iv.get(ref, []), pos):
            rp += 1

denom = tot - rna
out = {"srr": SRR, "total": tot, "rrna_reads": rna, "coding_reads": denom,
       "rp_reads": rp, "pctRP_coding": round(100 * rp / denom, 2) if denom else None,
       "rrna_frac": round(100 * rna / tot, 1) if tot else None,
       "n_rrna_feat": n_rrna, "n_rp_genes": n_rp}
json.dump(out, open(f"out/{SRR}_rp2.json", "w"), indent=1)
print(json.dumps(out))
