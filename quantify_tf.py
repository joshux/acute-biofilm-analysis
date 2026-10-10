#!/usr/bin/env python3
"""Independent growth proxy: %TF = translation/transcription-machinery share.

%TF = reads overlapping PAO1 cat=="growth" intervals (from gff/panels.json)
      MINUS ribosomal-protein genes (RP family excluded so the proxy is
      genuinely independent of the %RP metric)
      / (PAO1-mapping reads OUTSIDE rRNA/tRNA features)

The denominator is byte-for-byte the same as quantify_rp2.py (PAO1 reads
outside rRNA/tRNA), so %TF and %RP_coding are measured on the same base.

Mapping recipe identical to quantify_rp2.py:
  minimap2 -x sr -t 2 -a --secondary=no compfull_pseu.mmi <fastq>
  skip @ lines, flag&4, MAPQ<20, refs starting with "decoy|"
  rRNA/tRNA/tmRNA intervals from pao1.gff.gz excluded from denominator

usage: quantify_tf.py <SRR>
"""
import sys, os, json, collections, gzip, re

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
MAPQ = 20
SRR = sys.argv[1]
FASTQ = f"fastq/{SRR}.fastq.gz"

# --- RNA intervals (denominator) from the GFF, exactly as quantify_rp2.py ---
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")
rna_iv = collections.defaultdict(list)
n_rrna = 0
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

# --- TF intervals: panel cat=="growth" MINUS RP-family genes ----------------
panel = json.load(open("gff/panels.json"))["GCF_000006765.1"]
tf_iv = collections.defaultdict(list)
tf_genes, rp_excluded = [], []
for g in panel:
    if g.get("cat") != "growth":
        continue
    if RP_RE.match(g["gene"]):
        rp_excluded.append(g["gene"])
        continue
    tf_iv[g["seqid"]].append((int(g["start"]), int(g["end"])))
    tf_genes.append(g["gene"])
for k in tf_iv:
    tf_iv[k].sort()
n_tf = sum(len(v) for v in tf_iv.values())

if n_tf < 10:
    sys.stderr.write(f"STOP: only {n_tf} non-RP growth genes; proxy too thin\n")
    sys.exit(2)


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
tot = rna = tf = 0
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
        if overlaps(tf_iv.get(ref, []), pos):
            tf += 1

denom = tot - rna
out = {"srr": SRR, "total": tot, "rrna_reads": rna, "coding_reads": denom,
       "tf_reads": tf, "pctTF_coding": round(100 * tf / denom, 2) if denom else None,
       "rrna_frac": round(100 * rna / tot, 1) if tot else None,
       "n_rrna_feat": n_rrna, "n_tf_genes": n_tf,
       "tf_genes": sorted(tf_genes), "rp_excluded_from_tf": sorted(rp_excluded)}
json.dump(out, open(f"out/{SRR}_tf.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("tf_genes", "rp_excluded_from_tf")}))
