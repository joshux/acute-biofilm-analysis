#!/usr/bin/env python3
"""Feasibility probe: full read-composition breakdown for a metatranscriptome.

Answers the question "can %RP be measured in this dataset at all?" by splitting
every MAPQ>=20 alignment into:
    PAO1 rRNA/tRNA | PAO1 ribosomal-protein | PAO1 other coding | decoy (non-Pseudomonas)
plus unmapped. Reports the rRNA fraction and the coding-read count that would
form the %RP denominator.

usage: probe_composition.py <SRR> <fastq>
"""
import sys, os, json, collections, gzip, re

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
MAPQ = 20

SRR, FASTQ = sys.argv[1], sys.argv[2]

RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")
rna_iv = collections.defaultdict(list)
rp_iv = collections.defaultdict(list)
with gzip.open("pao1.gff.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        seqid, typ, s, e = c[0], c[2], int(c[3]), int(c[4])
        if typ in ("rRNA", "tRNA", "tmRNA"):
            rna_iv[seqid].append((s, e))
        elif typ == "CDS":
            m = re.search(r"gene=([^;]+)", c[8]) or re.search(r"Name=([^;]+)", c[8])
            g = m.group(1) if m else ""
            if RP_RE.match(g):
                rp_iv[seqid].append((s, e))
for d in (rna_iv, rp_iv):
    for k in d:
        d[k].sort()


def overlaps(ivs, pos):
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
aln = pseu = pseu_rna = pseu_rp = decoy = 0
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < MAPQ:
            continue
        aln += 1
        ref = f[2]
        if ref.startswith("decoy|"):
            decoy += 1
            continue
        pseu += 1
        pos = int(f[3])
        if overlaps(rna_iv.get(ref, []), pos):
            pseu_rna += 1
        elif overlaps(rp_iv.get(ref, []), pos):
            pseu_rp += 1

coding = pseu - pseu_rna
out = {
    "srr": SRR,
    "aligned_mapq20": aln,
    "pseu_total": pseu,
    "pseu_rrna": pseu_rna,
    "pseu_rp": pseu_rp,
    "pseu_other_coding": coding - pseu_rp,
    "decoy": decoy,
    "coding_reads": coding,
    "pctRP_coding": round(100 * pseu_rp / coding, 2) if coding else None,
    "rrna_frac_of_pseu": round(100 * pseu_rna / pseu, 1) if pseu else None,
}
json.dump(out, open(f"out/{SRR}_probe.json", "w"), indent=1)
print(json.dumps(out))
