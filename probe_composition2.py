#!/usr/bin/env python3
"""Feasibility probe v2 — full contingency table, no MAPQ pre-filter.

Reports, for each sample:
  alignments by MAPQ band (0 / 1-19 / >=20)
  x target (PAO1) class (rRNA-tRNA / ribosomal-protein / other coding) vs decoy

The point is to see how much of the library is bacterial at all, what the rRNA
burden is, and how much target-organism coding signal survives at MAPQ>=20.

usage: probe_composition2.py <SRR> <fastq>
"""
import sys, os, json, collections, gzip, re

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"

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


def band(q):
    return "0" if q == 0 else ("1-19" if q < 20 else "20+")


T = collections.Counter()          # (band, class) -> n
n_rec = 0
cmd = f"{MM2} -x sr -t 2 -a --secondary=no compfull_pseu.mmi {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7:
            continue
        n_rec += 1
        if int(f[1]) & 4:
            T[(band(0), "unmapped")] += 1
            continue
        b = band(int(f[4]))
        ref, pos = f[2], int(f[3])
        if ref.startswith("decoy|"):
            T[(b, "decoy")] += 1
        elif overlaps(rna_iv.get(ref, []), pos):
            T[(b, "PAO1_rRNA_tRNA")] += 1
        elif overlaps(rp_iv.get(ref, []), pos):
            T[(b, "PAO1_ribosomal_protein")] += 1
        else:
            T[(b, "PAO1_other_CDS")] += 1

bands = ["20+", "1-19", "0"]
classes = ["PAO1_ribosomal_protein", "PAO1_other_CDS", "PAO1_rRNA_tRNA", "decoy", "unmapped"]
out = {"srr": SRR, "records": n_rec}
print(f"=== {SRR}  ({n_rec:,} reads) ===")
print(f"{'class':26s}" + "".join(f"{b:>12s}" for b in bands) + f"{'total':>12s}")
for c in classes:
    row = [T[(b, c)] for b in bands]
    out[c] = dict(zip(bands, row))
    print(f"{c:26s}" + "".join(f"{v:>12,d}" for v in row) + f"{sum(row):>12,d}")

# headline numbers
t20 = T[("20+", "PAO1_ribosomal_protein")]
c20 = T[("20+", "PAO1_other_CDS")]
r20 = T[("20+", "PAO1_rRNA_tRNA")]
out["pctRP_coding_mapq20"] = round(100 * t20 / (t20 + c20), 2) if (t20 + c20) else None
out["pctRP_genome_mapq20"] = round(100 * t20 / (t20 + c20 + r20), 2) if (t20 + c20 + r20) else None
out["coding_reads_mapq20"] = t20 + c20
# permissive: any PAO1 hit
tp = sum(T[(b, "PAO1_ribosomal_protein")] for b in bands)
cp = sum(T[(b, "PAO1_other_CDS")] for b in bands)
rp_ = sum(T[(b, "PAO1_rRNA_tRNA")] for b in bands)
out["pctRP_permissive"] = round(100 * tp / (tp + cp), 2) if (tp + cp) else None
out["coding_reads_permissive"] = tp + cp
out["pao1_any"] = tp + cp + rp_
out["decoy_any"] = sum(T[(b, "decoy")] for b in bands)
json.dump(out, open(f"out/{SRR}_probe2.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if not isinstance(v, dict)}))
