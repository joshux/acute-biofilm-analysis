#!/usr/bin/env python3
"""Quantify Acinetobacter-dominant acute samples.

Matrix axis : pgaABCD (PNAG) + csuABCDE (Csu pili) gene detection
Growth axis : %RP = RP reads / coding reads (Gifford-faithful, rRNA excluded)
              against the Acinetobacter genome + decoys
"""
import os, sys, json, gzip, re, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
ACC = "GCF_026167805.1"
MATRIX = {"pgaA", "pgaB", "pgaC", "pgaD", "csuA", "csuAB", "csuB", "csuC", "csuD", "csuE"}

SRR = sys.argv[1]
FASTQ = f"fastq/{SRR}.fastq.gz"
os.makedirs("aci_out", exist_ok=True)
OUT = f"aci_out/{SRR}.json"
if os.path.exists(OUT):
    print("have", SRR); sys.exit(0)

# intervals
rna_iv = collections.defaultdict(list); rp_iv = collections.defaultdict(list)
mat_iv = collections.defaultdict(list)
gene_of = {}
with gzip.open("gff_acine.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        seqid, typ, s, e = c[0], c[2], int(c[3]), int(c[4])
        if typ in ("rRNA", "tRNA", "tmRNA"):
            rna_iv[seqid].append((s, e)); continue
        if typ != "CDS":
            continue
        m = re.search(r"gene=([^;]+)", c[8]) or re.search(r"Name=([^;]+)", c[8])
        g = m.group(1) if m else ""
        if g in MATRIX:
            mat_iv[seqid].append((s, e, g))
        elif re.match(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$", g):
            rp_iv[seqid].append((s, e))


def hit(ivs, pos):
    lo, hi = 0, len(ivs)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid][0], ivs[mid][1]
        if pos < s - 300:
            hi = mid
        elif pos > e + 300:
            lo = mid + 1
        else:
            return ivs[mid]
    return None


t = rna = rp = dec = 0
mat_genes = collections.Counter()
cmd = f"{MM2} -x sr -t 2 -a --secondary=no rep_acine2.mmi {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20:
            continue
        if f[2].startswith("decoy|"):
            dec += 1; continue
        t += 1
        pos = int(f[3]); ref = f[2]
        if hit(rna_iv.get(ref, []), pos):
            rna += 1; continue
        m = hit(mat_iv.get(ref, []), pos)
        if m:
            mat_genes[m[2]] += 1
        if hit(rp_iv.get(ref, []), pos):
            rp += 1

den = t - rna
out = {"srr": SRR, "org": "Acinetobacter", "total": t, "rrna": rna, "coding": den,
       "rp": rp, "decoy": dec,
       "pctRP": round(100 * rp / den, 2) if den else None,
       "matrix_genes": dict(mat_genes),
       "matrix_genes_detected": len(mat_genes),
       "matrix_reads": sum(mat_genes.values())}
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps(out))
