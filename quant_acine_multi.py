#!/usr/bin/env python3
"""Acinetobacter multi-strain quantification (genus-level, read-deduplicated).

Why: clinical A. baumannii strains diverge from ATCC 17978, so a single-strain
index loses most reads (pilot + this run: median ~6 coding reads). But a
multi-strain index collapses MAPQ for conserved genes — so instead of a MAPQ
filter we take each read's PRIMARY alignment once (minimap2 --secondary=no,
MAPQ>=0) and classify by the gene it lands on. Non-Acinetobacter reads are
absorbed by the 34 decoy genomes when the decoy is the better hit.

Axes:
  growth : %RP = RP reads / coding reads  (rRNA/tRNA excluded from denominator)
  matrix : pgaABCD (PNAG) + csuABCDE (Csu pili) gene detection,
           systems = {pga, csu}; matrix-ON = both detected (>=2 genes each)

usage: quant_acine_multi.py <SRR>
"""
import os, sys, json, gzip, re, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"

STRAINS = {
    "GCF_026167805.1": "gff_acine.gz",        # 17978
    "GCF_000241685.1": "gff_ab5075.gz",       # AB5075
    "GCF_057224575.1": "gff_nosocomialis.gz",
    "GCF_055475515.1": "gff_pittii.gz",
    "GCF_988236795.1": "gff_johnsonii.gz",
    "GCF_050703235.1": "gff_lwoffii.gz",
}
MATRIX = {"pgaA", "pgaB", "pgaC", "pgaD", "csuA", "csuAB", "csuB", "csuC", "csuD", "csuE"}
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")

# intervals: ref in index is ">{acc}|{cid}" so we key by (acc, cid)
rna_iv = collections.defaultdict(list)
rp_iv = collections.defaultdict(list)
mat_iv = collections.defaultdict(list)   # (acc,cid) -> [(s,e,gene)]
for acc, gff in STRAINS.items():
    with gzip.open(gff, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            cid, typ, s, e = c[0], c[2], int(c[3]), int(c[4])
            key = (acc, cid)
            if typ in ("rRNA", "tRNA", "tmRNA"):
                rna_iv[key].append((s, e))
            elif typ == "CDS":
                m = re.search(r"gene=([^;]+)", c[8]) or re.search(r"Name=([^;]+)", c[8])
                g = m.group(1) if m else ""
                if g in MATRIX:
                    mat_iv[key].append((s, e, g))
                elif RP_RE.match(g):
                    rp_iv[key].append((s, e))


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


SRR = sys.argv[1]
FASTQ = f"fastq/{SRR}.fastq.gz"
os.makedirs("aci_out2", exist_ok=True)
OUT = f"aci_out2/{SRR}.json"
if os.path.exists(OUT):
    print("have", SRR); sys.exit(0)

t = rna = rp = dec = 0
mat_genes = collections.Counter()
cmd = f"{MM2} -x sr -t 2 -a --secondary=no acine_multi.mmi {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4:
            continue
        ref = f[2]
        if ref.startswith("decoy|"):
            dec += 1
            continue
        acc, cid = ref.split("|", 1)
        key = (acc, cid)
        t += 1
        pos = int(f[3])
        if hit(rna_iv.get(key, []), pos):
            rna += 1
            continue
        m = hit(mat_iv.get(key, []), pos)
        if m:
            mat_genes[m[2]] += 1
        if hit(rp_iv.get(key, []), pos):
            rp += 1

den = t - rna
pga_det = sum(1 for g in ("pgaA", "pgaB", "pgaC", "pgaD") if mat_genes[g] > 0)
csu_det = sum(1 for g in ("csuA", "csuAB", "csuB", "csuC", "csuD", "csuE") if mat_genes[g] > 0)
matrix_on = (pga_det >= 2) and (csu_det >= 2)

out = {"srr": SRR, "org": "Acinetobacter (multi-strain)", "total": t, "rrna": rna,
       "coding": den, "rp": rp, "decoy": dec,
       "pctRP": round(100 * rp / den, 2) if den else None,
       "pga_genes": pga_det, "csu_genes": csu_det,
       "matrix_on": matrix_on, "matrix_genes": dict(mat_genes),
       "matrix_reads": sum(mat_genes.values())}
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "matrix_genes"}))
