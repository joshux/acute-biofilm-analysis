#!/usr/bin/env python3
"""Final quantification: matrix vs growth transcript investment per sample,
with a decoy-competition index and explicit metric variants.

Usage: quantify3.py <SRR> <pseu|kleb> [MAPQ]

Panels (single shared mapping run over comp_<org>.mmi = 34 decoy genera + target):
  matrix = matrix / matrix_prod genes (alg/pel/psl/csg/bcs/ica/cdr + EPS/capsule)
  rp     = ribosomal-protein genes        <- canonical growth marker (Gifford 2014)
  og     = other growth genes (replication/division/transcription/translation factors)

Metrics (all length-normalised as reads per kb of class panel):
  R_cls  = reads_cls / (sum of cls gene lengths / 1000)
  B_rp   = log(R_matrix / R_rp)      <- PRIMARY (both classes canonical)
  B_full = log(R_matrix / R_growth)  <- SENSITIVITY (growth = rp + og)
  f_cls  = reads_cls / total target reads   (compositional fractions)
Permutation CI resamples the observed reads across panel genes (multinomial)
and recomputes B, giving a read-level null for the class ratio.
"""
import sys, os, json, collections, math, random

ABF = "/tmp/abf"
os.chdir(ABF)

SRR = sys.argv[1]
ORG = sys.argv[2]
MAPQ = int(sys.argv[3]) if len(sys.argv) > 3 else 20
NPERM = int(os.environ.get("NPERM", "3000"))
PSEUDO = 0.5
random.seed(20261009)

MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
PANEL = f"comp_{ORG}.mmi"
FASTQ = f"fastq/{SRR}.fastq.gz"

refs = {}
with open(f"qp_{ORG}.fna") as f:
    h = None; L = 0
    for line in f:
        if line.startswith(">"):
            if h is not None:
                cls, acc, gene, i = h.split("|")
                refs[h] = (cls, acc, gene, L)
            h = line[1:].rstrip("\n"); L = 0
        else:
            L += len(line.rstrip("\n"))
    if h is not None:
        cls, acc, gene, i = h.split("|")
        refs[h] = (cls, acc, gene, L)

cls_of = {k: v[0] for k, v in refs.items()}
len_of = {k: v[3] for k, v in refs.items()}
gene_of = {k: v[2] for k, v in refs.items()}

counts = collections.Counter()
decoy_hits = 0
total = 0
nmap_all = 0
cmd = f"{MM2} -x sr -t 4 -a --secondary=no {PANEL} {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        total += 1
        f = line.split("\t")
        if len(f) < 11:
            continue
        if int(f[1]) & 4:
            continue
        if int(f[4]) < MAPQ:
            continue
        nmap_all += 1
        if "|" in f[2]:
            counts[f[2]] += 1
        else:
            decoy_hits += 1

by_cls = collections.defaultdict(list)
for ref in refs:
    by_cls[cls_of[ref]].append(ref)
matrix_refs = by_cls["matrix"]
rp_refs = by_cls["rp"]
og_refs = by_cls["og"]
growth_refs = rp_refs + og_refs

kb = {c: sum(len_of[r] for r in rl) / 1000.0 for c, rl in
      (("matrix", matrix_refs), ("rp", rp_refs), ("og", og_refs), ("growth", growth_refs))}
reads = {c: sum(counts[r] for r in rl) for c, rl in
         (("matrix", matrix_refs), ("rp", rp_refs), ("og", og_refs), ("growth", growth_refs))}
R = {c: (reads[c] + PSEUDO * len(by_cls[c if c != "growth" else "rp"])
         + (PSEUDO * len(og_refs) if c == "growth" else 0)) / kb[c] for c in reads}
R["matrix"] = (reads["matrix"] + PSEUDO * len(matrix_refs)) / kb["matrix"]
R["rp"] = (reads["rp"] + PSEUDO * len(rp_refs)) / kb["rp"]
R["og"] = (reads["og"] + PSEUDO * len(og_refs)) / kb["og"]
R["growth"] = (reads["growth"] + PSEUDO * len(growth_refs)) / kb["growth"]

B_rp = math.log(R["matrix"]) - math.log(R["rp"])
B_full = math.log(R["matrix"]) - math.log(R["growth"])

tot_target = reads["matrix"] + reads["rp"] + reads["og"]
f = {c: reads[c] / tot_target if tot_target else None for c in ("matrix", "rp", "og")}

# permutation null: resample reads across panel genes
read_list = []
for ref, c in counts.items():
    read_list.extend([ref] * c)
nreads = len(read_list)
perm_rp, perm_full = [], []
if nreads > 0:
    for _ in range(NPERM):
        cnt = collections.Counter(random.choices(read_list, k=nreads))
        m = (sum(cnt[r] for r in matrix_refs) + PSEUDO * len(matrix_refs)) / kb["matrix"]
        g1 = (sum(cnt[r] for r in rp_refs) + PSEUDO * len(rp_refs)) / kb["rp"]
        g2 = (sum(cnt[r] for r in growth_refs) + PSEUDO * len(growth_refs)) / kb["growth"]
        perm_rp.append(math.log(m) - math.log(g1))
        perm_full.append(math.log(m) - math.log(g2))
    perm_rp.sort(); perm_full.sort()

def ci(v, p):
    if not v: return None
    return round(v[min(len(v) - 1, int(p * len(v)))], 4)

detect = {
    "matrix": round(sum(1 for r in matrix_refs if counts[r] > 0) / len(matrix_refs), 4),
    "rp": round(sum(1 for r in rp_refs if counts[r] > 0) / len(rp_refs), 4),
    "og": round(sum(1 for r in og_refs if counts[r] > 0) / len(og_refs), 4),
}

gene_tab = sorted(
    [(gene_of[r], cls_of[r], refs[r][1], counts[r], len_of[r]) for r in refs if counts[r] > 0],
    key=lambda x: -x[3])

out = {
    "srr": SRR, "org": ORG, "mapq": MAPQ,
    "total_reads": total, "mapped_reads": nmap_all, "decoy_hits": decoy_hits,
    "pct_mapped": round(100.0 * nmap_all / total, 4) if total else 0,
    "target_reads": tot_target,
    "reads": reads, "kb": {k: round(v, 2) for k, v in kb.items()},
    "R": {k: round(v, 4) for k, v in R.items()},
    "fractions": {k: (round(v, 5) if v is not None else None) for k, v in f.items()},
    "B_rp": round(B_rp, 4), "B_rp_lo": ci(perm_rp, 0.025), "B_rp_hi": ci(perm_rp, 0.975),
    "B_full": round(B_full, 4), "B_full_lo": ci(perm_full, 0.025), "B_full_hi": ci(perm_full, 0.975),
    "detect": detect, "n_panel_genes": {"matrix": len(matrix_refs), "rp": len(rp_refs), "og": len(og_refs)},
    "top_genes": gene_tab[:60],
}
json.dump(out, open(f"{SRR}_{ORG}_q3.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "top_genes"}))
