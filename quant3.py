#!/usr/bin/env python3
"""quant3.py — two-axis quantification with the CORRECTED per-organism panels.

  growth %RP = RP reads / coding reads   (rRNA/tRNA excluded)
  matrix: named systems, each needing >= MIN_GENES genes detected
  matrix_ON_EPS  = >=2 EPS systems           (polysaccharide matrix)
  matrix_ON_PROT = >=2 PROTEIN systems       (proteinaceous fimbriae/adhesin)
  matrix_ON_any  = >=2 systems of either type

usage: quant3.py <SRR> <fastq_dir> <out_dir> [subsample_reads]
"""
import os, sys, json, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
MIN_GENUS_READS = 20
MIN_GENES = 2

GENUS2ORG = {"pseudomonas": "pseu", "klebsiella": "kleb", "escherichia": "ecoli",
             "staphylococcus": "staph", "acinetobacter": "acine",
             "haemophilus": "haemo", "streptococcus": "strep",
             "stenotrophomonas": "steno"}

SRR = sys.argv[1]
FDIR = sys.argv[2] if len(sys.argv) > 2 else "sub"
ODIR = sys.argv[3] if len(sys.argv) > 3 else "final3"
os.makedirs(ODIR, exist_ok=True)
OUT = f"{ODIR}/{SRR}.json"
if os.path.exists(OUT):
    print("have", SRR); sys.exit(0)

FQ1 = f"{FDIR}/{SRR}_1.fastq.gz"; FQ2 = f"{FDIR}/{SRR}_2.fastq.gz"
FASTQ = f"{FDIR}/{SRR}.fastq.gz"
if os.path.exists(FQ1):
    FASTQ = f"{FQ1} {FQ2}" if os.path.exists(FQ2) else FQ1
elif not os.path.exists(FASTQ):
    print("nofile", SRR); sys.exit(1)

IV = json.load(open("gff_intervals3.json"))
g2 = json.load(open("acc2genus.json")) if os.path.exists("acc2genus.json") else {}

# ---- dominant genus -------------------------------------------------------
c = collections.Counter()
cmd = f"{MM2} -x sr -t 2 -a --secondary=no panel_multi3.mmi {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@": continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20: continue
        if f[2].startswith("decoy|"): continue
        acc = f[2].split("|")[0]
        c[acc] += 1

TARGETS = {"GCF_000006765.1": "pseudomonas", "GCF_000240185.1": "klebsiella",
           "GCF_000005845.2": "escherichia", "GCF_000013425.1": "staphylococcus",
           "GCF_026167805.1": "acinetobacter", "GCF_000012185.1": "haemophilus",
           "GCF_000006885.1": "streptococcus", "GCF_000072485.1": "stenotrophomonas"}
gcount = collections.Counter()
for acc, n in c.items():
    gcount[TARGETS.get(acc, acc)] += n

dom = gcount.most_common(1)[0][0] if gcount else None
org = GENUS2ORG.get(dom)
if not org or gcount[dom] < MIN_GENUS_READS:
    json.dump({"srr": SRR, "dominant": dom, "org": None, "pctRP": None,
               "matrix_on": None, "note": "no dominant target organism",
               "genus_counts": dict(gcount)}, open(OUT, "w"), indent=1)
    print(json.dumps({"srr": SRR, "dominant": dom, "org": None})); sys.exit(0)

iv = IV.get(org, {})

def hit(lst, pos, width=300):
    lo, hi = 0, len(lst)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e = lst[mid][0], lst[mid][1]
        if pos < s - width: hi = mid
        elif pos > e + width: lo = mid + 1
        else: return lst[mid]
    return None

t = rna = rp = dec = 0
sysdet = collections.Counter()
cmd = f"{MM2} -x sr -t 2 -a --secondary=no gcomp3_{org}.mmi {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@": continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20: continue
        ref = f[2]
        if ref.startswith("decoy|"): dec += 1; continue
        d = iv.get(ref)
        if not d: continue
        t += 1
        pos = int(f[3])
        if hit(d["rna"], pos): rna += 1; continue
        m = hit(d["matrix"], pos)
        if m: sysdet[m[2]] += 1
        if hit(d["rp"], pos): rp += 1

# system types from the panel census
P3 = json.load(open("panels3.json"))
TYPES = {}
for o, v in P3.items():
    for s in v["systems"]:
        nm, ty = s.split("("); TYPES[nm] = ty.rstrip(")")

den = t - rna
on = sorted([s for s, n in sysdet.items() if n >= MIN_GENES])
eps_on = [s for s in on if TYPES.get(s) == "EPS"]
prot_on = [s for s in on if TYPES.get(s) == "PROTEIN"]
out = {"srr": SRR, "dominant": dom, "org": org,
       "target_reads": t, "rrna": rna, "coding": den, "rp": rp, "decoy": dec,
       "pctRP": round(100 * rp / den, 2) if den else None,
       "matrix_systems": dict(sysdet), "systems_on": on,
       "eps_on": eps_on, "prot_on": prot_on,
       "matrix_on": len(on) >= 2,
       "matrix_on_eps": len(eps_on) >= 2,
       "genus_counts": dict(gcount)}
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "genus_counts"}))
