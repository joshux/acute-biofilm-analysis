#!/usr/bin/env python3
"""Per-organism Gifford-faithful two-axis quantification.

  growth  %RP      = RP reads / coding reads (rRNA/tRNA excluded)
  matrix  ON       = >=2 matrix subsystems each with >=2 genes detected

Each sample is scored against its own dominant organism's full-genome index
(gcomp_<org>.mmi = full genome + other-genus decoys).

usage: quant_final.py <SRR> [fastq_dir] [out_dir]
"""
import os, sys, json, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
MIN_GENUS_READS = 20

GENUS2ORG = {"pseudomonas": "pseu", "klebsiella": "kleb", "escherichia": "ecoli",
             "staphylococcus": "staph", "acinetobacter": "acine",
             "haemophilus": "haemo", "streptococcus": "strep",
             "stenotrophomonas": "steno"}

SRR = sys.argv[1]
FDIR = sys.argv[2] if len(sys.argv) > 2 else "fastq"
ODIR = sys.argv[3] if len(sys.argv) > 3 else "final_out"
os.makedirs(ODIR, exist_ok=True)
OUT = f"{ODIR}/{SRR}.json"
if os.path.exists(OUT):
    print("have", SRR); sys.exit(0)
FASTQ = f"{FDIR}/{SRR}.fastq.gz"
if not os.path.exists(FASTQ):
    print("nofile", SRR); sys.exit(1)

IV = json.load(open("gff_intervals.json"))
g2 = json.load(open("acc2genus.json"))

# --- dominant organism ----------------------------------------------------
c = collections.Counter()
cmd = f"{MM2} -x sr -t 2 -a --secondary=no panel_multi.mmi {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20:
            continue
        if f[2].startswith("decoy|"):
            continue
        c[g2.get(f[2].split("|")[0], "?")] += 1

dom = c.most_common(1)[0][0] if c else None
org = GENUS2ORG.get(dom)
if not org or c[dom] < MIN_GENUS_READS:
    json.dump({"srr": SRR, "dominant": dom, "org": None, "pctRP": None,
               "matrix_on": None, "note": "no dominant target organism"},
              open(OUT, "w"), indent=1)
    print(json.dumps({"srr": SRR, "dominant": dom, "org": None}))
    sys.exit(0)

iv = IV.get(org, {})


def hit(lst, pos, width=300):
    lo, hi = 0, len(lst)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e = lst[mid][0], lst[mid][1]
        if pos < s - width:
            hi = mid
        elif pos > e + width:
            lo = mid + 1
        else:
            return lst[mid]
    return None


t = rna = rp = dec = 0
sysdet = collections.Counter()
cmd = f"{MM2} -x sr -t 2 -a --secondary=no gcomp_{org}.mmi {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20:
            continue
        ref = f[2]
        if ref.startswith("decoy|"):
            dec += 1; continue
        d = iv.get(ref)
        if not d:
            continue
        t += 1
        pos = int(f[3])
        if hit(d["rna"], pos):
            rna += 1; continue
        m = hit(d["matrix"], pos)
        if m:
            sysdet[m[2]] += 1
        if hit(d["rp"], pos):
            rp += 1

den = t - rna
systems_on = [s for s, n in sysdet.items() if n >= 2]
out = {"srr": SRR, "dominant": dom, "org": org,
       "target_reads": t, "rrna": rna, "coding": den, "rp": rp, "decoy": dec,
       "pctRP": round(100 * rp / den, 2) if den else None,
       "matrix_systems": dict(sysdet), "systems_on": systems_on,
       "matrix_on": len(systems_on) >= 2,
       "genus_counts": dict(c.most_common(6))}
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "genus_counts"}))
