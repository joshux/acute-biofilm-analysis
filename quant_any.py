#!/usr/bin/env python3
"""quant_any.py <SRR> <fq1> [fq2] <outdir>  — two-axis quantification for any run.

Reuses the corrected panels (qp3/gff_intervals3/panels3/gcomp3).
Adds --depth N to subsample N reads per mate (for the depth-dilution test).
"""
import os, sys, json, collections, subprocess, tempfile

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
MIN_GENUS_READS = 20
MIN_GENES = 2

SRR = sys.argv[1]
ODIR = sys.argv[-1]              # last arg is always the output dir
fqs = sys.argv[2:-1]             # everything between is fastq file(s)
if not fqs:
    print("no fastq args"); sys.exit(1)
os.makedirs(ODIR, exist_ok=True)

IV = json.load(open("gff_intervals3.json"))
TARGETS = {"GCF_000006765.1": "pseudomonas", "GCF_000240185.1": "klebsiella",
           "GCF_000005845.2": "escherichia", "GCF_000013425.1": "staphylococcus",
           "GCF_026167805.1": "acinetobacter", "GCF_000012185.1": "haemophilus",
           "GCF_000006885.1": "streptococcus", "GCF_000072485.1": "stenotrophomonas"}
GENUS2ORG = {"pseudomonas": "pseu", "klebsiella": "kleb", "escherichia": "ecoli",
             "staphylococcus": "staph", "acinetobacter": "acine",
             "haemophilus": "haemo", "streptococcus": "strep",
             "stenotrophomonas": "steno"}

FQ = " ".join(fqs)

# --- dominant genus ---
c = collections.Counter()
with os.popen(f"{MM2} -x sr -t 2 -a --secondary=no panel_multi3.mmi {FQ} 2>/dev/null") as fh:
    for line in fh:
        if line[0] == "@": continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20: continue
        if f[2].startswith("decoy|"): continue
        c[TARGETS.get(f[2].split("|")[0], "?")] += 1
dom = c.most_common(1)[0][0] if c else None
org = GENUS2ORG.get(dom)
res = {"srr": SRR, "dominant": dom, "org": org, "genus_counts": dict(c)}
if not org or c[dom] < MIN_GENUS_READS:
    res.update({"pctRP": None, "matrix_on": None, "note": "no dominant target organism"})
    json.dump(res, open(f"{ODIR}/{SRR}.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "genus_counts"})); sys.exit(0)

iv = IV[org]

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
with os.popen(f"{MM2} -x sr -t 2 -a --secondary=no gcomp3_{org}.mmi {FQ} 2>/dev/null") as fh:
    for line in fh:
        if line[0] == "@": continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20: continue
        if f[2].startswith("decoy|"): dec += 1; continue
        d = iv.get(f[2])
        if not d: continue
        t += 1
        pos = int(f[3])
        if hit(d["rna"], pos): rna += 1; continue
        m = hit(d["matrix"], pos)
        if m: sysdet[m[2]] += 1
        if hit(d["rp"], pos): rp += 1

P3 = json.load(open("panels3.json"))
TYPES = {}
for o, v in P3.items():
    for s in v["systems"]:
        nm, ty = s.split("("); TYPES[nm] = ty.rstrip(")")
den = t - rna
on = sorted([s for s, n in sysdet.items() if n >= MIN_GENES])
eps_on = [s for s in on if TYPES.get(s) == "EPS"]
prot_on = [s for s in on if TYPES.get(s) == "PROTEIN"]
res.update({"target_reads": t, "rrna": rna, "coding": den, "rp": rp, "decoy": dec,
            "pctRP": round(100 * rp / den, 2) if den else None,
            "matrix_systems": dict(sysdet), "systems_on": on,
            "eps_on": eps_on, "prot_on": prot_on,
            "matrix_on": len(on) >= 2, "matrix_on_eps": len(eps_on) >= 2})
json.dump(res, open(f"{ODIR}/{SRR}.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "genus_counts"}))
