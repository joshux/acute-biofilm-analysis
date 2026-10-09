#!/usr/bin/env python3
"""MAPQ diagnostic: are ribosomal-protein reads being dropped by the MAPQ filter?

Maps once at MAPQ>=0 and bins every aligned read by (gene class, MAPQ). If RP
reads pile up at MAPQ 0 while matrix/og reads do not, the %RP metric is being
deflated by multi-genome mapping competition.

Run for both the multi-genome panel (comp_pseu) and the single-genome panel
(comp1_pseu) so the difference is visible.
"""
import sys, os, json, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"


def load_refs(fn):
    refs = {}
    with open(fn) as f:
        for line in f:
            if line.startswith(">"):
                p = line[1:].strip().split("|")
                refs[line[1:].strip()] = p     # [cls, acc, gene, system, i]
    return refs


def run(panel_mmi, panel_fna, fastq):
    refs = load_refs(panel_fna)
    cls_of = {h: v[0] for h, v in refs.items()}
    bycls = collections.defaultdict(collections.Counter)   # cls -> {mapq: n}
    dec = collections.Counter()
    cmd = f"{MM2} -x sr -t 2 -a --secondary=no {panel_mmi} {fastq} 2>/dev/null"
    with os.popen(cmd) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.split("\t", 6)
            if len(f) < 7 or int(f[1]) & 4:
                continue
            mq = int(f[4])
            ref = f[2]
            if ref in refs:
                bycls[cls_of[ref]][mq] += 1
            else:
                dec[mq] += 1
    res = {}
    for thr in (0, 5, 10, 20, 30):
        rp = sum(c for m, c in bycls["rp"].items() if m >= thr)
        og = sum(c for m, c in bycls["og"].items() if m >= thr)
        mt = sum(c for m, c in bycls["matrix"].items() if m >= thr)
        tot = rp + og + mt
        res[thr] = {"rp": rp, "og": og, "matrix": mt,
                    "pctRP": round(100 * rp / tot, 2) if tot else None}
    return res, bycls


SAMPLES = [
    ("SRR27343249", "acute"),
    ("SRR27343409", "acute"),
    ("SRR27343405", "acute"),
    ("SRR27343410", "acute"),
    ("ERR2275089", "chronic-invivo"),
    ("ERR2275090", "chronic-invivo"),
    ("ERR2275092", "chronic-invivo"),
    ("ERR2591781", "chronic-exp"),
    ("ERR2591782", "chronic-exp"),
]

out = {}
for srr, grp in SAMPLES:
    fastq = f"fastq/{srr}.fastq.gz"
    if not os.path.exists(fastq):
        print(f"MISSING {srr}"); continue
    multi, bym = run("comp_pseu.mmi", "qp_pseu.fna", fastq)
    single, bys = run("comp1_pseu.mmi", "qp1_pseu.fna", fastq)
    # median MAPQ per class (multi panel)
    def med(cls):
        vals = sorted(m for m, c in bym[cls].items() for _ in range(c))
        return vals[len(vals) // 2] if vals else None
    out[srr] = {"group": grp, "multi": multi, "single": single,
                "medianMAPQ": {c: med(c) for c in ("rp", "og", "matrix")}}
    print(f"\n{srr}  [{grp}]")
    print(f"  median MAPQ   rp={med('rp')}  og={med('og')}  matrix={med('matrix')}")
    print(f"  {'thr':>4s} | {'MULTI-genome':^28s} | {'SINGLE-genome (PAO1)':^28s}")
    print(f"  {'':>4s} | {'rp':>5s}{'og':>6s}{'mat':>6s}{'%RP':>7s} | {'rp':>5s}{'og':>6s}{'mat':>6s}{'%RP':>7s}")
    for thr in (0, 5, 10, 20, 30):
        a = multi[thr]; b = single[thr]
        print(f"  {thr:>4d} | {a['rp']:>5d}{a['og']:>6d}{a['matrix']:>6d}{a['pctRP']:>7.2f} | "
              f"{b['rp']:>5d}{b['og']:>6d}{b['matrix']:>6d}{b['pctRP']:>7.2f}")

json.dump(out, open("diag_mapq.json", "w"), indent=1)
print("\nwrote diag_mapq.json")
