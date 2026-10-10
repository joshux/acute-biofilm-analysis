#!/usr/bin/env python3
"""Build Gifford-faithful per-organism indexes: FULL genome + other-genus decoys.

  gcomp_<org>.mmi  = full target genome + decoys of other genera
  gff_intervals.json = per organism: {(acc,cid): {'rna':[(s,e)],'rp':[(s,e)],'matrix':[(s,e,sys)]}}
"""
import os, gzip, json, re, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")
RP_PROD_RE = re.compile(r"(50S|30S|ribosomal protein)", re.I)

ORGS = {
    "pseu":  ("GCF_000006765.1", "PANELS",        "refs2/GCF_000006765.1.fna.gz", "pseudomonas"),
    "kleb":  ("GCF_000240185.1", "gff_kleb.gz",   "refs2/GCF_000240185.1.fna.gz", "klebsiella"),
    "ecoli": ("GCF_000005845.2", "gff_ecoli.gz",  "refs2/GCF_000005845.2.fna.gz", "escherichia"),
    "staph": ("GCF_000013425.1", "gff_strep.gz",  "refs2/GCF_000013425.1.fna.gz", "staphylococcus"),
    "acine": ("GCF_026167805.1", "gff_acine.gz",  "refs2/GCF_026167805.1.fna.gz", "acinetobacter"),
    "haemo": ("GCF_000012185.1", "gff_staph.gz",  "refs2/GCF_000012185.1.fna.gz", "haemophilus"),
    "strep": ("GCF_000006885.1", "gff_haemo.gz",  "refs2/GCF_000006885.1.fna.gz", "streptococcus"),
    "steno": ("GCF_000072485.1", "gff_steno.gz",  "refs2/GCF_000072485.1.fna.gz", "stenotrophomonas"),
}
RULES = {
    "pseu": [("alginate", r"^(alg[0-9A-Za-z]+|mucC|wspC)$", None),
             ("psl", r"^psl[A-Z]$", None), ("pel", r"^pel[A-Z]$", None),
             ("cdr", r"^cdr[AB]$", None)],
    "ecoli": [("pga", r"^pga[ABCD]$", None), ("csg", r"^csg[A-G]$", None), ("bcs", r"^bcs[A-Z]$", None)],
    "acine": [("pga", r"^pga[ABCD]$", None), ("csu", r"^csu(AB|[ABCDE])$", None)],
    "staph": [("ica", r"^ica[ABCDR]$", r"intercellular adhesion"), ("bap", r"^bap$", None)],
    "kleb":  [("mrk", r"^mrk[A-FJ]$", r"^Mrk[A-F] "), ("fim", r"^fim[A-Z]$", None),
              ("cps", r"^(wzi|wza|wzb|wzc|galF)$", None)],
    "haemo": [("hmw", r"^hmw[12][ABC]$", r"HMW[12]"), ("adhesin", r"^(hia|hap)$", None),
              ("lic", r"^lic[1-3][A-D]$", None)],
    "strep": [("cps", r"^(cps[0-9]?[A-Z]|wzg|wzh|wzd|wze|wchA)$", r"capsular polysaccharide"),
              ("cbp", r"^cbp[A-Z]$", None)],
    "steno": [("fimbriae", r"^smf-?1$", r"fimbrial"), ("qs", r"^rpf[FCB]$", None)],
}


def classify(org, gene, prod):
    for s, gre, pre in RULES.get(org, []):
        if gre and gene and re.match(gre, gene, re.I):
            return s
        if pre and prod and re.search(pre, prod, re.I):
            return s
    return None


def read_gz(p):
    with gzip.open(p, "rt") as f:
        h = None; buf = []
        for line in f:
            if line.startswith(">"):
                if h:
                    yield h, "".join(buf)
                h = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if h:
            yield h, "".join(buf)


decoys = json.load(open("decoy/manifest.json"))
GD = {d["acc"]: d["taxon"].split()[0].lower() for d in decoys}
panels = json.load(open("gff/panels.json"))
INTERVALS = {}

for org, (acc, gff, fna, tgenus) in ORGS.items():
    if gff == "PANELS":
        rows = [(g["seqid"], int(g["start"]), int(g["end"]), g["gene"], "") for g in panels[acc]]
    else:
        rows = []
        with gzip.open(gff, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 9 or c[2] not in ("CDS", "rRNA", "tRNA", "tmRNA"):
                    continue
                a = c[8]
                gm = re.search(r"(?:^|;)gene=([^;]+)", a)
                pm = re.search(r"(?:^|;)product=([^;]+)", a)
                rows.append((c[0], int(c[3]), int(c[4]), c[2],
                             gm.group(1) if gm else "", pm.group(1) if pm else ""))
    iv = collections.defaultdict(lambda: {"rna": [], "rp": [], "matrix": []})
    for r in rows:
        if gff == "PANELS":
            cid, s, e, g, prod = r
            typ = "CDS"
        else:
            cid, s, e, typ, g, prod = r
        key = f"{acc}|{cid}"
        if typ in ("rRNA", "tRNA", "tmRNA"):
            iv[key]["rna"].append((s, e))
        elif typ == "CDS":
            sysname = classify(org, g, prod)
            if sysname:
                iv[key]["matrix"].append((s, e, sysname))
            elif RP_RE.match(g) or RP_PROD_RE.search(prod):
                iv[key]["rp"].append((s, e))
    INTERVALS[org] = {k: v for k, v in iv.items()}
    # genome index
    with open(f"gcomp_{org}.fna", "w") as out:
        for cid, seq in read_gz(fna):
            out.write(f"{acc}|{cid}\n{seq}\n".replace(f"{acc}|{cid}\n", f">{acc}|{cid}\n", 1)
                      if False else f">{acc}|{cid}\n{seq}\n")
        for d in decoys:
            if GD.get(d["acc"]) == tgenus:
                continue
            p = f"decoy/{d['acc']}.fna.gz"
            if not os.path.exists(p):
                continue
            for cid, seq in read_gz(p):
                out.write(f">decoy|{cid}\n{seq}\n")
    os.system(f"{MM2} -d gcomp_{org}.mmi gcomp_{org}.fna 2>/dev/null")
    nm = sum(len(v["matrix"]) for v in INTERVALS[org].values())
    nrp = sum(len(v["rp"]) for v in INTERVALS[org].values())
    nrna = sum(len(v["rna"]) for v in INTERVALS[org].values())
    print(f"{org:6s} intervals: rna={nrna:3d} rp={nrp:3d} matrix={nm:3d}")

json.dump(INTERVALS, open("gff_intervals.json", "w"))
print("\nwrote gff_intervals.json; built all gcomp_<org>.mmi")
