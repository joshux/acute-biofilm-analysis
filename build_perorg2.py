#!/usr/bin/env python3
"""Per-organism two-axis panels — corrected.

Two bugs fixed from the first attempt:
  1. genome<->GFF pairing was ROTATED. Corrected by contig identity:
       GCF_000013425.1 = S. aureus   (gff_strep.gz, NC_007795.1)
       GCF_000012185.1 = H. influenzae (gff_staph.gz, NC_007146.2)
       GCF_000006885.1 = S. pneumoniae (gff_haemo.gz, NC_003028.3)
  2. several GFFs have no gene= attribute; matrix genes must be matched by
     product= text as well (e.g. S. aureus icaA-D = "intercellular adhesion
     protein B/C").

Matrix gene matching is per organism via (gene_regex, product_regex) pairs.
"""
import os, gzip, json, re, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")

# organism -> (genome accession, gff, genome fna, target genus)
ORGS = {
    "pseu":  ("GCF_000006765.1", "PANELS",         "refs2/GCF_000006765.1.fna.gz", "pseudomonas"),
    "kleb":  ("GCF_000240185.1", "gff_kleb.gz",    "refs2/GCF_000240185.1.fna.gz", "klebsiella"),
    "ecoli": ("GCF_000005845.2", "gff_ecoli.gz",   "refs2/GCF_000005845.2.fna.gz", "escherichia"),
    "staph": ("GCF_000013425.1", "gff_strep.gz",   "refs2/GCF_000013425.1.fna.gz", "staphylococcus"),
    "acine": ("GCF_026167805.1", "gff_acine.gz",   "refs2/GCF_026167805.1.fna.gz", "acinetobacter"),
    "haemo": ("GCF_000012185.1", "gff_staph.gz",   "refs2/GCF_000012185.1.fna.gz", "haemophilus"),
    "strep": ("GCF_000006885.1", "gff_haemo.gz",   "refs2/GCF_000006885.1.fna.gz", "streptococcus"),
    "steno": ("GCF_000072485.1", "gff_steno.gz",   "refs2/GCF_000072485.1.fna.gz", "stenotrophomonas"),
}

# organism -> list of (system_name, gene_regex, product_regex)
RULES = {
    "pseu": [("alginate", r"^(alg[0-9A-Za-z]+|mucC|wspC)$", None),
             ("psl", r"^psl[A-Z]$", None),
             ("pel", r"^pel[A-Z]$", None),
             ("cdr", r"^cdr[AB]$", None)],
    "ecoli": [("pga", r"^pga[ABCD]$", None),
              ("csg", r"^csg[A-G]$", None),
              ("bcs", r"^bcs[A-Z]$", None)],
    "acine": [("pga", r"^pga[ABCD]$", None),
              ("csu", r"^csu(AB|[ABCDE])$", None)],
    "staph": [("ica", r"^ica[ABCDR]$", r"intercellular adhesion"),
              ("bap", r"^bap$", None)],
    "kleb":  [("mrk", r"^mrk[A-FJ]$", r"^Mrk[A-F] "),
              ("fim", r"^fim[A-Z]$", None),
              ("cps", r"^(wzi|wza|wzb|wzc|galF)$", None)],
    "haemo": [("hmw", r"^hmw[12][ABC]$", r"HMW[12]"),
              ("adhesin", r"^(hia|hap)$", None),
              ("lic", r"^lic[1-3][A-D]$", None)],
    "strep": [("cps", r"^(cps[0-9]?[A-Z]|wzg|wzh|wzd|wze|wchA)$", r"capsular polysaccharide"),
              ("cbp", r"^cbp[A-Z]$", None)],
    "steno": [("fimbriae", r"^smf-?1$", r"fimbrial"),
              ("qs", r"^rpf[FCB]$", None)],
}


def classify(org, gene, product):
    for sysname, gre, pre in RULES.get(org, []):
        if gre and gene and re.match(gre, gene, re.I):
            return sysname
        if pre and product and re.search(pre, product, re.I):
            return sysname
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
GENUS_OF_DECOY = {d["acc"]: d["taxon"].split()[0].lower() for d in decoys}
panels = json.load(open("gff/panels.json"))

summary = {}
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
                if len(c) < 9 or c[2] != "CDS":
                    continue
                a = c[8]
                gm = re.search(r"(?:^|;)gene=([^;]+)", a)
                pm = re.search(r"(?:^|;)product=([^;]+)", a)
                rows.append((c[0], int(c[3]), int(c[4]),
                             gm.group(1) if gm else "",
                             pm.group(1) if pm else ""))
    contigs = dict(read_gz(fna))
    recs = []; i = 0
    syscnt = collections.Counter(); clscnt = collections.Counter(); missing = 0
    for cid, s, e, g, prod in rows:
        seq = contigs.get(cid)
        if seq is None:
            for k in contigs:
                if k == cid or k.startswith(cid) or cid.startswith(k):
                    seq = contigs[k]; break
        if seq is None or e > len(seq) or s < 1:
            missing += 1; continue
        sub = seq[s-1:e]
        if len(sub) < 30:
            continue
        sysname = classify(org, g, prod)
        if sysname:
            cls = "matrix"; syscnt[sysname] += 1
        elif RP_RE.match(g):
            cls = "rp"
        else:
            cls = "og"
        clscnt[cls] += 1
        recs.append((f"{cls}|{acc}|{g or 'na'}|{sysname or '-'}|{i}", sub))
        i += 1
    with open(f"qp_{org}.fna", "w") as f:
        for h, s in recs:
            f.write(f">{h}\n{s}\n")
    with open(f"comp_{org}.fna", "w") as out:
        out.write(open(f"qp_{org}.fna").read())
        for d in decoys:
            if GENUS_OF_DECOY.get(d["acc"]) == tgenus:
                continue
            p = f"decoy/{d['acc']}.fna.gz"
            if not os.path.exists(p):
                continue
            for cid, seq in read_gz(p):
                out.write(f">decoy|{cid}\n{seq}\n")
    summary[org] = {"acc": acc, "classes": dict(clscnt), "systems": dict(syscnt)}
    print(f"{org:6s} {dict(clscnt)}  systems={dict(syscnt)}  (skipped {missing})")

json.dump(summary, open("perorg_panels.json", "w"), indent=1)
for org in ORGS:
    os.system(f"{MM2} -d comp_{org}.mmi comp_{org}.fna 2>/dev/null")
print("\nbuilt all comp_<org>.mmi")
