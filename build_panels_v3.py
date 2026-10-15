#!/usr/bin/env python3
"""build_panels_v3.py — corrected per-organism two-axis panels.

Fixes the Klebsiella matrix bug: the committed qp_kleb.fna tagged all 62
matrix-ish CDS as a single 'other_matrix' system, so the >=2-systems rule
could never fire and Klebsiella-dominated samples were scored matrix-OFF
by construction.

This version assigns every matrix CDS to a NAMED system, and additionally
tags each system with a TYPE so cross-genus calls can be made commensurable:

  EPS        extracellular polysaccharide matrix (alginate/psl/pel/pga/bcs/cps)
  PROTEIN    proteinaceous fimbriae/adhesin matrix (csg/mrk/fim/csu/cbp/hmw/ecp)
  REG        regulator/signalling only (not counted as a structural system)

System rules are (regex on gene, regex on product). Gene names win when present.

Outputs:
  qp3_<org>.fna          sequences, header: {cls}|{acc}|{gene}|{system}|{i}
  gff_intervals3.json    intervals per (acc|cid): rna / rp / matrix(system,type)
  panels3.json           per-organism census
  systems3.json          (acc|gene) -> system
"""
import os, gzip, json, re, collections, urllib.parse

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")
RP_PROD_RE = re.compile(r"(50S|30S|ribosomal protein)", re.I)

ORGS = {
    "pseu":  ("GCF_000006765.1", "genomes/GCF_000006765.1.fna.gz", "pseudomonas"),
    "kleb":  ("GCF_000240185.1", "genomes/GCF_000240185.1.fna.gz", "klebsiella"),
    "ecoli": ("GCF_000005845.2", "genomes/GCF_000005845.2.fna.gz", "escherichia"),
    "staph": ("GCF_000013425.1", "genomes/GCF_000013425.1.fna.gz", "staphylococcus"),
    "acine": ("GCF_026167805.1", "genomes/GCF_026167805.1.fna.gz", "acinetobacter"),
    "haemo": ("GCF_000012185.1", "genomes/GCF_000012185.1.fna.gz", "haemophilus"),
    "strep": ("GCF_000006885.1", "genomes/GCF_000006885.1.fna.gz", "streptococcus"),
    "steno": ("GCF_000072485.1", "genomes/GCF_000072485.1.fna.gz", "stenotrophomonas"),
}

# (system, type, gene_regex, product_regex)
RULES = {
    "pseu": [
        ("alginate", "EPS",     r"^alg[0-9A-Za-z]+$", r"alginate"),
        ("psl",      "EPS",     r"^psl[A-Z]$",        None),
        ("pel",      "EPS",     r"^pel[A-Z]$",        None),
        ("cdr",      "PROTEIN", r"^cdr[AB]$",         None),
    ],
    "ecoli": [
        ("pga",  "EPS",     r"^pga[ABCD]$",  r"poly-beta-1,6|PGA"),
        ("bcs",  "EPS",     r"^bcs[A-Z]$",   r"cellulose"),
        ("csg",  "PROTEIN", r"^csg[A-G]$",   r"curli"),
        ("fim",  "PROTEIN", r"^fim[A-Z]$",   r"type 1 fimbri"),
        ("ecp",  "PROTEIN", r"^ecp[A-Z]$",   r"common pilus"),
    ],
    "kleb": [
        ("mrk",   "PROTEIN", r"^mrk[A-FJ]$", r"\bMrk[A-F]\b|mannose-resistant|type 3 fimbri"),
        ("fim",   "PROTEIN", r"^fim[A-Z]$",  r"type 1 fimbri"),
        ("ecp",   "PROTEIN", r"^ecp[A-Z]$",  r"common pilus|fimbrillin"),
        ("pga",   "EPS",     r"^pga[ABCD]$", r"poly-beta-1,6|PGA synthesis"),
        ("bcs",   "EPS",     r"^bcs[A-Z]$",  r"cellulose"),
        ("csg",   "PROTEIN", r"^csg[A-G]$",  r"curli"),
        ("capsule", "EPS",   r"^(wzi|wza|wzb|wzc|galF|wcaJ)$",
                             r"capsule assembly|polysaccharide export|Wzi|Wza|Wzb|Wzc"),
        ("cps",   "EPS",     r"^(man[ABC]|gmd|rml[A-D]|ugd|nudK|wca[A-Z]|cps[A-Z])$",
                             r"GDP-mannose|dTDP-rhamnose|UDP-glucose 6-dehydrogenase|capsular polysaccharide"),
    ],
    "staph": [
        ("ica",  "EPS",     r"^ica[ABCDR]$", r"intercellular adhesion|PNAG"),
        ("bap",  "PROTEIN", r"^bap$",        r"biofilm-associated protein"),
        ("fnb",  "PROTEIN", r"^(fnb[AB]|clf[AB])$", r"fibronectin|clumping factor"),
    ],
    "acine": [
        ("pga",  "EPS",     r"^pga[ABCD]$",  r"poly-beta-1,6|PGA"),
        ("csu",  "PROTEIN", r"^csu(AB|[ABCDE])$", r"csu|chaperone-usher"),
        ("bap",  "PROTEIN", r"^bap$",        r"biofilm-associated protein"),
    ],
    "haemo": [
        ("hmw",     "PROTEIN", r"^hmw[12][ABC]$", r"HMW[12]"),
        ("adhesin", "PROTEIN", r"^(hia|hap)$",    r"adhesin|Hia|Hap"),
        ("lic",     "PROTEIN", r"^lic[1-3][A-D]$", None),
    ],
    "strep": [
        ("cbp",  "PROTEIN", r"^cbp[A-Z]$",  r"choline-binding"),
        ("cps",  "EPS",     r"^(cps[0-9]?[A-Z]|wzg|wzh|wzd|wze|wchA)$", r"capsular polysaccharide"),
        ("psa",  "PROTEIN", r"^(psa[ABC]|lyt[ABC]|pce)$", r"pneumococcal surface|choline"),
    ],
    "steno": [
        ("fimbriae", "PROTEIN", r"^smf-?1$", r"fimbrial"),
        ("qs",       "REG",     r"^rpf[FCB]$", None),
    ],
}


def read_gz(p):
    with gzip.open(p, "rt", errors="replace") as f:
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


def classify(org, gene, prod):
    """-> (system, type) or None"""
    for sysname, typ, gre, pre in RULES.get(org, []):
        if gre and gene and re.match(gre, gene, re.I):
            return sysname, typ
    for sysname, typ, gre, pre in RULES.get(org, []):
        if pre and prod and re.search(pre, prod, re.I):
            return sysname, typ
    return None


decoys = json.load(open("decoy/manifest.json"))
GD = {d["acc"]: d["taxon"].split()[0].lower() for d in decoys}
INTERVALS = {}
SYSMAP = {}
summary = {}

for org, (acc, fna, tgenus) in ORGS.items():
    gff = f"genomes/{acc}.gff.gz"
    rows = []
    with gzip.open(gff, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("CDS", "rRNA", "tRNA", "tmRNA"):
                continue
            a = urllib.parse.unquote(c[8])   # GFF attributes are URL-encoded (%2C = comma)
            gm = re.search(r"(?:^|;)gene=([^;]+)", a)
            pm = re.search(r"(?:^|;)product=([^;]+)", a)
            rows.append((c[0], int(c[3]), int(c[4]), c[2],
                         (gm.group(1) if gm else "").strip(),
                         (pm.group(1) if pm else "").strip()))
    iv = collections.defaultdict(lambda: {"rna": [], "rp": [], "matrix": []})
    recs = []; i = 0
    syscnt = collections.Counter(); clscnt = collections.Counter()
    contigs = dict(read_gz(fna))
    for cid, s, e, typ, g, prod in rows:
        key = f"{acc}|{cid}"
        if typ in ("rRNA", "tRNA", "tmRNA"):
            iv[key]["rna"].append((s, e)); continue
        hit = classify(org, g, prod)
        if hit:
            sysname, styp = hit
            iv[key]["matrix"].append((s, e, sysname, styp))
            syscnt[f"{sysname}({styp})"] += 1
            cls = "matrix"
            SYSMAP[f"{acc}|{g or prod[:30]}"] = sysname
        elif RP_RE.match(g) or RP_PROD_RE.search(prod):
            iv[key]["rp"].append((s, e)); cls = "rp"
        else:
            cls = "og"
        clscnt[cls] += 1
        seq = contigs.get(cid)
        if seq is None:
            for k in contigs:
                if k == cid or k.startswith(cid) or cid.startswith(k):
                    seq = contigs[k]; break
        if seq is None or e > len(seq) or s < 1:
            continue
        sub = seq[s-1:e]
        if len(sub) < 30:
            continue
        recs.append((f"{cls}|{acc}|{g or 'na'}|{hit[0] if hit else '-'}|{i}", sub))
        i += 1
    INTERVALS[org] = {k: v for k, v in iv.items()}
    with open(f"qp3_{org}.fna", "w") as f:
        for h, s in recs:
            f.write(f">{h}\n{s}\n")
    with open(f"comp3_{org}.fna", "w") as out:
        out.write(open(f"qp3_{org}.fna").read())
        for d in decoys:
            if GD.get(d["acc"]) == tgenus:
                continue
            p = f"decoy/{d['acc']}.fna.gz"
            if not os.path.exists(p):
                continue
            for cid, seq in read_gz(p):
                out.write(f">decoy|{cid}\n{seq}\n")
    nm = sum(len(v["matrix"]) for v in iv.values())
    summary[org] = {"acc": acc, "classes": dict(clscnt),
                    "systems": dict(syscnt),
                    "n_systems": len({s.split("(")[0] for s in syscnt}),
                    "n_eps_systems": len({s.split("(")[0] for s in syscnt if "(EPS)" in s}),
                    "matrix_genes": nm}
    print(f"{org:6s} {dict(clscnt)}")
    print(f"        systems={dict(syscnt)}")

json.dump(summary, open("panels3.json", "w"), indent=1)
json.dump(INTERVALS, open("gff_intervals3.json", "w"))
json.dump(SYSMAP, open("systems3.json", "w"), indent=1)
print("\nwrote panels3.json, gff_intervals3.json, systems3.json")
for org in ORGS:
    os.system(f"{MM2} -d comp3_{org}.mmi comp3_{org}.fna 2>/dev/null")
print("built all comp3_<org>.mmi")
