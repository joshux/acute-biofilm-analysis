#!/usr/bin/env python3
"""Rebuild reference + quantification panels from committed coordinates.

Inputs (committed to the repo, nothing guessed):
  gff/panels.json   accession -> [{cat, gene, seqid, start, end, ...}]
  refs2/*.fna.gz    panel genomes
  decoy/*.fna.gz    one genome per co-infecting genus (MAPQ competition)

Outputs:
  qp_<org>.fna        >cls|acc|gene|system|i   cls in {matrix,rp,og}
  comp_<org>.fna      qp_<org>.fna + every decoy contig  (single mapping index)
  systems.json        gene -> matrix system map (alginate/psl/pel/other)

The matrix axis is per-system so the quadrant call (>=2 of 3 systems with
>=3 genes each) can be made for any sample.
"""
import json, os, gzip, re, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)

PSEU = ["GCF_000006765.1", "GCF_045571375.1", "GCF_047118485.1"]
KLEB = ["GCA_021989295.1", "GCF_000240185.1", "GCF_988229555.1"]
ORGS = {"pseu": PSEU, "kleb": KLEB}

panels = json.load(open("gff/panels.json"))
decoy_man = json.load(open("decoy/manifest.json"))

RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")


def system_of(gene, cat):
    """map a matrix gene name to its system"""
    g = gene.lower()
    if cat == "matrix" or cat == "matrix_prod":
        if re.match(r"^alg(\d|[a-z])", g) or g in ("mucC", "wspC", "wspA", "wspB",
                                                    "wspD", "wspE", "wspF", "algU"):
            return "alginate"
        if g.startswith("psl"):
            return "psl"
        if g.startswith("pel"):
            return "pel"
        if g.startswith("cdr"):
            return "cdr"
        return "other_matrix"
    return None


def read_genome(acc, sub):
    """yield (contig_id, seq) from a gzipped fna"""
    path = f"{sub}/{acc}.fna.gz"
    with gzip.open(path, "rt") as f:
        h, buf = None, []
        for line in f:
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(buf)
                h = line[1:].split()[0].rstrip()
                buf = []
            else:
                buf.append(line.strip())
        if h is not None:
            yield h, "".join(buf)


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


systems = {}
for org, accs in ORGS.items():
    recs = []
    i = 0
    syscount = collections.Counter()
    clscount = collections.Counter()
    for acc in accs:
        # load genome contigs for this accession
        contigs = dict(read_genome(acc, "refs2"))
        for gene_rec in panels[acc]:
            cat = gene_rec["cat"]
            gene = gene_rec["gene"]
            seqid = gene_rec["seqid"]
            s, e = gene_rec["start"], gene_rec["end"]
            seq = contigs.get(seqid)
            if seq is None:
                # try prefix match
                for k in contigs:
                    if k == seqid or k.startswith(seqid) or seqid.startswith(k):
                        seq = contigs[k]; break
            if seq is None or e > len(seq) or s < 1:
                continue
            sub = seq[s - 1:e]
            if len(sub) < 30:
                continue
            if cat.startswith("matrix"):
                cls = "matrix"
                sysname = system_of(gene, cat)
                systems[f"{acc}|{gene}"] = sysname
                syscount[sysname] += 1
            elif RP_RE.match(gene):
                cls = "rp"; sysname = "-"
            else:
                cls = "og"; sysname = "-"
            clscount[cls] += 1
            recs.append((f"{cls}|{acc}|{gene}|{sysname}|{i}", sub))
            i += 1
    with open(f"qp_{org}.fna", "w") as f:
        for h, s in recs:
            f.write(f">{h}\n{s}\n")
    print(f"{org}: {dict(clscount)}  matrix systems: {dict(syscount)}  total {len(recs)}")
    json.dump({"clscount": dict(clscount), "syscount": dict(syscount),
               "n": len(recs)}, open(f"qp_{org}_meta.json", "w"), indent=1)

json.dump(systems, open("systems.json", "w"), indent=1)
print(f"\nwrote systems.json ({len(systems)} matrix genes labelled)")

# --- concatenate decoys into the competition index -----------------------
decoy_accs = sorted({d["acc"] for d in decoy_man})
n_dec = 0
for org in ORGS:
    with open(f"comp_{org}.fna", "w") as out:
        with open(f"qp_{org}.fna") as f:
            out.write(f.read())
        for acc in decoy_accs:
            p = f"decoy/{acc}.fna.gz"
            if not os.path.exists(p):
                continue
            for cid, seq in read_genome(acc, "decoy"):
                out.write(f">decoy|{acc}|{cid}\n{seq}\n")
                n_dec += 1
    print(f"comp_{org}.fna written (panel + {len(decoy_accs)} decoy genomes)")
print(f"decoy contigs appended: {n_dec}")
