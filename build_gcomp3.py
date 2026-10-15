#!/usr/bin/env python3
"""Build gcomp3_<org>.mmi = FULL target genome + other-genus decoys.

The two-axis metric maps to the full genome and classifies each read by
interval (rna / rp / matrix). Intervals come from gff_intervals3.json.
"""
import os, gzip, json

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"

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


decoys = json.load(open("decoy/manifest.json"))
GD = {d["acc"]: d["taxon"].split()[0].lower() for d in decoys}

for org, (acc, fna, tgenus) in ORGS.items():
    out = f"gcomp3_{org}.fna"
    with open(out, "w") as o:
        for cid, seq in read_gz(fna):
            o.write(f">{acc}|{cid}\n{seq}\n")
        for d in decoys:
            if GD.get(d["acc"]) == tgenus:
                continue
            p = f"decoy/{d['acc']}.fna.gz"
            if not os.path.exists(p):
                continue
            for cid, seq in read_gz(p):
                o.write(f">decoy|{cid}\n{seq}\n")
    os.system(f"{MM2} -d gcomp3_{org}.mmi {out} 2>/dev/null")
    print(f"{org:6s} built gcomp3_{org}.mmi  ({os.path.getsize(out)/1e6:.1f} MB fasta)")
print("done")
