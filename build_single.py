#!/usr/bin/env python3
"""Build a SINGLE-genome (PAO1) quantification panel + decoys.

Purpose: test whether the multi-genome panel's MAPQ competition is dropping
conserved ribosomal-protein reads (which would deflate %RP).
"""
import json, os, gzip

ABF = "/agent/workspace/abf"
os.chdir(ABF)
ACC = "GCF_000006765.1"

keep = []
with open("qp_pseu.fna") as f:
    cur = None
    for line in f:
        if line.startswith(">"):
            h = line[1:].strip()
            cur = [h, []] if h.split("|")[1] == ACC else None
            if cur:
                keep.append(cur)
        elif cur is not None:
            cur[1].append(line.strip())

with open("qp1_pseu.fna", "w") as f:
    for h, seq in keep:
        f.write(">" + h + "\n" + "".join(seq) + "\n")
print("single-genome (PAO1) panel genes:", len(keep))


def read_genome(path):
    with gzip.open(path, "rt") as f:
        h = None
        buf = []
        for line in f:
            if line.startswith(">"):
                if h:
                    yield h, "".join(buf)
                h = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if h:
            yield h, "".join(buf)


decoy = [d["acc"] for d in json.load(open("decoy/manifest.json"))]
with open("comp1_pseu.fna", "w") as out:
    out.write(open("qp1_pseu.fna").read())
    n = 0
    for acc in decoy:
        p = f"decoy/{acc}.fna.gz"
        if not os.path.exists(p):
            continue
        for cid, seq in read_genome(p):
            out.write(f">decoy|{acc}|{cid}\n{seq}\n")
            n += 1
print("decoy contigs:", n)
