#!/usr/bin/env python3
"""Build panel_multi3.mmi = all 8 target-organism genomes (+ decoys) for genus assignment."""
import os, gzip, json
os.chdir("/agent/workspace/abf")
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"

TARGETS = {
    "GCF_000006765.1": "pseudomonas", "GCF_000240185.1": "klebsiella",
    "GCF_000005845.2": "escherichia", "GCF_000013425.1": "staphylococcus",
    "GCF_026167805.1": "acinetobacter", "GCF_000012185.1": "haemophilus",
    "GCF_000006885.1": "streptococcus", "GCF_000072485.1": "stenotrophomonas",
}

def read_gz(p):
    with gzip.open(p, "rt", errors="replace") as f:
        h = None; buf = []
        for line in f:
            if line.startswith(">"):
                if h: yield h, "".join(buf)
                h = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
        if h: yield h, "".join(buf)

with open("panel_multi3.fna", "w") as o:
    for acc in TARGETS:
        for cid, seq in read_gz(f"genomes/{acc}.fna.gz"):
            o.write(f">{acc}|{cid}\n{seq}\n")
os.system(f"{MM2} -d panel_multi3.mmi panel_multi3.fna 2>/dev/null")
print("built panel_multi3.mmi", os.path.getsize("panel_multi3.fna")/1e6, "MB")
