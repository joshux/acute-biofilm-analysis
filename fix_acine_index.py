#!/usr/bin/env python3
"""Rebuild the Acinetobacter quantification index WITHOUT the duplicate
A. baumannii decoy (GCF_009035845.1) — a target-genome decoy collapses MAPQ
for every real Acinetobacter read (the same multi-genome artifact as before,
self-inflicted). Also removes duplicate Pseudomonas/Klebsiella-adjacent decoys
per genus where a target genome exists in the index.
"""
import os, gzip, json

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"

TARGET = "GCF_026167805.1"          # A. baumannii ATCC 17978 (the quantification reference)
DROP = {"GCF_009035845.1"}          # the duplicate A. baumannii decoy


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


with open("rep_acine2.fna", "w") as out:
    for cid, seq in read_gz(f"refs2/{TARGET}.fna.gz"):
        out.write(f">{cid}\n{seq}\n")
    n = 0
    for d in json.load(open("decoy/manifest.json")):
        if d["acc"] in DROP:
            continue
        p = f"decoy/{d['acc']}.fna.gz"
        if not os.path.exists(p):
            continue
        for cid, seq in read_gz(p):
            out.write(f">decoy|{cid}\n{seq}\n")
            n += 1
print(f"target + {n} decoy contigs")
os.system(f"{MM2} -d rep_acine2.mmi rep_acine2.fna 2>/dev/null")
print("built rep_acine2.mmi")
