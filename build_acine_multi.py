#!/usr/bin/env python3
"""Multi-strain Acinetobacter panel, per-organism MAPQ-safe.

Genus assignment index  : all 6 Acinetobacter genomes + decoys (multi-strain,
                           validated earlier: raises assignment 6.3x)
Quantification index     : ONE strain (ATCC 17978) + non-Acinetobacter decoys
                           (single target genome = conserved RP reads keep MAPQ)
Matrix/RP intervals      : union of pga/csu annotations from all 6 strains,
                           lifted onto the 17978 coordinates where the locus
                           is present in 17978 itself (quantification ref).
"""
import os, gzip, json, re, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"

ACC = "GCF_026167805.1"     # 17978 — quantification reference
OTHERS = ["GCF_000241685.1", "GCF_057224575.1", "GCF_055475515.1",
          "GCF_988236795.1", "GCF_050703235.1"]
DROP_DECOY = {"GCF_009035845.1"}
MATRIX = {"pgaA", "pgaB", "pgaC", "pgaD", "csuA", "csuAB", "csuB", "csuC", "csuD", "csuE"}


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


# --- 1. genus-assignment index: 6 Acinetobacter + decoys ------------------
with open("acine_multi.fna", "w") as out:
    for a in [ACC] + OTHERS:
        p = f"refs2/{a}.fna.gz"
        for cid, seq in read_gz(p):
            out.write(f">{a}|{cid}\n{seq}\n")
    for d in json.load(open("decoy/manifest.json")):
        if d["acc"] in DROP_DECOY or not os.path.exists(f"decoy/{d['acc']}.fna.gz"):
            continue
        for cid, seq in read_gz(f"decoy/{d['acc']}.fna.gz"):
            out.write(f">decoy|{cid}\n{seq}\n")
os.system(f"{MM2} -d acine_multi.mmi acine_multi.fna 2>/dev/null")
print("built acine_multi.mmi (genus assignment)")

# --- 2. quantify how much of each sample's reads map to 17978 vs others ---
# (the pilot showed single-strain 17978 yields poorly for clinical strains;
#  the multi index tells us WHICH strain dominates when 17978 does not)
json.dump({"quant_ref": ACC, "others": OTHERS,
           "matrix_genes": sorted(MATRIX)},
          open("acine_panel.json", "w"), indent=1)
print("wrote acine_panel.json")
