#!/usr/bin/env python3
"""item4_positive_control.py — run the two-axis assay on known-positive datasets.

Cornforth 2018 PRJNA437845:
  - in vitro PA14 coculture biofilm  (expected: matrix-ON)
  - in vivo human CF sputum, Pseudomonas-dominant (expected: matrix-ON, SLOW growth)

Writes results to item4/.
"""
import os, json, subprocess, sys

os.chdir("/agent/workspace/abf")
os.makedirs("item4", exist_ok=True)

# (label, srr) chosen from the ENA listing: in-vitro biofilm vs in-vivo sputum
TARGETS = [
    ("invitro_biofilm_SRR6833312", "SRR6833312"),
    ("invitro_biofilm_SRR6833317", "SRR6833317"),
    ("invivo_sputum_SRR6833345",   "SRR6833345"),
    ("invivo_sputum_SRR6833347",   "SRR6833347"),
    ("invivo_sputum_SRR6833349",   "SRR6833349"),
]

rows = [l.split("\t") for l in open("recon/cornforth.tsv").read().rstrip("\n").split("\n")]
h = rows[0]
URL = {dict(zip(h, r))["run_accession"]: dict(zip(h, r)).get("fastq_ftp", "") for r in rows[1:]}

for label, srr in TARGETS:
    if os.path.exists(f"item4/{srr}.json"):
        print("have", srr, flush=True); continue
    ftp = URL.get(srr, "")
    urls = ["https://" + u for u in ftp.split(";") if u.strip()]
    if not urls:
        print("nourl", srr, flush=True); continue
    fqs = []
    for i, u in enumerate(urls, 1):
        out = f"item4/{srr}_{i}.fastq.gz"
        if not (os.path.exists(out) and os.path.getsize(out) > 1000):
            # 500k reads per mate
            subprocess.run(f"curl -sfL --retry 3 '{u}' | zcat 2>/dev/null | head -n 2000000 | gzip -c > {out}",
                           shell=True, timeout=5400)
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            fqs.append(out)
    if not fqs:
        print("dlfail", srr, flush=True); continue
    subprocess.run(f"python3 -u quant_any.py {srr} {' '.join(fqs)} item4", shell=True)
    if os.path.exists(f"item4/{srr}.json"):
        d = json.load(open(f"item4/{srr}.json"))
        print(f"{label}: dom={d.get('dominant')} pctRP={d.get('pctRP')} "
              f"target={d.get('target_reads')} sys={d.get('systems_on')}", flush=True)
print("ITEM4 DONE", flush=True)
