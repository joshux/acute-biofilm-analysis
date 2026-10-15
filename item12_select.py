#!/usr/bin/env python3
"""item12_quantify.py — two-axis quantification on the matrix-capable no-BP samples
that the OLD mrk-only Klebsiella panel could never score, plus matched 2BP comparators.

Selects, per patient, the run with the most matrix-capable reads at the 40k screen,
then quantifies at 500k reads with the corrected panels.
"""
import os, json, subprocess, collections, sys

os.chdir("/agent/workspace/abf")
os.makedirs("item12", exist_ok=True)

S = json.load(open("screen3_ALL.json"))
J = json.load(open("recon/join.json"))
URL = {}
for p, v in J["by_patient"].items():
    for r in v["runs"]:
        URL[r["run_accession"]] = r.get("fastq_ftp", "")

MATRIX_CAPABLE = {"pseudomonas", "klebsiella", "escherichia", "staphylococcus", "acinetobacter"}

# best run per patient = most matrix-capable reads
best = {}
for srr, r in S.items():
    if "error" in r: continue
    n = sum(v for k, v in r["counts"].items() if k in MATRIX_CAPABLE)
    p = r["patient"]
    if p not in best or n > best[p][0]:
        best[p] = (n, srr, r["status"], r["counts"])

# priority: no-BP with a matrix-capable organism at depth; then 2BP comparators
sel = []
for p, (n, srr, st, cnt) in best.items():
    if n >= 100:
        sel.append((n, p, srr, st, cnt))
sel.sort(key=lambda x: -x[0])
print(f"patients with >=100 matrix-capable reads: {len(sel)}")
for n, p, srr, st, cnt in sel:
    print(f"  {p:7s} {st:6s} {srr:13s} matcap={n:>6} {sorted(cnt.items(), key=lambda x:-x[1])[:3]}")

open("recon/item12_srrs.txt", "w").write("\n".join(x[2] for x in sel))
print(f"\nwrote {len(sel)} runs to recon/item12_srrs.txt")
