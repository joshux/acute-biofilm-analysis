#!/usr/bin/env python3
"""item4b_rossi.py — Rossi 2018 PRJEB24688 CF sputum (RiboZero) as a deep, independent
positive control for the matrix axis. Also tests the depth-dilution behaviour (Item 5b):
quantify the same deep run at several read depths and watch %RP stability.
"""
import os, json, subprocess

os.chdir("/agent/workspace/abf")
os.makedirs("item4b", exist_ok=True)
os.makedirs("item5", exist_ok=True)

rows = [l.split("\t") for l in open("recon/rossi.tsv").read().rstrip("\n").split("\n")]
h = rows[0]
recs = [dict(zip(h, r)) for r in rows[1:]]
# in-vivo sputum runs are the INVIVO / sputum-titled ones; skip EXP/3H/STAT/24H lab cultures
invivo = [r for r in recs if "INVIVO" in (r.get("sample_title", "") or "").upper()
          or "SPUTUM" in (r.get("sample_title", "") or "").upper()]
print("invivo candidates:", len(invivo))
for r in invivo[:6]:
    print(" ", r["run_accession"], f"{int(r['read_count']):>12,}", r.get("sample_title", "")[:46])

# take the two deepest invivo
invivo.sort(key=lambda r: -int(r["read_count"] or 0))
picks = invivo[:2] if invivo else recs[:2]
print("picks:", [p["run_accession"] for p in picks])

for r in picks:
    srr = r["run_accession"]
    if os.path.exists(f"item4b/{srr}.json"):
        print("have", srr, flush=True); continue
    urls = ["https://" + u for u in (r.get("fastq_ftp", "") or "").split(";") if u.strip()]
    fqs = []
    for i, u in enumerate(urls, 1):
        out = f"item4b/{srr}_{i}.fastq.gz"
        if not (os.path.exists(out) and os.path.getsize(out) > 1000):
            subprocess.run(f"curl -sfL --retry 3 '{u}' | zcat 2>/dev/null | head -n 2000000 | gzip -c > {out}",
                           shell=True, timeout=5400)
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            fqs.append(out)
    if not fqs:
        print("dlfail", srr, flush=True); continue
    subprocess.run(f"python3 -u quant_any.py {srr} {' '.join(fqs)} item4b", shell=True)
print("ITEM4B DONE", flush=True)
