#!/usr/bin/env python3
"""dl_and_quant3.py — download (subsampled) runs using CACHED urls from recon/join.json
and quantify each with quant3.py. Avoids the ENA portal API (currently 500).

usage: dl_and_quant3.py <mode> <SRR> [<SRR> ...]
   mode = 500k | full
"""
import os, sys, json, subprocess, gzip, io

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MODE = sys.argv[1]
SRRS = sys.argv[2:]
NREADS = 500_000 if MODE == "500k" else None

J = json.load(open("recon/join.json"))
URL = {}
for p, v in J["by_patient"].items():
    for r in v["runs"]:
        URL[r["run_accession"]] = r.get("fastq_ftp", "")

os.makedirs("sub", exist_ok=True)
os.makedirs("final3", exist_ok=True)

for srr in SRRS:
    if os.path.exists(f"final3/{srr}.json"):
        print("have", srr); continue
    ftp = URL.get(srr, "")
    if not ftp:
        print("nourl", srr); continue
    urls = ["https://" + u for u in ftp.split(";") if u.strip()]
    if MODE == "500k":
        ok = []
        for i, u in enumerate(urls, 1):
            out = f"sub/{srr}_{i}.fastq.gz"
            if os.path.exists(out) and os.path.getsize(out) > 1000:
                ok.append(out); continue
            cmd = (f"curl -sfL --retry 3 '{u}' | zcat 2>/dev/null | "
                   f"head -n {4*NREADS} | gzip -c > {out}")
            subprocess.run(cmd, shell=True, timeout=3600)
            if os.path.exists(out) and os.path.getsize(out) > 1000:
                ok.append(out)
            else:
                os.path.exists(out) and os.remove(out)
        if not ok:
            print("dlfail", srr); continue
        subprocess.run(f"python3 -u quant3.py {srr} sub final3", shell=True)
    else:
        for i, u in enumerate(urls, 1):
            out = f"fastq/{srr}_{i}.fastq.gz"
            if not (os.path.exists(out) and os.path.getsize(out) > 1000):
                subprocess.run(f"curl -sfL --retry 3 -o '{out}' '{u}'",
                               shell=True, timeout=7200)
        subprocess.run(f"python3 -u quant3.py {srr} fastq final3", shell=True)
