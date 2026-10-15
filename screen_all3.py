#!/usr/bin/env python3
"""screen_all3.py — 40k-read pan-genus screen of ALL PRJNA1033689 runs.

Rebuilds the composition table with the corrected 8-genus panel so the
committed claims (e.g. "HS243 Esch 7.2k") can be independently re-tested.

Reads URLs from recon/join.json (ENA portal API is 500 right now).
usage: screen_all3.py [nreads]
"""
import os, sys, json, subprocess, collections

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
NREADS = int(sys.argv[1]) if len(sys.argv) > 1 else 40000

TARGETS = {"GCF_000006765.1": "pseudomonas", "GCF_000240185.1": "klebsiella",
           "GCF_000005845.2": "escherichia", "GCF_000013425.1": "staphylococcus",
           "GCF_026167805.1": "acinetobacter", "GCF_000012185.1": "haemophilus",
           "GCF_000006885.1": "streptococcus", "GCF_000072485.1": "stenotrophomonas"}

J = json.load(open("recon/join.json"))
jobs = []
for p, v in J["by_patient"].items():
    for r in v["runs"]:
        jobs.append((r["run_accession"], p, v["status"], r.get("fastq_ftp", ""),
                     int(r.get("read_count") or 0)))

os.makedirs("screen3", exist_ok=True)
os.makedirs("screenfq", exist_ok=True)
out = {}
for i, (srr, pat, st, ftp, rc) in enumerate(sorted(jobs, key=lambda x: -x[4]), 1):
    f = f"screen3/{srr}.json"
    if os.path.exists(f):
        out[srr] = json.load(open(f)); continue
    urls = ["https://" + u for u in ftp.split(";") if u.strip()]
    if not urls:
        out[srr] = {"error": "nourl"}; continue
    # take first NREADS from mate1 only (matches the committed screen design)
    fq = f"screenfq/{srr}.fq"
    if not (os.path.exists(fq) and os.path.getsize(fq) > 0):
        cmd = (f"curl -sfL --retry 3 '{urls[0]}' | zcat 2>/dev/null | "
               f"head -n {4*NREADS} > {fq}")
        subprocess.run(cmd, shell=True, timeout=1800)
    if not (os.path.exists(fq) and os.path.getsize(fq) > 0):
        out[srr] = {"error": "dlfail"}; continue
    c = collections.Counter()
    cmd = f"{MM2} -x sr -t 2 -a --secondary=no panel_multi3.mmi {fq} 2>/dev/null"
    with os.popen(cmd) as fh:
        for line in fh:
            if line[0] == "@": continue
            f2 = line.split("\t", 6)
            if len(f2) < 7 or int(f2[1]) & 4 or int(f2[4]) < 20: continue
            c[TARGETS.get(f2[2].split("|")[0], "?")] += 1
    res = {"srr": srr, "patient": pat, "status": st, "read_count": rc,
           "counts": dict(c), "total": sum(c.values())}
    json.dump(res, open(f, "w"), indent=1)
    out[srr] = res
    os.remove(fq)
    print(f"[{i}/{len(jobs)}] {srr} {pat} {st} total={res['total']} "
          f"top={c.most_common(3)}", flush=True)

json.dump(out, open("screen3_ALL.json", "w"), indent=1)
print("SCREEN DONE", flush=True)
