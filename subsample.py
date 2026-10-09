#!/usr/bin/env python3
"""Stream-subsample FASTQ from ENA: download only the first N reads of each run.

Never stores the full file. For paired runs, takes N/2 reads from each mate and
concatenates (mapping is per-read, so mate order is irrelevant).

usage: subsample.py <runs.txt> <reads_cap> <outdir>
"""
import sys, os, subprocess, urllib.request, concurrent.futures as cf

ABF = "/agent/workspace/abf"
os.chdir(ABF)
runs_file, cap, outdir = sys.argv[1], int(sys.argv[2]), sys.argv[3]
os.makedirs(outdir, exist_ok=True)


def filereport(run):
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run}"
           f"&result=read_run&fields=run_accession,library_layout,read_count,"
           f"fastq_ftp&format=tsv")
    txt = urllib.request.urlopen(url, timeout=60).read().decode().strip().split("\n")
    h = txt[0].split("\t")
    return dict(zip(h, txt[1].split("\t")))


def one(run):
    dest = f"{outdir}/{run}.fastq.gz"
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return f"have {run}"
    try:
        row = filereport(run)
    except Exception as e:
        return f"FAIL {run} meta {e}"
    urls = ["https://" + u for u in row["fastq_ftp"].split(";") if u]
    layout = row.get("library_layout", "SINGLE")
    if layout == "PAIRED" and len(urls) >= 2:
        per = cap // 2
        urls = urls[:2]
    else:
        per = cap
        urls = urls[:1]
    parts = []
    for i, u in enumerate(urls):
        p = f"{dest}.m{i}"
        cmd = (f"curl -sfL --retry 3 '{u}' | zcat 2>/dev/null | head -n {4*per} "
               f"| gzip -c > {p}")
        subprocess.run(cmd, shell=True, timeout=2400, check=False)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            parts.append(p)
    if not parts:
        return f"FAIL {run}"
    if len(parts) == 1:
        os.rename(parts[0], dest)
    else:
        with open(dest, "wb") as out:
            for p in parts:
                with open(p, "rb") as fh:
                    out.write(fh.read())
                os.remove(p)
    return f"OK {run} {os.path.getsize(dest)/1e6:.1f}MB ({layout})"


runs = [l.strip() for l in open(runs_file) if l.strip()]
with cf.ThreadPoolExecutor(max_workers=2) as ex:
    for res in ex.map(one, runs):
        print(res, flush=True)
