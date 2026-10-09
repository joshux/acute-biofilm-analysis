#!/usr/bin/env python3
"""Download FASTQ from ENA using the paths ENA itself reports.

usage: get_fastq.py <accession|SRR...> [outdir] [max_gb_per_run]
 - if given a study accession, fetches every RNA-Seq run
 - streams only the first N reads when --subsample N is given (healthy controls)
"""
import sys, os, json, urllib.request, subprocess, gzip, concurrent.futures as cf

ABF = "/agent/workspace/abf"
os.chdir(ABF)
os.makedirs("fastq", exist_ok=True)


def filereport(acc):
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}"
           f"&result=read_run&fields=run_accession,library_strategy,library_layout,"
           f"read_count,fastq_ftp,fastq_bytes&format=tsv")
    txt = urllib.request.urlopen(url, timeout=60).read().decode().strip().split("\n")
    hdr = txt[0].split("\t")
    rows = [dict(zip(hdr, l.split("\t"))) for l in txt[1:] if l.strip()]
    return rows


def ena_urls(row):
    return ["https://" + u for u in row["fastq_ftp"].split(";") if u]


def subsample_download(urls, dest, nreads):
    """stream first nreads from the first URL (single-end) or both mates"""
    parts = []
    for i, u in enumerate(urls):
        p = f"{dest}.part{i}"
        # zcat | head -n (4*nreads) -> gzip
        cmd = (f"curl -sfL --retry 3 '{u}' | zcat 2>/dev/null | head -n {4*nreads} "
               f"| gzip -c > {p}")
        subprocess.run(cmd, shell=True, timeout=1200, check=False)
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            return False
        parts.append(p)
    return parts


def fetch(spec):
    run, urls, dest, nreads = spec
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return f"have {run}"
    if nreads:
        # for paired runs, concatenate both mates' subsamples into one file
        tmp = []
        for i, u in enumerate(urls):
            p = f"{dest}.m{i}.gz"
            cmd = (f"curl -sfL --retry 3 '{u}' | zcat 2>/dev/null | head -n {4*nreads} "
                   f"| gzip -c > {p}")
            subprocess.run(cmd, shell=True, timeout=1800, check=False)
            if os.path.exists(p) and os.path.getsize(p) > 0:
                tmp.append(p)
        if not tmp:
            return f"FAIL {run}"
        if len(tmp) == 1:
            os.rename(tmp[0], dest)
        else:
            # merge mate subsamples (interleaved-ish: just cat the gz streams)
            with open(dest, "wb") as out:
                for p in tmp:
                    with open(p, "rb") as fh:
                        out.write(fh.read())
                    os.remove(p)
        return f"OK {run} {os.path.getsize(dest)/1e6:.1f}MB (sub {nreads})"
    # full download
    cmd = f"curl -sfL --retry 3 -o {dest} '{urls[0]}'"
    subprocess.run(cmd, shell=True, timeout=1800, check=False)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return f"OK {run} {os.path.getsize(dest)/1e6:.1f}MB"
    return f"FAIL {run}"


if __name__ == "__main__":
    acc = sys.argv[1]
    sub_n = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    rows = filereport(acc)
    # single run accession?
    if len(rows) == 1 and rows[0]["run_accession"] == acc:
        runs = rows
    else:
        runs = [r for r in rows if r.get("library_strategy") == "RNA-Seq"]
    print(f"{acc}: {len(runs)} RNA-Seq runs")
    specs = []
    for r in runs:
        dest = f"fastq/{r['run_accession']}.fastq.gz"
        specs.append((r["run_accession"], ena_urls(r), dest, sub_n))
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        for res in ex.map(fetch, specs):
            print(res, flush=True)
