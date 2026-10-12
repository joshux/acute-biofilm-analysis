#!/usr/bin/env python3
"""Screen all PRJNA1033689 runs (subsampled to 40k reads) against panel_multi.

For each run: genus-level MAPQ>=20 counts among the 8 two-axis organisms +
all decoys. Flags any run with >=20 reads on a matrix-capable organism.
usage: screen_nobp.py <runs_json> [nreads]
"""
import json, subprocess, os, sys, collections

ABF = "/tmp/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
g2 = json.load(open("acc2genus.json"))
ref2g = json.load(open("ref2genus.json"))
runs = json.load(open(sys.argv[1]))
NREADS = int(sys.argv[2]) if len(sys.argv) > 2 else 40000

AXIS_ORGS = {"Pseudomonas","Klebsiella","Escherichia","Staphylococcus",
             "Acinetobacter","Haemophilus","Streptococcus","Stenotrophomonas"}

os.makedirs("fastq", exist_ok=True)
results = {}
for i, r in enumerate(runs):
    srr = r['srr']
    outf = f"screen_out/{srr}.json"
    os.makedirs("screen_out", exist_ok=True)
    if os.path.exists(outf):
        results[srr] = json.load(open(outf)); continue
    fq = f"fastq/{srr}.fastq.gz"
    if not (os.path.exists(fq) and os.path.getsize(fq) > 1000):
        urls = []
        # get fastq urls
        tsv = subprocess.run(['curl','-s', f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={srr}&result=read_run&fields=fastq_ftp&format=tsv"], capture_output=True, text=True).stdout
        lines = tsv.strip().split('\n')
        if len(lines) >= 2:
            cols = lines[1].split('\t')
            ftp = cols[1] if len(cols) > 1 else ''
            urls = ['https://' + u for u in ftp.split(';') if u.strip()]
        if not urls:
            results[srr] = {"error": "no_url"}; continue
        # subsample first NREADS reads from mate1+mate2
        parts = []
        for j, u in enumerate(urls[:2]):
            p = f"{fq}.p{j}"
            cmd = f"curl -sfL --retry 3 '{u}' | zcat 2>/dev/null | head -n {4*NREADS} | gzip -c > {p}"
            subprocess.run(cmd, shell=True, timeout=1200)
            if os.path.exists(p) and os.path.getsize(p) > 0: parts.append(p)
        if not parts or all(os.path.getsize(p) < 1000 for p in parts):
            results[srr] = {"error": "dl_fail"}
            for p in parts: os.path.exists(p) and os.remove(p)
            continue
        if len(parts) == 1: os.rename(parts[0], fq)
        else:
            with open(fq,'wb') as o:
                for p in parts: o.write(open(p,'rb').read()); os.remove(p)
    # map
    c = collections.Counter()
    cmd = f"{MM2} -x sr -t 4 -a --secondary=no panel_multi.mmi {fq} 2>/dev/null"
    with os.popen(cmd) as fh:
        for line in fh:
            if line[0] == '@': continue
            f = line.split('\t', 6)
            if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20: continue
            ref = f[2].split()[0]
            acc = ref.split('|')[0] if '|' in ref else ref
            g = g2.get(acc) or ref2g.get(ref)
            c[g or '?'] += 1
    res = {"axis_counts": {g: n for g, n in c.items() if g in AXIS_ORGS},
           "top": c.most_common(5), "total_axis": sum(c[g] for g in AXIS_ORGS)}
    json.dump(res, open(outf, 'w'))
    results[srr] = res
    flag = {g: n for g, n in res['axis_counts'].items() if n >= 20}
    print(f"{i+1}/{len(runs)} {srr} {r['spec']}/{r['label']:5s} axis={res['total_axis']:5d} flag={flag}", flush=True)

json.dump(results, open("screen_nobp_results.json","w"), indent=0)
print("DONE")
