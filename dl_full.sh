#!/bin/bash
# download full run (both mates) for one SRR
S=$1
D=/tmp/abf/fastq
mkdir -p $D
python3 - <<PY
import subprocess, os
srr = "$S"
tsv = subprocess.run(['curl','-s', f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={srr}&result=read_run&fields=fastq_ftp&format=tsv"], capture_output=True, text=True).stdout
lines = tsv.strip().split('\n')
urls = []
if len(lines) >= 2:
    cols = lines[1].split('\t')
    ftp = cols[1] if len(cols) > 1 else ''
    urls = ['https://' + u for u in ftp.split(';') if u.strip()]
ok = []
for j, u in enumerate(urls):
    out = f"fastq/{srr}_{j+1}.fastq.gz"
    if not (os.path.exists(out) and os.path.getsize(out) > 100000):
        subprocess.run(['curl','-sfL','--retry','3','-o',out, u], timeout=3000)
    ok.append((out, os.path.getsize(out) if os.path.exists(out) else 0))
print('DL', srr, ok)
PY
