#!/bin/bash
cd /tmp/abf
MM2=bin/minimap2-2.28_x64-linux/minimap2
S=$1; P=$2   # sample, panel tag (pseu|aci)
F=fastq/${S}.fastq.gz
$MM2 -x sr -t 4 -a --secondary=no panel_${P}_matrix.mmi $F 2>/dev/null | python3 -c "
import sys,json,collections
c=collections.Counter()
for line in sys.stdin:
    if line.startswith('@'): continue
    f=line.split('\t',6)
    if len(f)<7 or int(f[1])&4: continue
    if int(f[4])>=20: c[f[2].split('|')[2]]+=1
json.dump(dict(c), open('${S}_${P}_matrix.json','w'))
"
$MM2 -x sr -t 4 -a --secondary=no panel_${P}_rp.mmi $F 2>/dev/null | python3 -c "
import sys,json,collections
c=collections.Counter()
for line in sys.stdin:
    if line.startswith('@'): continue
    f=line.split('\t',6)
    if len(f)<7 or int(f[1])&4: continue
    if int(f[4])>=20: c[f[2].split('|')[2]]+=1
json.dump(dict(c), open('${S}_${P}_rp.json','w'))
"
echo "$S $P done"
