#!/bin/bash
# usage: run_sample.sh SRRxxxxx
cd /tmp/abf
S=$1
MM2=bin/minimap2-2.28_x64-linux/minimap2
F=fastq/${S}.fastq.gz
[ -s "$F" ] || { echo "$S MISSING"; exit 1; }
# map to multi-strain genome panel (for taxonomy/genus) + matrix + growth + RP panels
$MM2 -x sr -t 4 -a --secondary=no panel_multi.mmi $F 2>/dev/null | python3 -c "
import sys,pysam,collections,json
import io
# read SAM from stdin, count MAPQ>=20 by ref
counts=collections.Counter(); tot=0
for line in sys.stdin:
    if line.startswith('@'): continue
    tot+=1
    f=line.split('\t',6)
    if len(f)<7: continue
    flag=int(f[1]); rn=f[2]; mq=int(f[4])
    if flag&4: continue
    if mq>=20: counts[rn]+=1
json.dump({'total':tot,'mq20_by_ref':dict(counts)}, open('${S}_genome.json','w'))
"
$MM2 -x sr -t 4 -a --secondary=no panel_matrix.mmi $F 2>/dev/null | python3 -c "
import sys,json,collections
c=collections.Counter(); tot=0
for line in sys.stdin:
    if line.startswith('@'): continue
    tot+=1; f=line.split('\t',6)
    if len(f)<7: continue
    if int(f[1])&4: continue
    if int(f[4])>=20:
        p=f[2].split('|'); c[p[2]]+=1
json.dump({'total':tot,'matrix_genes':dict(c)}, open('${S}_matrix.json','w'))
"
$MM2 -x sr -t 4 -a --secondary=no panel_rp.mmi $F 2>/dev/null | python3 -c "
import sys,json,collections
c=collections.Counter(); tot=0
for line in sys.stdin:
    if line.startswith('@'): continue
    tot+=1; f=line.split('\t',6)
    if len(f)<7: continue
    if int(f[1])&4: continue
    if int(f[4])>=20:
        p=f[2].split('|'); c[p[1]]+=1
json.dump({'total':tot,'rp_by_acc':dict(c)}, open('${S}_rp.json','w'))
"
echo "$S OK"
