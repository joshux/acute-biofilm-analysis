#!/bin/bash
cd /tmp/abf
MM2=bin/minimap2-2.28_x64-linux/minimap2
S=$1
F=fastq/${S}.fastq.gz
# matrix at MAPQ>=0 and >=20
$MM2 -x sr -t 4 -a --secondary=no panel_matrix.mmi $F 2>/dev/null | python3 -c "
import sys,json,collections
c0=collections.Counter(); c20=collections.Counter()
for line in sys.stdin:
    if line.startswith('@'): continue
    f=line.split('\t',6)
    if len(f)<7 or int(f[1])&4: continue
    g=f[2].split('|')[2]
    c0[g]+=1
    if int(f[4])>=20: c20[g]+=1
json.dump({'mq0':dict(c0),'mq20':dict(c20)}, open('${S}_mat2.json','w'))
"
$MM2 -x sr -t 4 -a --secondary=no panel_rp.mmi $F 2>/dev/null | python3 -c "
import sys,json,collections
c0=collections.Counter(); c20=collections.Counter()
for line in sys.stdin:
    if line.startswith('@'): continue
    f=line.split('\t',6)
    if len(f)<7 or int(f[1])&4: continue
    a=f[2].split('|')[1]
    c0[a]+=1
    if int(f[4])>=20: c20[a]+=1
json.dump({'mq0':dict(c0),'mq20':dict(c20)}, open('${S}_rp2.json','w'))
"
echo "$S rescored"
