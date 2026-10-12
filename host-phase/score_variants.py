#!/usr/bin/env python3
"""Diagnose why the GSE246795 baseline differs: test score-normalization variants."""
import json, gzip, csv, math, random
from collections import defaultdict
random.seed(7)

sig = json.load(open('mick108_raw.json'))
BACT = [g['ens'] for g in sig if g['logFC'] < 0]
VIR  = [g['ens'] for g in sig if g['logFC'] > 0]

def read_counts(path):
    with gzip.open(path,'rt') as f:
        r=csv.reader(f); h=next(r); samples=h[1:]; genes={}
        for row in r:
            if not row: continue
            g=row[0].strip().strip('"')
            try: genes[g]=[float(x) for x in row[1:]]
            except ValueError: continue
    return samples, genes

def median(x):
    s=sorted(x); n=len(s); return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2.0

def mwu_p(a,b):
    n1,n2=len(a),len(b); allv=sorted(a+b); rank={}; i=0
    while i<len(allv):
        j=i
        while j+1<len(allv) and allv[j+1]==allv[i]: j+=1
        rank[allv[i]]=(i+j)/2.0+1; i=j+1
    U=sum(rank[x] for x in a)-n1*(n1+1)/2.0
    mu=n1*n2/2.0
    from collections import Counter
    cnt=Counter(allv); ties=sum(c**3-c for c in cnt.values()); N=n1+n2
    sigma=math.sqrt(n1*n2/12.0*((N+1)-ties/(N*(N-1)))) if N>1 else 0
    z=(U-mu)/sigma if sigma>0 else 0
    return math.erfc(abs(z)/math.sqrt(2))

lab={}
with gzip.open('GSE246795_series_matrix.txt.gz','rt') as f:
    titles=None
    for line in f:
        if line.startswith('!Sample_title'):
            titles=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        if line.startswith('!Sample_characteristics_ch1') and '2bp_status' in line:
            vals=[x.strip().strip('"').replace('2bp_status:','').strip() for x in line.rstrip('\n').split('\t')[1:]]
            lab=dict(zip(titles,vals)); break

samples, genes = read_counts('GSE246795_host_counts.csv.gz')
C=len(samples)
lib=[sum(genes[g][i] for g in genes) for i in range(C)]
l2c={g:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(genes[g])] for g in genes}

def zc(vals):  # z across samples
    m=sum(vals)/len(vals); sd=math.sqrt(sum((x-m)**2 for x in vals)/len(vals))
    return [(x-m)/sd for x in vals] if sd>0 else [0.0]*len(vals)

present_b=[g for g in BACT if g in l2c]; present_v=[g for g in VIR if g in l2c]

# V1: mean of per-gene across-sample z   (mine)
zb=[zc(l2c[g]) for g in present_b]; zv=[zc(l2c[g]) for g in present_v]
V1=[sum(zb[j][i] for j in range(len(zb)))/len(zb) - sum(zv[j][i] for j in range(len(zv)))/len(zv) for i in range(C)]

# V2: per-sample within-108-gene z  (centred within sample)
allg=present_b+present_v
V2=[]
for i in range(C):
    vals=[l2c[g][i] for g in allg]
    m=sum(vals)/len(vals); sd=math.sqrt(sum((x-m)**2 for x in vals)/len(vals)) or 1
    V2.append(sum((l2c[g][i]-m)/sd for g in present_b)/len(present_b)
              - sum((l2c[g][i]-m)/sd for g in present_v)/len(present_v))

# V3: z after averaging within program
mb=[sum(l2c[g][i] for g in present_b)/len(present_b) for i in range(C)]
mv=[sum(l2c[g][i] for g in present_v)/len(present_v) for i in range(C)]
V3=[a-b for a,b in zip(zc(mb),zc(mv))]

# V4: no z at all, mean log2CPM difference
V4=[mb[i]-mv[i] for i in range(C)]

for name,V in [('V1 mean-of-across-sample-z',V1),
               ('V2 within-sample-z',V2),
               ('V3 z-after-program-mean',V3),
               ('V4 raw log2CPM diff',V4)]:
    g=defaultdict(list)
    for s,v in zip(samples,V): g[lab.get(s,'?')].append(v)
    a,b=g.get('2BP',[]),g.get('No-BP',[])
    print(f"{name:28s}  2BP med={median(a):+.3f} (n={len(a)})   No-BP med={median(b):+.3f} (n={len(b)})   "
          f"diff={median(a)-median(b):+.3f}   MWU p={mwu_p(a,b):.3f}")
