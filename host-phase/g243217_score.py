#!/usr/bin/env python3
"""GSE243217: 35 COVID vs 22 bacterial sepsis vs 15 healthy (PBMC).
Direct test of whether the 'bacterial' program is higher in COVID than in true bacterial infection."""
import gzip, csv, glob, math, os
from collections import defaultdict, Counter
sig=json.load(open('mick108_raw.json')) if False else None
import json
sig=json.load(open('mick108_raw.json'))
BACT=[g['ens'] for g in sig if g['logFC']<0]; VIR=[g['ens'] for g in sig if g['logFC']>0]
ENS2SYM={}
with open('ped_counts.csv') as f:
    r=csv.reader(f); next(r)
    for row in r:
        if row and row[0].startswith('ENSG'): ENS2SYM[row[0].strip().strip('"')]=row[1].strip().strip('"')
BACT_S=[ENS2SYM[g] for g in BACT if g in ENS2SYM]; VIR_S=[ENS2SYM[g] for g in VIR if g in ENS2SYM]
def mean(x): return sum(x)/len(x)
def median(x):
    s=sorted(x); n=len(s); return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2.0
def mwu_p(a,b):
    n1,n2=len(a),len(b)
    if not a or not b: return float('nan')
    allv=sorted(a+b); rank={}; i=0
    while i<len(allv):
        j=i
        while j+1<len(allv) and allv[j+1]==allv[i]: j+=1
        rank[allv[i]]=(i+j)/2.0+1; i=j+1
    U=sum(rank[x] for x in a)-n1*(n1+1)/2.0; mu=n1*n2/2.0
    cnt=Counter(allv); ties=sum(c**3-c for c in cnt.values()); N=n1+n2
    sigma=math.sqrt(n1*n2/12.0*((N+1)-ties/(N*(N-1)))) if N>1 else 0
    return math.erfc(abs((U-mu)/sigma)/math.sqrt(2)) if sigma>0 else 1.0
def zc(v):
    m=mean(v); sd=math.sqrt(sum((x-m)**2 for x in v)/len(v))
    return [(x-m)/sd for x in v] if sd>0 else [0.0]*len(v)

# build matrix
files=sorted(glob.glob('g243217/*.gz'))
cols={}; names=[]
for fp in files:
    with gzip.open(fp,'rt',errors='replace') as f:
        r=csv.reader(f, delimiter='\t'); h=next(r)
        nm=h[1].replace('_raw_count','')
        names.append(nm); vals={}
        for row in r:
            if len(row)<2: continue
            try: vals[row[0].strip()]=float(row[1])
            except ValueError: continue
        cols[nm]=vals
print(f"loaded {len(names)} samples; genes/file ~{len(cols[names[0]])}")

# disease state from series
grp={}
with gzip.open('GSE243217_series.txt.gz','rt',errors='replace') as f:
    titles=None; chars=[]
    for line in f:
        if line.startswith('!Sample_title'): titles=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'): chars.append([x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]])
    for ch in chars:
        if any('disease state' in x for x in ch):
            for i,t in enumerate(titles):
                if ':' in ch[i]:
                    nm=t.split(' ')[0]
                    grp[nm]=ch[i].split(':',1)[1].strip()
print("groups:", Counter(grp.values()))

# matrix by symbol
allsym=set()
for nm in names: allsym |= set(cols[nm].keys())
allsym=sorted(allsym)
M={s:[cols[nm].get(s,0.0) for nm in names] for s in allsym}
C=len(names); lib=[sum(M[s][i] for s in M) for i in range(C)]
l2c={s:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(M[s])] for s in M}
pb=[s for s in BACT_S if s in l2c]; pv=[s for s in VIR_S if s in l2c]
print(f"coverage: bacterial {len(pb)}/{len(BACT_S)}  viral {len(pv)}/{len(VIR_S)}")
zb=[zc(l2c[s]) for s in pb]; zv=[zc(l2c[s]) for s in pv]
sc=[mean([zb[j][i] for j in range(len(zb))])-mean([zv[j][i] for j in range(len(zv))]) for i in range(C)]
g=defaultdict(list)
for nm,v in zip(names,sc):
    k=grp.get(nm,'?')
    if k!='?': g[k].append(v)
print("\nGSE243217 (PBMC) score by disease state:")
for k in ['COVID-19','Sepsis','healthy donor']:
    if k in g: print(f"  {k:14s} n={len(g[k]):3d}  median={median(g[k]):+.3f}  mean={mean(g[k]):+.3f}")
for a,b in [('COVID-19','Sepsis'),('Sepsis','healthy donor'),('COVID-19','healthy donor')]:
    if a in g and b in g:
        print(f"  {a} vs {b}: p={mwu_p(g[a],g[b]):.4f}  median diff={median(g[a])-median(g[b]):+.3f}")
# component check
print("\n  components:")
for comp,ids in [('bacterial',pb),('viral',pv)]:
    zz=[zc(l2c[s]) for s in ids]
    line=f"    {comp:10s}: "
    for k in ['COVID-19','Sepsis','healthy donor']:
        if k in g:
            vals=[mean([zz[j][i] for j in range(len(zz))]) for i,nm in enumerate(names) if grp.get(nm)==k]
            line+=f"{k}={mean(vals):+.3f} "
    print(line)
