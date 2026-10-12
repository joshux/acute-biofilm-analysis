#!/usr/bin/env python3
"""Test the misalignment hypothesis: counts-file column order != series-matrix sample order.
If the earlier run joined labels POSITIONALLY, labels are scrambled -> null."""
import json, gzip, csv, math
from collections import defaultdict, Counter
sig=json.load(open('mick108_raw.json'))
BACT=[g['ens'] for g in sig if g['logFC']<0]; VIR=[g['ens'] for g in sig if g['logFC']>0]
def mean(x): return sum(x)/len(x)
def median(x):
    s=sorted(x); n=len(s); return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2.0
def mwu_p(a,b):
    n1,n2=len(a),len(b); allv=sorted(a+b); rank={}; i=0
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

# series matrix order + labels
titles=None; vals=None
with gzip.open('GSE246795_series_matrix.txt.gz','rt') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            titles=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        if line.startswith('!Sample_characteristics_ch1') and '2bp_status' in line:
            vals=[x.strip().strip('"').replace('2bp_status:','').strip() for x in line.rstrip('\n').split('\t')[1:]]
            break

genes={}; cols=None
with gzip.open('GSE246795_host_counts.csv.gz','rt') as f:
    r=csv.reader(f); h=next(r); cols=h[1:]
    for row in r:
        if not row: continue
        g=row[0].strip().strip('"')
        try: genes[g]=[float(x) for x in row[1:]]
        except ValueError: continue

print(f"counts columns: {len(cols)}   series-matrix samples: {len(titles)}")
print(f"orders identical? {cols == titles}")
print(f"first 4 counts cols : {cols[:4]}")
print(f"first 4 matrix cols : {titles[:4]}")

C=len(cols); lib=[sum(genes[g][i] for g in genes) for i in range(C)]
l2c={g:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(genes[g])] for g in genes}
pb=[g for g in BACT if g in l2c]; pv=[g for g in VIR if g in l2c]
zb=[zc(l2c[g]) for g in pb]; zv=[zc(l2c[g]) for g in pv]
sc=[mean([zb[j][i] for j in range(len(zb))])-mean([zv[j][i] for j in range(len(zv))]) for i in range(C)]

print("\n--- (A) JOIN BY NAME (correct) ---")
by_name=defaultdict(list)
for s,v in zip(cols,sc): by_name[dict(zip(titles,vals)).get(s,'?')].append(v)
a,b=by_name['2BP'],by_name['No-BP']
print(f"  2BP mean={mean(a):+.3f}  No-BP mean={mean(b):+.3f}  MWU p={mwu_p(a,b):.4f}")

print("\n--- (B) JOIN POSITIONALLY (misalignment hypothesis) ---")
by_pos=defaultdict(list)
for i,v in enumerate(sc): by_pos[vals[i]].append(v)   # label from matrix order, value from counts order
a2,b2=by_pos['2BP'],by_pos['No-BP']
print(f"  2BP mean={mean(a2):+.3f}  No-BP mean={mean(b2):+.3f}  MWU p={mwu_p(a2,b2):.4f}")
print(f"  (committed null was: 2BP +0.243  No-BP +0.192  p=0.87)")

# how scrambled? match fraction
print(f"\n  label agreement between the two joins: {sum(1 for i in range(C) if dict(zip(titles,vals))[cols[i]]==vals[i])}/{C}")
