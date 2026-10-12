#!/usr/bin/env python3
"""Exact positive control: the signature was derived on pediatric Definite
bacterial (n=7) vs purely-viral (n=36). Verify it separates them here."""
import json, csv, math
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

import openpyxl
wb=openpyxl.load_workbook('sdd3.xlsx',read_only=True); ws=wb['Supp_Data_1']
rows=list(ws.iter_rows(values_only=True)); h=rows[0]
itype={r[0]:r[4] for r in rows[1:]}
grp={r[0]:r[3] for r in rows[1:]}

meta={}
with open('ped_meta.csv') as f:
    for r in csv.DictReader(f): meta[r['sample_name']]=r
genes={}; samples=None
with open('ped_counts.csv') as f:
    r=csv.reader(f); hd=next(r); samples=hd[2:]
    for row in r:
        if not row: continue
        g=row[0].strip().strip('"')
        try: genes[g]=[float(x) for x in row[2:]]
        except ValueError: continue
keep=[i for i,s in enumerate(samples) if s in meta and meta[s]['is_water']=='FALSE']
samples=[samples[i] for i in keep]; genes={g:[v[i] for i in keep] for g,v in genes.items()}
C=len(samples); lib=[sum(genes[g][i] for g in genes) for i in range(C)]
l2c={g:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(genes[g])] for g in genes}
pb=[g for g in BACT if g in l2c]; pv=[g for g in VIR if g in l2c]
zb=[zc(l2c[g]) for g in pb]; zv=[zc(l2c[g]) for g in pv]
sc=[mean([zb[j][i] for j in range(len(zb))])-mean([zv[j][i] for j in range(len(zv))]) for i in range(C)]

g=defaultdict(list)
for s,v in zip(samples,sc):
    t=itype.get(s,'NA')
    if grp.get(s)=='Definite': g[t].append(v)
print("PEDIATRIC Definite arm by infection type (exact derivation contrast):")
for k in ['bacterial','viral','coinfection','NA']:
    if k in g: print(f"  {k:12s} n={len(g[k]):3d}  mean={mean(g[k]):+.3f}  median={median(g[k]):+.3f}")
if 'bacterial' in g and 'viral' in g:
    print(f"  bacterial vs purely-viral: diff(mean)={mean(g['bacterial'])-mean(g['viral']):+.3f}  "
          f"MWU p={mwu_p(g['bacterial'],g['viral']):.4f}  (n=7 vs 36, derivation contrast)")
    # also coinfection
    print(f"  coinfection vs purely-viral: diff(mean)={mean(g['coinfection'])-mean(g['viral']):+.3f}  "
          f"MWU p={mwu_p(g['coinfection'],g['viral']):.2e}")

# GSE161731 hospitalized COVID timecourse
print("\nGSE161731 hospitalized COVID timecourse check:")
import gzip
key={}
with gzip.open('GSE161731_counts_key.csv.gz','rt') as f:
    for r in csv.DictReader(f): key[r['rna_id']]=r
with gzip.open('GSE161731_counts.csv.gz','rt') as f:
    r=csv.reader(f); hd=next(r); s2=hd[1:]
    gg={}
    for row in r:
        if not row: continue
        gi=row[0].strip().strip('"')
        try: gg[gi]=[float(x) for x in row[1:]]
        except ValueError: continue
C2=len(s2); lib2=[sum(gg[x][i] for x in gg) for i in range(C2)]
l2c2={x:[math.log2(v/lib2[i]*1e6+1.0) for i,v in enumerate(gg[x])] for x in gg}
pb2=[x for x in BACT if x in l2c2]; pv2=[x for x in VIR if x in l2c2]
zb2=[zc(l2c2[x]) for x in pb2]; zv2=[zc(l2c2[x]) for x in pv2]
sc2=[mean([zb2[j][i] for j in range(len(zb2))])-mean([zv2[j][i] for j in range(len(zv2))]) for i in range(C2)]
ht=defaultdict(list)
for s,v in zip(s2,sc2):
    k=key.get(s)
    if k and k['cohort']=='COVID-19' and k['hospitalized']=='Yes' and k['time_since_onset'] in ('early','middle','late'):
        ht[k['time_since_onset']].append((k['subject_id'],v))
for t in ['early','middle','late']:
    if t in ht:
        vals=[x[1] for x in ht[t]]; subs=set(x[0] for x in ht[t])
        print(f"  {t:7s} n={len(vals):2d} samples / {len(subs)} patients  mean={mean(vals):+.3f} median={median(vals):+.3f}")
subj=defaultdict(dict)
for t in ['early','middle','late']:
    for sid,v in ht.get(t,[]): subj[sid][t]=v
pairs=[(s,d['early'],d['late'],d['late']-d['early']) for s,d in subj.items() if 'early' in d and 'late' in d]
print(f"  within-patient early->late pairs n={len(pairs)}: {[round(p[3],2) for p in pairs]}")
