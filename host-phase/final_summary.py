#!/usr/bin/env python3
"""Final consolidated summary for the write-up."""
import json, gzip, csv, math
from collections import defaultdict, Counter
sig=json.load(open('mick108_raw.json'))
BACT=[g['ens'] for g in sig if g['logFC']<0]; VIR=[g['ens'] for g in sig if g['logFC']>0]
SYM={g['ens']:g['sym'] for g in sig}
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
def score(genes, samples, gids_b, gids_v):
    C=len(samples); lib=[sum(genes[g][i] for g in genes) for i in range(C)]
    l2c={g:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(genes[g])] for g in genes}
    pb=[g for g in gids_b if g in l2c]; pv=[g for g in gids_v if g in l2c]
    zb=[zc(l2c[g]) for g in pb]; zv=[zc(l2c[g]) for g in pv]
    sc=[mean([zb[j][i] for j in range(len(zb))])-mean([zv[j][i] for j in range(len(zv))]) for i in range(C)]
    return sc, len(pb), len(pv)

print("="*78); print("CONSOLIDATED HOST-PHASE RESULTS"); print("="*78)

# ---- GSE161731 ----
key={}
with gzip.open('GSE161731_counts_key.csv.gz','rt') as f:
    for r in csv.DictReader(f): key[r['rna_id']]=r
with gzip.open('GSE161731_counts.csv.gz','rt') as f:
    r=csv.reader(f); hd=next(r); s2=hd[1:]; gg={}
    for row in r:
        if not row: continue
        gi=row[0].strip().strip('"')
        try: gg[gi]=[float(x) for x in row[1:]]
        except ValueError: continue
sc2,nb,nv=score(gg,s2,BACT,VIR)
rows=[]
for s,v in zip(s2,sc2):
    k=key.get(s)
    if k: rows.append((k['cohort'],k['time_since_onset'],k['hospitalized'],v))
g=defaultdict(list)
for c,t,h,v in rows: g[c].append(v)
print("\n[GSE161731] adult whole blood — score by cohort")
for c in ['Bacterial','Influenza','CoV other','COVID-19','healthy']:
    if c in g:
        lo=median(g[c]); print(f"   {c:12s} n={len(g[c]):3d}  median={lo:+.3f}  mean={mean(g[c]):+.3f}")
print("   Bacterial vs COVID-19 : p=%.2e" % mwu_p(g['Bacterial'],g['COVID-19']))
print("   Bacterial vs healthy  : p=%.2e" % mwu_p(g['Bacterial'],g['healthy']))
print("   COVID-19  vs healthy  : p=%.4f" % mwu_p(g['COVID-19'],g['healthy']))

print("\n[GSE161731] COVID-19 by time since onset (all + hospitalized-only)")
for label,flt in [('all',None),('hospitalized',{'Yes'})]:
    gt=defaultdict(list)
    for c,t,h,v in rows:
        if c=='COVID-19' and t in ('early','middle','late') and (flt is None or h in flt):
            gt[t].append(v)
    line=f"   {label:13s}: "
    for t in ['early','middle','late']:
        if t in gt: line+=f"{t}={mean(gt[t]):+.3f}(n{len(gt[t])}) "
    print(line)

# ---- GSE246795 ----
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
        gi=row[0].strip().strip('"')
        try: genes[gi]=[float(x) for x in row[1:]]
        except ValueError: continue
sc,nb,nv=score(genes,cols,BACT,VIR)
lmap=dict(zip(titles,vals))
by_name=defaultdict(list); by_pos=defaultdict(list)
for i,v in enumerate(sc):
    by_name[lmap[cols[i]]].append(v); by_pos[vals[i]].append(v)
print(f"\n[GSE246795] 2BP-labelled COVID tracheal aspirate (coverage {nb}/83, {nv}/25)")
print("   JOIN BY NAME  (correct):   2BP mean=%+.3f  No-BP mean=%+.3f  p=%.4f" % (
    mean(by_name['2BP']),mean(by_name['No-BP']),mwu_p(by_name['2BP'],by_name['No-BP'])))
print("   JOIN POSITIONALLY (bug):   2BP mean=%+.3f  No-BP mean=%+.3f  p=%.4f" % (
    mean(by_pos['2BP']),mean(by_pos['No-BP']),mwu_p(by_pos['2BP'],by_pos['No-BP'])))
print("   COMMITTED (RESULTS-HOST-PHASE.md): 2BP +0.243  No-BP +0.192  p=0.87")
print("   label agreement between joins: %d/56" % sum(1 for i in range(len(cols)) if lmap[cols[i]]==vals[i]))

# top contributing genes in the corrected 2BP contrast
print("\n   top genes driving corrected 2BP>no-BP (mean log2CPM diff):")
C=len(cols); lib=[sum(genes[g][i] for g in genes) for i in range(C)]
l2c={g:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(genes[g])] for g in genes}
ib=[i for i,s in enumerate(cols) if lmap[s]=='2BP']; iv=[i for i,s in enumerate(cols) if lmap[s]=='No-BP']
d=[]
for gid in BACT:
    if gid in l2c: d.append((mean([l2c[gid][i] for i in ib])-mean([l2c[gid][i] for i in iv]), SYM.get(gid,gid)))
for x,s in sorted(d,key=lambda z:-z[0])[:10]: print(f"     {s:10s} {x:+.2f}")
print("   viral program top (negative = lower in 2BP):")
d=[]
for gid in VIR:
    if gid in l2c: d.append((mean([l2c[gid][i] for i in ib])-mean([l2c[gid][i] for i in iv]), SYM.get(gid,gid)))
for x,s in sorted(d,key=lambda z:z[0])[:8]: print(f"     {s:10s} {x:+.2f}")
