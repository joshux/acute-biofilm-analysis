#!/usr/bin/env python3
"""
Josh's question: take COVID data that was NEVER split into 2bp/no-bp, score it on the
classifier as a whole (ignoring timing), and see what it reads.

Cohorts (all whole blood / blood, none 2BP-stratified):
  GSE157103  100 COVID + 26 non-COVID, severity gradient (ICU, APACHE II, SOFA, CRP, procalcitonin)
  GSE152641   62 COVID + 24 healthy
  GSE243217   35 COVID + 22 bacterial sepsis + 15 healthy  (PBMC - neutrophil-depleted)

Also: is there a bacterial-READING SUBSET inside COVID (bimodality), not just a mean?
"""
import gzip, csv, math, json, tarfile, io, os
from collections import defaultdict, Counter
random_seed = 20261012

sig=json.load(open('mick108_raw.json'))
BACT=[g['ens'] for g in sig if g['logFC']<0]; VIR=[g['ens'] for g in sig if g['logFC']>0]
SYM={g['ens']:g['sym'] for g in sig}

# ensembl -> symbol map from the pediatric file (has both)
ENS2SYM={}
with open('ped_counts.csv') as f:
    r=csv.reader(f); next(r)
    for row in r:
        if row and row[0].startswith('ENSG'):
            ENS2SYM[row[0].strip().strip('"')]=row[1].strip().strip('"')
BACT_S=[ENS2SYM[g] for g in BACT if g in ENS2SYM]
VIR_S =[ENS2SYM[g] for g in VIR  if g in ENS2SYM]
print(f"signature: {len(BACT)} bacterial ({len(BACT_S)} with symbols), {len(VIR)} viral ({len(VIR_S)} with symbols)")

def mean(x): return sum(x)/len(x)
def median(x):
    s=sorted(x); n=len(s); return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2.0
def quant(x,q):
    s=sorted(x); i=int(q*(len(s)-1)); return s[i]
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
def spearman(x,y):
    pts=[(a,b) for a,b in zip(x,y) if a is not None and b is not None]
    if len(pts)<4: return float('nan'), len(pts)
    x=[p[0] for p in pts]; y=[p[1] for p in pts]; n=len(x)
    def rk(v):
        order=sorted(range(n),key=lambda i:v[i]); r=[0]*n
        for pos,i in enumerate(order): r[i]=pos+1
        return r
    rx,ry=rk(x),rk(y); mx,my=mean(rx),mean(ry)
    num=sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    den=math.sqrt(sum((rx[i]-mx)**2 for i in range(n))*sum((ry[i]-my)**2 for i in range(n)))
    return (num/den if den else float('nan')), n

def score_symbols(genes_by_symbol, samples, bsyms, vsyms):
    """genes_by_symbol: {symbol: [values per sample]}"""
    C=len(samples)
    lib=[sum(genes_by_symbol[s][i] for s in genes_by_symbol) for i in range(C)]
    l2c={s:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(genes_by_symbol[s])] for s in genes_by_symbol}
    pb=[s for s in bsyms if s in l2c]; pv=[s for s in vsyms if s in l2c]
    zb=[zc(l2c[s]) for s in pb]; zv=[zc(l2c[s]) for s in pv]
    sc=[mean([zb[j][i] for j in range(len(zb))])-mean([zv[j][i] for j in range(len(zv))]) for i in range(C)]
    return sc, len(pb), len(pv)

# ---------------- GSE161731 anchors (whole blood, has bacterial arm) ----------------
print("\n"+"="*78); print("ANCHORS — GSE161731 (whole blood; has a real bacterial arm)"); print("="*78)
key={}
with gzip.open('GSE161731_counts_key.csv.gz','rt') as f:
    for r in csv.DictReader(f): key[r['rna_id']]=r
gb={}; s731=None
with gzip.open('GSE161731_counts.csv.gz','rt') as f:
    r=csv.reader(f); h=next(r); s731=h[1:]
    for row in r:
        if not row: continue
        g=row[0].strip().strip('"')
        try: gb[g]=[float(x) for x in row[1:]]
        except ValueError: continue
C=len(s731); lib=[sum(gb[g][i] for g in gb) for i in range(C)]
l2c={g:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(gb[g])] for g in gb}
pb=[g for g in BACT if g in l2c]; pv=[g for g in VIR if g in l2c]
zb=[zc(l2c[g]) for g in pb]; zv=[zc(l2c[g]) for g in pv]
sc731=[mean([zb[j][i] for j in range(len(zb))])-mean([zv[j][i] for j in range(len(zv))]) for i in range(C)]
g731=defaultdict(list)
for s,v in zip(s731,sc731):
    k=key.get(s)
    if k: g731[k['cohort']].append(v)
for c in ['Bacterial','Influenza','CoV other','COVID-19','healthy']:
    if c in g731:
        v=g731[c]; print(f"  {c:12s} n={len(v):3d}  median={median(v):+.3f}  mean={mean(v):+.3f}  Q25={quant(v,.25):+.3f} Q75={quant(v,.75):+.3f}")
bac_lo=quant(g731['Bacterial'],.25); bac_hi=quant(g731['Bacterial'],.75)
print(f"  --> bacterial-arm IQR: [{bac_lo:+.3f}, {bac_hi:+.3f}]")
for c in ['COVID-19','Influenza','healthy']:
    if c in g731:
        frac=sum(1 for x in g731[c] if x>=bac_lo)/len(g731[c])
        print(f"      fraction of {c} at/above bacterial Q25: {frac:.2f}")

# ---------------- GSE157103 ----------------
print("\n"+"="*78); print("GSE157103 — 100 COVID + 26 non-COVID, whole blood, severity gradient"); print("="*78)
# series metadata
sm={}
with gzip.open('GSE157103_series.txt.gz','rt',errors='replace') as f:
    titles=None; chars=[]
    for line in f:
        if line.startswith('!Sample_title'):
            titles=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            chars.append([x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]])
    n=len(titles)
    meta=[{} for _ in range(n)]
    for ch in chars:
        for i in range(n):
            if ':' in ch[i]:
                k,v=ch[i].split(':',1); meta[i][k.strip().lower()]=v.strip()
    for i,t in enumerate(titles): meta[i]['title']=t
# expression (symbol x C<n>)
gs={}; cols=[]
with gzip.open('GSE157103_ec.tsv.gz','rt',errors='replace') as f:
    r=csv.reader(f, delimiter='\t'); h=next(r)
    cols=h[1:]
    for row in r:
        if not row or len(row)<2: continue
        s=row[0].strip().strip('"')
        try: gs[s]=[float(x) for x in row[1:]]
        except ValueError: continue
# map columns C<n> -> COVID sample n-1 ; NC<n> -> NONCOVID sample n-1 (after the COVID block)
import re as _re
def col_to_idx(c):
    m=_re.match(r'^([A-Za-z]+)(\d+)$', c.strip())
    if not m: return None
    pre,num=m.group(1).upper(), int(m.group(2))
    if pre=='C':  return num-1
    if pre=='NC': return 100+(num-1)
    return None
idx=[col_to_idx(c) for c in cols]
bad=[c for c,i in zip(cols,idx) if i is None or not (0<=i<len(titles))]
print(f"expression: {len(gs)} genes x {len(cols)} cols; sample titles {len(titles)}; unmapped cols: {bad[:5]}")
keep=[j for j,i in enumerate(idx) if 0<=i<len(titles)]
gs2={s:[v[j] for j in keep] for s,v in gs.items()}
samp=[meta[idx[j]] for j in keep]
sc,nb,nv=score_symbols(gs2, samp, BACT_S, VIR_S)
print(f"coverage: bacterial {nb}/{len(BACT_S)}  viral {nv}/{len(VIR_S)}  samples {len(samp)}")
by=defaultdict(list)
for m,v in zip(samp,sc): by[m.get('disease state','?')].append(v)
for k in by: print(f"  {k:14s} n={len(by[k]):3d}  median={median(by[k]):+.3f}  mean={mean(by[k]):+.3f}")
if 'COVID-19' in by and 'non-COVID-19' in by:
    print(f"  COVID vs non-COVID: p={mwu_p(by['COVID-19'],by['non-COVID-19']):.4f}")

cov=[(m,v) for m,v in zip(samp,sc) if m.get('disease state')=='COVID-19']
print("\n  -- within COVID: score vs severity (the confound test) --")
def num(x):
    try: return float(x)
    except: return None
for field in ['icu','apacheii','sofa','crp (mg/l)','procalcitonin (ng/ml)','ferritin (ng/ml)','ddimer (mg/l_feu)','mechanical ventilation']:
    if field=='icu' or field=='mechanical ventilation':
        yes=[v for m,v in cov if m.get(field,'').lower()=='yes']
        no=[v for m,v in cov if m.get(field,'').lower()=='no']
        if yes and no:
            print(f"    {field:26s} yes={median(yes):+.3f}(n{len(yes)})  no={median(no):+.3f}(n{len(no)})  p={mwu_p(yes,no):.4f}")
    else:
        rho,nn=spearman([num(m.get(field)) for m,v in cov],[v for m,v in cov])
        print(f"    {field:26s} Spearman rho={rho:+.3f}  (n={nn})")

# bimodality / subset test
v=by['COVID-19']
print(f"\n  -- distribution shape of COVID arm (bacterial-reading subset?) --")
print(f"    min={min(v):+.3f}  Q10={quant(v,.10):+.3f}  Q25={quant(v,.25):+.3f}  median={median(v):+.3f}  Q75={quant(v,.75):+.3f}  Q90={quant(v,.90):+.3f}  max={max(v):+.3f}")
print(f"    fraction of COVID at/above GSE161731 bacterial-arm Q25 ({bac_lo:+.3f}): {sum(1 for x in v if x>=bac_lo)/len(v):.2f}")
print(f"    fraction of COVID at/above GSE161731 bacterial-arm median: {sum(1 for x in v if x>=median(g731['Bacterial']))/len(v):.2f}")

# ---------------- GSE152641 ----------------
print("\n"+"="*78); print("GSE152641 — 62 COVID + 24 healthy, whole blood"); print("="*78)
titles=None; chars=[]
with gzip.open('GSE152641_series.txt.gz','rt',errors='replace') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            titles=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            chars.append([x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]])
n=len(titles); meta=[{} for _ in range(n)]
for ch in chars:
    for i in range(n):
        if ':' in ch[i]:
            k,vv=ch[i].split(':',1); meta[i][k.strip().lower()]=vv.strip()
for i,t in enumerate(titles): meta[i]['title']=t
# counts are by Entrez -> need symbol. Try: the file's first column is Entrez id
with gzip.open('GSE152641_counts.csv.gz','rt',errors='replace') as f:
    r=csv.reader(f); h=next(r); s2=h[1:]
print(f"counts: {len(s2)} samples, first col values e.g.", end=' ')
row1=None
with gzip.open('GSE152641_counts.csv.gz','rt',errors='replace') as f:
    r=csv.reader(f); next(r); row1=next(r)
print(row1[0], "(Entrez id - needs mapping; skipped in favour of symbol-based cohorts)")
