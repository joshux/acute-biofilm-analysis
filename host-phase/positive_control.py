#!/usr/bin/env python3
"""
Resolve the conflict: my score separates 2BP vs no-BP in GSE246795 at p=0.001,
but the committed result said p=0.87. Run (a) a POSITIVE CONTROL in the derivation
cohort (pediatric, Definite vs No Evidence), and (b) confound checks in GSE246795
(depth, patient identity, gene coverage, program behaviour).
"""
import json, gzip, csv, math, random
from collections import defaultdict, Counter
random.seed(11)

sig = json.load(open('mick108_raw.json'))
BACT=[g['ens'] for g in sig if g['logFC']<0]
VIR =[g['ens'] for g in sig if g['logFC']>0]

def median(x):
    s=sorted(x); n=len(s); return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2.0
def mwu_p(a,b):
    n1,n2=len(a),len(b)
    if n1==0 or n2==0: return float('nan')
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
    m=sum(v)/len(v); sd=math.sqrt(sum((x-m)**2 for x in v)/len(v))
    return [(x-m)/sd for x in v] if sd>0 else [0.0]*len(v)

def score_matrix(genes, samples, gids_b, gids_v):
    C=len(samples)
    lib=[sum(genes[g][i] for g in genes) for i in range(C)]
    l2c={g:[math.log2(v/lib[i]*1e6+1.0) for i,v in enumerate(genes[g])] for g in genes}
    pb=[g for g in gids_b if g in l2c]; pv=[g for g in gids_v if g in l2c]
    zb=[zc(l2c[g]) for g in pb]; zv=[zc(l2c[g]) for g in pv]
    sc=[sum(zb[j][i] for j in range(len(zb)))/len(zb)
        - sum(zv[j][i] for j in range(len(zv)))/len(zv) for i in range(C)]
    return sc, len(pb), len(pv), l2c, lib

# ============ POSITIVE CONTROL: pediatric derivation cohort ============
print("="*70); print("POSITIVE CONTROL — pediatric cohort (signature derived here)"); print("="*70)
meta={}
with open('ped_meta.csv') as f:
    for r in csv.DictReader(f): meta[r['sample_name']]=r
# counts: col0 = ensembl, col1 = symbol, rest = samples
genes={}; samples=None
with open('ped_counts.csv') as f:
    r=csv.reader(f); h=next(r); samples=h[2:]
    for row in r:
        if not row: continue
        g=row[0].strip().strip('"')
        try: genes[g]=[float(x) for x in row[2:]]
        except ValueError: continue
keep=[i for i,s in enumerate(samples) if s in meta and meta[s]['is_water']=='FALSE']
samples=[samples[i] for i in keep]
genes={g:[v[i] for i in keep] for g,v in genes.items()}
print(f"pediatric: {len(samples)} non-water samples, {len(genes)} genes")
sc, nb, nv, l2c, lib = score_matrix(genes, samples, BACT, VIR)
print(f"gene coverage: bacterial {nb}/83  viral {nv}/25")
g=defaultdict(list)
for s,v in zip(samples,sc): g[meta[s]['LRTI_adjudication']].append(v)
for k in ['Definite','Suspected','Indeterminate','No Evidence']:
    if k in g: print(f"  {k:14s} n={len(g[k]):3d}  median={median(g[k]):+.3f}")
if 'Definite' in g and 'No Evidence' in g:
    print(f"  Definite vs No Evidence: diff={median(g['Definite'])-median(g['No Evidence']):+.3f} "
          f"MWU p={mwu_p(g['Definite'],g['No Evidence']):.2e}")

# ============ CONFOUND CHECKS: GSE246795 ============
print("\n"+"="*70); print("CONFOUND CHECKS — GSE246795"); print("="*70)
lab={}
with gzip.open('GSE246795_series_matrix.txt.gz','rt') as f:
    titles=None
    for line in f:
        if line.startswith('!Sample_title'):
            titles=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        if line.startswith('!Sample_characteristics_ch1') and '2bp_status' in line:
            lab=dict(zip(titles,[x.strip().strip('"').replace('2bp_status:','').strip()
                                 for x in line.rstrip('\n').split('\t')[1:]])); break
genes={}; samples=None
with gzip.open('GSE246795_host_counts.csv.gz','rt') as f:
    r=csv.reader(f); h=next(r); samples=h[1:]
    for row in r:
        if not row: continue
        g=row[0].strip().strip('"')
        try: genes[g]=[float(x) for x in row[1:]]
        except ValueError: continue
sc, nb, nv, l2c, lib = score_matrix(genes, samples, BACT, VIR)
print(f"gene coverage: bacterial {nb}/83  viral {nv}/25")
meta2=[{'s':s,'grp':lab.get(s,'?'),'score':v,'lib':lib[i]} for i,(s,v) in enumerate(zip(samples,sc))]
a=[m['score'] for m in meta2 if m['grp']=='2BP']; b=[m['score'] for m in meta2 if m['grp']=='No-BP']

# depth confound
print("\n-- library size (total host counts) by group --")
for grp in ['2BP','No-BP']:
    L=[m['lib'] for m in meta2 if m['grp']==grp]
    print(f"  {grp:6s} median={median(L):,.0f}")
la=[m['lib'] for m in meta2 if m['grp']=='2BP']; lb=[m['lib'] for m in meta2 if m['grp']=='No-BP']
print(f"  depth MWU p={mwu_p(la,lb):.3f}")

# rank corr score vs depth
def spearman(x,y):
    def rk(v):
        s=sorted(range(len(v)),key=lambda i:v[i]); r=[0]*len(v)
        for pos,i in enumerate(s): r[i]=pos+1
        return r
    rx,ry=rk(x),rk(y); n=len(x)
    mx,my=sum(rx)/n,sum(ry)/n
    num=sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    den=math.sqrt(sum((rx[i]-mx)**2 for i in range(n))*sum((ry[i]-my)**2 for i in range(n)))
    return num/den if den else 0
print(f"\n-- score vs library size: Spearman rho={spearman([m['score'] for m in meta2],[m['lib'] for m in meta2]):+.3f}")

# per-gene: fraction of 83 bacterial genes with higher mean in 2BP; and top genes
print("\n-- per-gene direction (2BP vs No-BP), top 12 by |diff| --")
ib=[i for i,m in enumerate(meta2) if m['grp']=='2BP']; iv=[i for i,m in enumerate(meta2) if m['grp']=='No-BP']
diffs=[]
for gid in BACT:
    if gid not in l2c: continue
    mb=sum(l2c[gid][i] for i in ib)/len(ib); mv=sum(l2c[gid][i] for i in iv)/len(iv)
    diffs.append((mb-mv,gid))
up=sum(1 for d,_ in diffs if d>0)
print(f"  bacterial genes higher in 2BP: {up}/{len(diffs)}")
for d,gid in sorted(diffs,key=lambda x:-abs(x[0]))[:12]: print(f"    {gid}  diff={d:+.3f}")
print("  -- viral program --")
dv=[]
for gid in VIR:
    if gid not in l2c: continue
    mb=sum(l2c[gid][i] for i in ib)/len(ib); mv=sum(l2c[gid][i] for i in iv)/len(iv)
    dv.append((mb-mv,gid))
print(f"  viral genes higher in 2BP: {sum(1 for d,_ in dv if d>0)}/{len(dv)}")
for d,gid in sorted(dv,key=lambda x:-abs(x[0]))[:8]: print(f"    {gid}  diff={d:+.3f}")

# patient identity (multi-sample patients?)
pat=defaultdict(list)
for m in meta2:
    hs=m['s'].split('-')[1] if '-' in m['s'] else m['s']
    pat[hs].append(m)
multi={k:len(v) for k,v in pat.items() if len(v)>1}
print(f"\n-- patients: {len(pat)} unique; with >1 sample: {len(multi)} {multi}")
