#!/usr/bin/env python3
"""
Host-side viral->bacterial phase test using the Mick/Tsitsiklis 108-gene signature.

Score = mean( z(log2 CPM) over 83 bacterial-high genes )
      - mean( z(log2 CPM) over 25 viral-high genes )

Baseline: GSE246795 (Spottiswoode/Langelier 2BP-labelled COVID ETA) -- reproduce the committed null.
Test:     GSE161731 (McClain/Duke whole blood: COVID-19 + bacterial pneumonia + influenza + other CoV + healthy)
          - transferability: does the signature separate ADULT bacterial pneumonia from viral?
          - phase: does the bacterial program rise (and viral fall) across COVID time-since-onset?
"""
import json, gzip, csv, math, random, sys
from collections import defaultdict

random.seed(20261012)

# ---------- signature ----------
sig = json.load(open('mick108_raw.json'))
BACT = [g['ens'] for g in sig if g['logFC'] < 0]   # 83 bacterial-high (down in purely-viral)
VIR  = [g['ens'] for g in sig if g['logFC'] > 0]   # 25 viral-high
print(f"signature: {len(BACT)} bacterial-high, {len(VIR)} viral-high")

# ---------- helpers ----------
def read_counts(path):
    """Return (genes, samples, matrix dict gene->list) from a GEO counts csv (first col = gene id)."""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        r = csv.reader(f)
        header = next(r)
        samples = header[1:]
        genes = {}
        for row in r:
            if not row: continue
            g = row[0].strip().strip('"')
            try:
                genes[g] = [float(x) for x in row[1:]]
            except ValueError:
                continue
    return samples, genes

def zscore_log2cpm(genes, samples, gene_set):
    """mean z(log2 CPM) over gene_set. z computed per gene across samples."""
    cols = len(samples)
    # library sizes
    lib = [0.0]*cols
    for g, vals in genes.items():
        for i, v in enumerate(vals):
            lib[i] += v
    present = [g for g in gene_set if g in genes]
    acc = [0.0]*cols
    for g in present:
        vals = genes[g]
        y = [math.log2((v/lib[i]*1e6) + 1.0) for i, v in enumerate(vals)]
        m = sum(y)/cols
        sd = math.sqrt(sum((x-m)**2 for x in y)/cols) if cols > 1 else 0.0
        if sd == 0: continue
        for i in range(cols):
            acc[i] += (y[i]-m)/sd
    n = len(present) if present else 1
    return [a/n for a in acc], present

def mwu_p(a, b):
    """two-sided Mann-Whitney U, normal approx with tie correction."""
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0: return float('nan'), float('nan')
    allv = sorted(a+b)
    # ranks
    rank = {}
    i = 0
    while i < len(allv):
        j = i
        while j+1 < len(allv) and allv[j+1] == allv[i]: j += 1
        r = (i+j)/2.0 + 1
        rank[allv[i]] = r
        i = j+1
    R1 = sum(rank[x] for x in a)
    U1 = R1 - n1*(n1+1)/2.0
    mu = n1*n2/2.0
    # tie-corrected sigma
    _, cnt = zip(*[(v, allv.count(v)) for v in set(allv)]) if allv else (None, None)
    ties = sum(c**3 - c for c in cnt)
    N = n1+n2
    sigma = math.sqrt(n1*n2/12.0 * ((N+1) - ties/(N*(N-1)))) if N > 1 else 0
    z = (U1-mu)/sigma if sigma > 0 else 0.0
    p = math.erfc(abs(z)/math.sqrt(2))
    return p, U1

def perm_p(a, b, n=20000):
    """two-sided permutation p on difference of medians."""
    obs = median(a) - median(b)
    pool = a+b
    k = len(a)
    cnt = 0
    for _ in range(n):
        random.shuffle(pool)
        d = median(pool[:k]) - median(pool[k:])
        if abs(d) >= abs(obs): cnt += 1
    return (cnt+1)/(n+1), obs

def median(x):
    s = sorted(x); n = len(s)
    return s[n//2] if n % 2 else (s[n//2-1]+s[n//2])/2.0

def boot_ci(a, n=5000):
    s=[]
    for _ in range(n):
        s.append(median([random.choice(a) for _ in a]))
    s.sort()
    return s[int(0.025*n)], s[int(0.975*n)]

def report(name, groups, order):
    print(f"\n### {name}")
    for k in order:
        v = groups.get(k, [])
        if v:
            lo, hi = boot_ci(v)
            print(f"  {k:12s} n={len(v):3d}  median={median(v):+.3f}  [95% CI {lo:+.3f},{hi:+.3f}]  mean={sum(v)/len(v):+.3f}")

# ============================================================
# PART 1 -- baseline: GSE246795 (2BP-labelled COVID)
# ============================================================
print("\n" + "="*70)
print("PART 1  BASELINE: GSE246795 (2BP-labelled COVID tracheal aspirate)")
print("="*70)

# labels from series matrix
lab = {}
with gzip.open('GSE246795_series_matrix.txt.gz','rt') as f:
    titles=None
    for line in f:
        if line.startswith('!Sample_title'):
            titles=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        if line.startswith('!Sample_characteristics_ch1') and '2bp_status' in line:
            vals=[x.strip().strip('"').replace('2bp_status:','').strip() for x in line.rstrip('\n').split('\t')[1:]]
            lab=dict(zip(titles,vals))
            break

samples795, genes795 = read_counts('GSE246795_host_counts.csv.gz')
sc_b, pres_b = zscore_log2cpm(genes795, samples795, BACT)
sc_v, pres_v = zscore_log2cpm(genes795, samples795, VIR)
score795 = [b - v for b, v in zip(sc_b, sc_v)]
print(f"genes present: bacterial {len(pres_b)}/83, viral {len(pres_v)}/25 ; samples {len(samples795)}")

g795 = defaultdict(list)
for s, sc in zip(samples795, score795):
    g795[lab.get(s, '?')].append(sc)
report('GSE246795  bacterial-minus-viral score by 2BP status', g795, ['2BP','No-BP'])
p, U = mwu_p(g795['2BP'], g795['No-BP'])
pp, d = perm_p(g795['2BP'], g795['No-BP'])
print(f"  MWU p={p:.3f}  permutation p={pp:.3f}  median diff={d:+.3f}")

# ============================================================
# PART 2 -- GSE161731 (adult whole blood, multi-cohort)
# ============================================================
print("\n" + "="*70)
print("PART 2  GSE161731 (adult whole blood: COVID / bacterial / influenza / other-CoV / healthy)")
print("="*70)

key={}
with gzip.open('GSE161731_counts_key.csv.gz','rt') as f:
    for row in csv.DictReader(f):
        key[row['rna_id']] = row

samples, genes = read_counts('GSE161731_counts.csv.gz')
sc_b, pres_b = zscore_log2cpm(genes, samples, BACT)
sc_v, pres_v = zscore_log2cpm(genes, samples, VIR)
print(f"genes present: bacterial {len(pres_b)}/83, viral {len(pres_v)}/25 ; samples {len(samples)}")

rows=[]
for s, b, v in zip(samples, sc_b, sc_v):
    k = key.get(s)
    if not k: continue
    rows.append({'rna':s,'subj':k['subject_id'],'cohort':k['cohort'],
                 'time':k['time_since_onset'],'hosp':k['hospitalized'],
                 'bact':b,'vir':v,'score':b-v})
print(f"samples joined to metadata: {len(rows)}")
cohorts = defaultdict(int)
for r in rows: cohorts[r['cohort']] += 1
print("cohorts:", dict(cohorts))

# ---- 2a transferability: bacterial vs viral infection in adults ----
g = defaultdict(list)
for r in rows: g[r['cohort']].append(r['score'])
report('score by cohort (adult whole blood)', g,
       ['Bacterial','Influenza','CoV other','COVID-19','Healthy'])

for a,b in [('Bacterial','Influenza'),('Bacterial','COVID-19'),('Bacterial','Healthy'),
            ('COVID-19','Healthy'),('Influenza','Healthy')]:
    if a in g and b in g:
        p,U = mwu_p(g[a],g[b]); pp,d = perm_p(g[a],g[b])
        print(f"  {a:10s} vs {b:10s}: median diff={d:+.3f}  MWU p={p:.4f}  perm p={pp:.4f}")

# component check
print("\n  -- component means (bacterial-program vs viral-program separately) --")
for comp in ['bact','vir']:
    gg=defaultdict(list)
    for r in rows: gg[r['cohort']].append(r[comp])
    line="  "+comp+": "
    for c in ['Bacterial','Influenza','COVID-19','Healthy']:
        if c in gg: line+=f"{c}={sum(gg[c])/len(gg[c]):+.2f}(n{len(gg[c])})  "
    print(line)

# ---- 2b longitudinal: COVID time-since-onset ----
print("\n  -- COVID-19 by time since onset --")
gt=defaultdict(list)
for r in rows:
    if r['cohort']=='COVID-19' and r['time'] not in (None,'NA',''):
        gt[r['time']].append(r)
report('COVID score by time-since-onset', {k:[x['score'] for x in v] for k,v in gt.items()},
       ['early','middle','late'])
for comp in ['bact','vir']:
    print(f"    {comp}:", {k:round(sum(x[comp] for x in v)/len(v),3) for k,v in gt.items()})

# within-patient paired (early -> late)
by_subj=defaultdict(dict)
for r in rows:
    if r['cohort']=='COVID-19' and r['time'] in ('early','middle','late'):
        by_subj[r['subj']].setdefault(r['time'],[]).append(r)
pairs=[]
for s,d in by_subj.items():
    if 'early' in d and 'late' in d:
        e=sum(x['score'] for x in d['early'])/len(d['early'])
        l=sum(x['score'] for x in d['late'])/len(d['late'])
        pairs.append((s,e,l,l-e))
print(f"\n  within-patient early->late pairs: n={len(pairs)}")
if pairs:
    d=[x[3] for x in pairs]
    print(f"    median within-patient change={median(d):+.3f}  mean={sum(d)/len(d):+.3f}")
    pos=sum(1 for x in d if x>0)
    # sign test
    n=len(d); k=max(pos,n-pos)
    from math import comb
    pv=min(1.0, 2*sum(comb(n,i) for i in range(k,n+1))/2**n)
    print(f"    patients with increase: {pos}/{n}  sign-test p={pv:.3f}")

# ---- 2c severity proxy: hospitalisation ----
print("\n  -- by hospitalisation (severity proxy) --")
gh=defaultdict(list)
for r in rows:
    if r['hosp'] not in (None,'NA',''):
        gh[r['hosp']].append(r['score'])
report('score by hospitalised', gh, ['Yes','No'])
if 'Yes' in gh and 'No' in gh:
    p,U=mwu_p(gh['Yes'],gh['No']); print(f"    MWU p={p:.4f}")

json.dump({'rows':rows,
           'gse246795':{k:v for k,v in g795.items()},
           'pairs':pairs}, open('host_phase_results.json','w'), indent=1)
print("\nsaved host_phase_results.json")
