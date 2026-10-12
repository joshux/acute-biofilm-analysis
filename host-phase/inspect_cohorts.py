#!/usr/bin/env python3
"""Inspect arm structure + column mapping for the three replication cohorts."""
import gzip, re
from collections import Counter, defaultdict

def series(path):
    sm={'titles':[],'gsm':[],'chars':[],'src':[]}
    with gzip.open(path,'rt',errors='replace') as f:
        for line in f:
            if line.startswith('!Sample_title'): sm['titles']=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_geo_accession'): sm['gsm']=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'): sm['chars'].append([x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]])
            elif line.startswith('!Sample_source_name_ch1'): sm['src']=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
    return sm

for name,path in [('GSE157103','GSE157103_series.txt.gz'),
                  ('GSE152641','GSE152641_series.txt.gz'),
                  ('GSE243217','GSE243217_series.txt.gz')]:
    sm=series(path)
    n=len(sm['titles'])
    print(f"\n{'='*70}\n{name}: {n} samples")
    print("  first titles:", sm['titles'][:3])
    print("  last titles :", sm['titles'][-3:])
    # find disease/state field
    for ch in sm['chars']:
        keys=set(x.split(':')[0].strip().lower() for x in ch if ':' in x)
        vals=Counter(x.split(':',1)[1].strip() for x in ch if ':' in x)
        if len(vals)<=12:
            print(f"  field {list(keys)[0]!r:26s}: {dict(vals)}")
        else:
            print(f"  field {list(keys)[0]!r:26s}: {len(vals)} distinct (sample: {list(vals)[:3]})")
