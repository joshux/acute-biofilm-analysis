#!/usr/bin/env python3
"""Derive genus composition directly from RNA reads (panel_multi index) and
compare against the deposited composition table. Defines the cohort by what the
RNA actually contains, not by the table."""
import os, glob, json, collections, csv

ABF = "/tmp/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
ref2g = json.load(open("ref2genus.json"))

out = {}
for fq in sorted(glob.glob("fastq/*.fastq.gz")):
    srr = os.path.basename(fq)[:-9]
    cmd = f"{MM2} -x sr -t 4 -a --secondary=no panel_multi.mmi {fq} 2>/dev/null"
    c = collections.Counter()
    tot = 0
    with os.popen(cmd) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.split("\t")
            if len(f) < 11 or int(f[1]) & 4:
                continue
            if int(f[4]) < 20:
                continue
            c[ref2g.get(f[2].split()[0], "?")] += 1
            tot += 1
    out[srr] = {"total": tot, "counts": dict(c)}
    top = c.most_common(3)
    tops = "  ".join(f"{g}:{n}" for g, n in top)
    print(f"{srr}  n={tot:6d}  {tops}", flush=True)

json.dump(out, open("rna_genus_mapped.json", "w"), indent=1)

# compare to table
comp = {}
with open("/tmp/rna_genus.txt") as f:
    hdr = f.readline().rstrip("\n").split("\t")[1:]
    for line in f:
        p = line.rstrip("\n").split("\t")
        for pid, v in zip(hdr, p[1:]):
            try: comp.setdefault(pid, {})[p[0]] = float(v)
            except: pass
rows = list(csv.DictReader(open("/tmp/prjna1056765_full_mapping.csv")))
srr2pid = {r["sra_run_accession"]: r["patient_id"] for r in rows}

print("\n=== mapped-dominant vs table-dominant ===")
agree = 0; tot = 0
for srr, d in out.items():
    if d["total"] < 50:
        continue
    tot += 1
    md = max(d["counts"].items(), key=lambda x: x[1])[0]
    c = comp.get(srr2pid.get(srr, ""), {})
    td = max(c.items(), key=lambda x: x[1])[0] if c else "?"
    ok = (md == td)
    agree += ok
    print(f"  {srr}  mapped={md:16s} ({max(d['counts'].values())/d['total']*100:4.0f}%)   "
          f"table={td:16s} ({max(c.values())*100:4.0f}%)   {'AGREE' if ok else 'DIFFER'}")
print(f"\nagreement: {agree}/{tot}")
