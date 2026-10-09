#!/usr/bin/env python3
"""Strain-robustness check for the growth axis.

If the acute Pseudomonas is not PAO1-like, conserved RP reads could map while
non-RP reads do not, inflating %RP for acute only. Test: recompute %RP against
three different single-genome references (PAO1, PA14, P. putida). If the
acute > chronic contrast holds with all three, strain divergence is not driving it.
"""
import sys, os, json, collections, gzip, re

ABF = "/agent/workspace/abf"
os.chdir(ABF)
MM2 = "bin/minimap2-2.28_x64-linux/minimap2"
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")

REFS = {
    "PAO1":   ("GCF_000006765.1", "refs2/GCF_000006765.1.fna.gz"),
    "PA14":   ("GCF_047118485.1", "refs2/GCF_047118485.1.fna.gz"),
    "putida": ("GCF_045571375.1", "refs2/GCF_045571375.1.fna.gz"),
}
panels = json.load(open("gff/panels.json"))


def rp_intervals(acc):
    iv = collections.defaultdict(list)
    for g in panels.get(acc, []):
        if g["cat"].startswith("growth") and RP_RE.match(g["gene"]):
            iv[g["seqid"]].append((g["start"], g["end"]))
    return iv


def build_ref(name):
    acc, gz = REFS[name]
    mmi = f"compR_{name}.mmi"
    if os.path.exists(mmi):
        return mmi
    with open(f"compR_{name}.fna", "w") as out:
        with gzip.open(gz, "rt") as f:
            h = None
            buf = []
            for line in f:
                if line.startswith(">"):
                    if h:
                        out.write(f">{h}\n{''.join(buf)}\n")
                    h = line[1:].split()[0]; buf = []
                else:
                    buf.append(line.strip())
            if h:
                out.write(f">{h}\n{''.join(buf)}\n")
        for d in json.load(open("decoy/manifest.json")):
            p = f"decoy/{d['acc']}.fna.gz"
            if not os.path.exists(p):
                continue
            with gzip.open(p, "rt") as f:
                h = None
                buf = []
                for line in f:
                    if line.startswith(">"):
                        if h:
                            out.write(f">decoy|{h}\n{''.join(buf)}\n")
                        h = line[1:].split()[0]; buf = []
                    else:
                        buf.append(line.strip())
                if h:
                    out.write(f">decoy|{h}\n{''.join(buf)}\n")
    os.system(f"{MM2} -d {mmi} compR_{name}.fna 2>/dev/null")
    return mmi


def pct_rp(name, srr):
    acc, _ = REFS[name]
    mmi = f"compR_{name}.mmi"
    iv = rp_intervals(acc)
    tot = rp = 0
    cmd = f"{MM2} -x sr -t 2 -a --secondary=no {mmi} fastq/{srr}.fastq.gz 2>/dev/null"
    with os.popen(cmd) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.split("\t", 6)
            if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < 20:
                continue
            if f[2].startswith("decoy|"):
                continue
            tot += 1
            pos = int(f[3])
            ivs = iv.get(f[2], [])
            lo, hi = 0, len(ivs)
            while lo < hi:
                mid = (lo + hi) // 2
                s, e = ivs[mid]
                if pos < s - 300:
                    hi = mid
                elif pos > e + 300:
                    lo = mid + 1
                else:
                    rp += 1; break
    return tot, rp, (round(100 * rp / tot, 2) if tot else None)


if __name__ == "__main__":
    samples = sys.argv[1:]
    out = {}
    for srr in samples:
        row = {}
        for name in REFS:
            build_ref(name)
            tot, rp, pc = pct_rp(name, srr)
            row[name] = {"total": tot, "rp": rp, "pctRP": pc}
        out[srr] = row
        print(f"{srr}: " + "  ".join(f"{k}={v['pctRP']}%" for k, v in row.items()), flush=True)
    json.dump(out, open("strain_robust.json", "w"), indent=1)
