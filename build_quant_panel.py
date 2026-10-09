#!/usr/bin/env python3
"""Build per-organism quantification panels with a single shared mapping run.

Goal: matrix vs growth must be measured in ONE minimap2 run over ONE reference,
so that both arms share identical mapping conditions (read length, scoring,
alignment filters). Previously matrix and growth were mapped separately, which
conflates panel composition with mapping behaviour.

Output header: >{cls}|{acc}|{gene}|{i}   cls in {matrix, rp, og}
  matrix = matrix / matrix_prod genes
  rp     = ribosomal-protein growth genes (classical growth marker)
  og     = other growth genes (replication / division / transcription / translation factors)
"""
import sys, os, json, collections

ABF = "/tmp/abf"
os.chdir(ABF)

ORGANISMS = {
    "pseu": ["GCF_000006765.1", "GCF_045571375.1", "GCF_047118485.1"],
    "kleb": ["GCA_021989295.1", "GCF_000240185.1", "GCF_988229555.1"],
}


def read_fasta(path):
    """yield (header, seq)"""
    h = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(buf)
                h = line[1:].rstrip("\n")
                buf = []
            else:
                buf.append(line.rstrip("\n"))
    if h is not None:
        yield h, "".join(buf)


# --- load panels ---------------------------------------------------------
mat = {}   # (acc,gene) -> seq  ; from panel_matrix.fna
for h, s in read_fasta("panel_matrix.fna"):
    cat, acc, gene = h.split("|")[:3]
    mat[(acc, gene)] = s

rp = set()  # (acc,gene) set of ribosomal-protein genes
for h, s in read_fasta("panel_rp.fna"):
    _, acc, gene = h.split("|")[:3]
    rp.add((acc, gene))

grw = {}   # (acc,gene) -> seq ; from panel_growth.fna
for h, s in read_fasta("panel_growth.fna"):
    cat, acc, gene = h.split("|")[:3]
    grw[(acc, gene)] = s

print(f"matrix genes: {len(mat)}")
print(f"growth genes: {len(grw)}   of which RP: {len(rp & set(grw))}")

summary = {}
for org, accs in ORGANISMS.items():
    out = []
    counts = collections.Counter()
    i = 0
    for acc in accs:
        # matrix
        for (a, g), s in sorted(mat.items()):
            if a != acc:
                continue
            out.append((f"matrix|{a}|{g}|{i}", s)); i += 1; counts["matrix"] += 1
        # growth split into rp / og
        for (a, g), s in sorted(grw.items()):
            if a != acc:
                continue
            cls = "rp" if (a, g) in rp else "og"
            out.append((f"{cls}|{a}|{g}|{i}", s)); i += 1; counts[cls] += 1
    fn = f"qp_{org}.fna"
    with open(fn, "w") as f:
        for h, s in out:
            f.write(">" + h + "\n" + s + "\n")
    # per-gene length stats for the null model
    lens = collections.defaultdict(list)
    for h, s in out:
        cls = h.split("|")[0]
        lens[cls].append(len(s))
    summary[org] = {"n": dict(counts), "file": fn,
                    "median_len": {k: sorted(v)[len(v)//2] for k, v in lens.items()}}
    print(f"{org}: {dict(counts)}  -> {fn}")

json.dump(summary, open("qp_panels.json", "w"), indent=1)
print("\nwrote qp_panels.json")
