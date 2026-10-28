#!/usr/bin/env python3
"""recount_overlap.py — recount the rp2 growth axis with TRUE read-gene overlap.

Purpose (2026-10-28 implementation recount):
  quantify_rp2.py (and quantify_tf.py / quant_final.py / quant_any.py) classify a
  read as hitting a feature when its alignment START (SAM POS) falls within
  [feature_start-300, feature_end+300] — a +/-300 bp window that never looks at
  the aligned blocks (CIGAR). This script maps each sentinel sample ONCE and
  classifies every read under BOTH conventions on the SAME alignments:

    legacy   : quantify_rp2.py's exact binary search, POS only, +/-300 window
    overlap  : aligned reference blocks (CIGAR M/=/X; D/N advance ref) must
               intersect the feature interval by >= 1 aligned base

  Output: pctRP_coding under both, disagreement counts, and which genes the
  window-only "RP" reads actually land in (the spill diagnosis).

Usage: recount_overlap.py <SRR> <fastq.gz> <pao1.gff.gz> <mm2> <index.mmi> [out.json]

Self-tests: python3 recount_overlap.py --selftest
"""
import sys, os, json, gzip, re, collections

# ---------------------------------------------------------------- interval load
RP_RE = re.compile(r"^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$")

def load_features(gff_path):
    rna_iv = collections.defaultdict(list)   # seqid -> [(s,e)]
    rp_iv  = collections.defaultdict(list)   # seqid -> [(s,e)]
    cds_iv = collections.defaultdict(list)   # seqid -> [(s,e,gene)] for diagnostics
    n_rrna = n_rp = 0
    op = gzip.open if gff_path.endswith(".gz") else open
    with op(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            seqid, typ, s, e = c[0], c[2], int(c[3]), int(c[4])
            if typ in ("rRNA", "tRNA", "tmRNA"):
                rna_iv[seqid].append((s, e)); n_rrna += 1
            elif typ == "CDS":
                m = re.search(r"gene=([^;]+)", c[8]) or re.search(r"Name=([^;]+)", c[8])
                g = m.group(1) if m else ""
                cds_iv[seqid].append((s, e, g))
                if RP_RE.match(g):
                    rp_iv[seqid].append((s, e)); n_rp += 1
    for d in (rna_iv, rp_iv):
        for k in d:
            d[k].sort()
    for k in cds_iv:
        cds_iv[k].sort()
    return rna_iv, rp_iv, cds_iv, n_rrna, n_rp

# ---------------------------------------------------------------- CIGAR parsing
def cigar_blocks(pos, cigar):
    """Return (blocks, ref_span_end). blocks = [(rs, re_inclusive)] of ALIGNED
    ref bases (M/=/X). D/N advance reference but are not aligned bases."""
    blocks, rpos, cur_start, cur_len = [], pos, pos, 0
    for num, op_ in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        num = int(num)
        if op_ in "M=X":
            if cur_len == 0:
                cur_start = rpos
            cur_len += num
            rpos += num
        elif op_ in "DN":
            if cur_len:
                blocks.append((cur_start, cur_start + cur_len - 1)); cur_len = 0
            rpos += num
        # I, S consume query only; H, P consume neither
    if cur_len:
        blocks.append((cur_start, cur_start + cur_len - 1))
    return blocks, rpos - 1  # last ref base touched (incl D/N skips)

# ---------------------------------------------------------------- hit logic
def hit_window(ivs, pos, width=300):
    """quantify_rp2.py's EXACT logic (verbatim), POS-only +/-300."""
    if not ivs:
        return False
    lo, hi = 0, len(ivs)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid]
        if pos < s - 300:
            hi = mid
        elif pos > e + 300:
            lo = mid + 1
        else:
            return True
    return False

def hit_overlap(ivs, blocks):
    """True if any aligned block intersects any interval by >=1 base."""
    if not ivs or not blocks:
        return False
    for bs, be in blocks:
        lo, hi = 0, len(ivs)
        # first interval with e >= bs
        while lo < hi:
            mid = (lo + hi) // 2
            if ivs[mid][1] < bs:
                lo = mid + 1
            else:
                hi = mid
        i = lo
        while i < len(ivs) and ivs[i][0] <= be:
            return True   # intersecting pair found (first check suffices)
        # continue with next block
    return False

def containing_gene(cds_list, blocks):
    for bs, be in blocks:
        lo, hi = 0, len(cds_list)
        while lo < hi:
            mid = (lo + hi) // 2
            if cds_list[mid][1] < bs:
                lo = mid + 1
            else:
                hi = mid
        i = lo
        while i < len(cds_list) and cds_list[i][0] <= be:
            return cds_list[i][2]
    return None

# ---------------------------------------------------------------- self-tests
def selftest():
    b, end = cigar_blocks(100, "50M10D50M")
    assert b == [(100, 149), (160, 209)], b
    assert end == 209
    b, end = cigar_blocks(100, "20S80M")
    assert b == [(100, 179)] and end == 179, (b, end)
    b, end = cigar_blocks(100, "10I90M5H")
    assert b == [(100, 189)] and end == 189
    b, end = cigar_blocks(100, "30M5N40M")
    assert b == [(100, 129), (135, 174)] and end == 174
    ivs = [(1000, 1100), (2000, 2100)]
    assert hit_window(ivs, 700)          # 300 before gene start -> True (legacy)
    assert not hit_overlap(ivs, [(650, 749)])
    assert hit_overlap(ivs, [(1050, 1090)])
    assert hit_overlap(ivs, [(1095, 1150)])   # partial overlap counts
    assert not hit_overlap(ivs, [(1101, 1200)])
    # window counts a read landing entirely in the NEXT gene 250bp after rp end:
    assert hit_window(ivs, 2350) is True if False else hit_window(ivs, 2350)
    print("selftest OK")
    sys.exit(0)

if "--selftest" in sys.argv:
    selftest()

# ---------------------------------------------------------------- main
SRR   = sys.argv[1]
FASTQ = sys.argv[2]
GFF   = sys.argv[3]
MM2   = sys.argv[4]
INDEX = sys.argv[5]
OUT   = sys.argv[6] if len(sys.argv) > 6 else f"{SRR}_recount.json"
MAPQ  = 20

rna_iv, rp_iv, cds_iv, n_rrna, n_rp = load_features(GFF)

tot = rna_w = rp_w = rna_o = rp_o = decoy = 0
dis = collections.Counter()
spill_genes = collections.Counter()
rp_genes_w = collections.Counter()
rp_genes_o = collections.Counter()

# gene names for rp intervals (for distinct-gene counts)
rp_named = collections.defaultdict(list)
for seqid, ivs in rp_iv.items():
    for (s, e) in ivs:
        for (cs, ce, g) in cds_iv.get(seqid, []):
            if cs == s and ce == e:
                rp_named[seqid].append((s, e, g))
                break

cmd = f"{MM2} -x sr -t 2 -a --secondary=no {INDEX} {FASTQ} 2>/dev/null"
with os.popen(cmd) as fh:
    for line in fh:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        if len(f) < 7 or int(f[1]) & 4 or int(f[4]) < MAPQ:
            continue
        ref = f[2]
        if ref.startswith("decoy|"):
            decoy += 1
            continue
        tot += 1
        pos = int(f[3])
        cigar = f[5] if len(f) > 5 else "*"
        blocks, span_end = (cigar_blocks(pos, cigar) if cigar != "*" else ([], pos))

        # legacy (verbatim rp2 logic)
        w_rna = hit_window(rna_iv.get(ref, []), pos)
        w_rp  = (not w_rna) and hit_window(rp_iv.get(ref, []), pos)
        # true overlap
        o_rna = hit_overlap(rna_iv.get(ref, []), blocks)
        o_rp  = (not o_rna) and hit_overlap(rp_iv.get(ref, []), blocks)

        rna_w += w_rna; rp_w += w_rp
        rna_o += o_rna; rp_o += o_rp

        if w_rna and not o_rna: dis["rna_window_only"] += 1
        if o_rna and not w_rna: dis["rna_overlap_only"] += 1
        if w_rp and not o_rp:
            dis["rp_window_only"] += 1
            g = containing_gene(cds_iv.get(ref, []), blocks)
            spill_genes[g or "intergenic/no-CDS"] += 1
        if o_rp and not w_rp: dis["rp_overlap_only"] += 1
        if w_rp:
            for (cs, ce, g) in rp_named.get(ref, []):
                if hit_window([(cs, ce)], pos):
                    rp_genes_w[g] += 1
                    break
        if o_rp:
            for (cs, ce, g) in rp_named.get(ref, []):
                if hit_overlap([(cs, ce)], blocks):
                    rp_genes_o[g] += 1
                    break

den_w = tot - rna_w
den_o = tot - rna_o
out = {
    "srr": SRR, "mapq": MAPQ, "n_rna_feat": n_rrna, "n_rp_genes": n_rp,
    "tot_pao1": tot, "decoy_reads": decoy,
    "legacy": {"rrna": rna_w, "coding": den_w, "rp": rp_w,
               "pctRP_coding": round(100 * rp_w / den_w, 2) if den_w else None},
    "overlap": {"rrna": rna_o, "coding": den_o, "rp": rp_o,
                "pctRP_coding": round(100 * rp_o / den_o, 2) if den_o else None},
    "delta_pctRP": round(100 * rp_o / den_o - 100 * rp_w / den_w, 2) if (den_o and den_w) else None,
    "disagreements": dict(dis),
    "rp_window_only_actually_in": dict(spill_genes.most_common(12)),
    "distinct_rp_genes_window": len(rp_genes_w),
    "distinct_rp_genes_overlap": len(rp_genes_o),
}
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps(out))
