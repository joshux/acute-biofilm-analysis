#!/usr/bin/env python3
"""FINAL analysis.

Cohort is defined by the RNA reads themselves (rna_genus_mapped.json), not by the
deposited composition table, which disagrees for half the samples. A sample enters
the matrix-growth test only if its Pseudomonas fraction of mapped RNA reads is
>= 0.5 and it has >= 100 target reads in the quantification panel.

Statistic: B = log(R_matrix) - log(R_growth), R = length-normalised reads/kb.
  B < 0  matrix is a smaller share of the transcriptome than growth machinery
  B > 0  matrix is a larger share
Kolpen's claim (matrix-ON together with growth-FAST) predicts a POSITIVE
association between matrix and growth investment across samples.
"""
import json, csv, os, math, statistics as st, collections

ABF = "/tmp/abf"
os.chdir(ABF)

rna = json.load(open("rna_genus_mapped.json"))
q3 = {os.path.basename(f)[:-8].replace("_q3", ""): json.load(open(f))
      for f in os.listdir(".") if f.endswith("_q3.json")}
# fix key extraction: files are <SRR>_<org>_q3.json
q3 = {}
for f in os.listdir("."):
    if f.endswith("_q3.json"):
        d = json.load(open(f))
        q3[(d["srr"], d["org"])] = d

man = {r["srr"]: r for r in csv.DictReader(open("/tmp/bac_manifest.csv"))}


def spearman(x, y):
    def rank(a):
        s = sorted(range(len(a)), key=lambda i: a[i]); rk = [0] * len(a)
        for pos, i in enumerate(s): rk[i] = pos
        return rk
    rx, ry = rank(x), rank(y); n = len(x)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def sign_test(v):
    """two-sided sign test against 0"""
    from math import comb
    v = [x for x in v if x != 0]
    n = len(v); k = sum(1 for x in v if x > 0)
    p = sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n
    return min(1.0, 2 * p)


rows = []
for (srr, org), d in sorted(q3.items()):
    r = rna.get(srr, {"total": 0, "counts": {}})
    tot = r["total"]
    ps = r["counts"].get("Pseudomonas", 0)
    ac = r["counts"].get("Acinetobacter", 0)
    pf = ps / tot if tot else 0
    af = ac / tot if tot else 0
    m = man.get(srr, {})
    rows.append({
        "srr": srr, "org": org, "rna_total": tot, "pseud_frac": round(pf, 3),
        "aci_frac": round(af, 3), "rna_dom": max(r["counts"].items(), key=lambda x: x[1])[0] if r["counts"] else "?",
        "table_dom": m.get("dom", "?"), "micro": int(m.get("micro", 0)),
        "target": d["target_reads"], "decoy": d["decoy_hits"],
        "n_mat": d["reads"]["matrix"], "n_rp": d["reads"]["rp"], "n_og": d["reads"]["og"],
        "R_mat": d["R"]["matrix"], "R_rp": d["R"]["rp"], "R_og": d["R"]["og"], "R_grw": d["R"]["growth"],
        "f_mat": d["fractions"]["matrix"], "f_rp": d["fractions"]["rp"], "f_og": d["fractions"]["og"],
        "B_rp": d["B_rp"], "B_rp_lo": d["B_rp_lo"], "B_rp_hi": d["B_rp_hi"],
        "B_full": d["B_full"], "B_full_lo": d["B_full_lo"], "B_full_hi": d["B_full_hi"],
        "det_mat": d["detect"]["matrix"],
    })

rows.sort(key=lambda r: -r["pseud_frac"])
with open("results_final.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print("ALL samples with q3 (sorted by Pseudomonas RNA fraction)")
hdr = f"{'sample':13s}{'org':5s}{'rnaN':>7s}{'pseu%':>7s}{'aci%':>6s}{'target':>7s}{'mat':>5s}{'rp':>5s}{'og':>6s}{'B_rp':>7s}{'B_full':>8s}"
print(hdr); print("-" * len(hdr))
for r in rows:
    print(f"{r['srr']:13s}{r['org']:5s}{r['rna_total']:>7,}{r['pseud_frac']*100:>6.0f}%{r['aci_frac']*100:>5.0f}%"
          f"{r['target']:>7,}{r['n_mat']:>5}{r['n_rp']:>5}{r['n_og']:>6}{r['B_rp']:>7.2f}{r['B_full']:>8.2f}")

# --- cohort: Pseudomonas-dominant in the RNA, adequate depth ------------
coh = [r for r in rows if r["pseud_frac"] >= 0.5 and r["target"] >= 100 and r["org"] == "pseu"]
print(f"\n=== COHORT: Pseudomonas RNA fraction >= 0.5 and target >= 100  (n={len(coh)}) ===")
print(f"{'sample':13s}{'pseu%':>7s}{'target':>7s}{'mat':>5s}{'rp':>5s}{'og':>6s}{'B_rp':>7s}{'95% CI':>18s}{'B_full':>8s}{'f_mat':>7s}")
for r in coh:
    ci = f"[{r['B_rp_lo']:+.2f},{r['B_rp_hi']:+.2f}]"
    print(f"{r['srr']:13s}{r['pseud_frac']*100:>6.0f}%{r['target']:>7,}{r['n_mat']:>5}{r['n_rp']:>5}{r['n_og']:>6}"
          f"{r['B_rp']:>7.2f}{ci:>18s}{r['B_full']:>8.2f}{r['f_mat']:>7.3f}")

for key in ["B_rp", "B_full"]:
    v = [r[key] for r in coh]
    print(f"\n{key}: median {st.median(v):+.3f}  mean {st.mean(v):+.3f}  range [{min(v):+.2f},{max(v):+.2f}]"
          f"  n>0={sum(1 for x in v if x>0)}/{len(v)}  sign-test p={sign_test(v):.4f}")

x = [math.log(r["R_mat"]) for r in coh]; y = [math.log(r["R_grw"]) for r in coh]
print(f"\nCoupling: Spearman(log R_matrix, log R_growth) = {spearman(x, y):+.3f}")
z = [math.log(r["target"]) for r in coh]
print(f"Depth confound: Spearman(log target, B_full) = {spearman(z, [r['B_full'] for r in coh]):+.3f}")
print(f"Depth confound: Spearman(log target, B_rp)   = {spearman(z, [r['B_rp'] for r in coh]):+.3f}")
print(f"Pseudomonas fraction vs B_full: Spearman = {spearman([r['pseud_frac'] for r in coh], [r['B_full'] for r in coh]):+.3f}")

json.dump({"cohort": coh, "all": rows}, open("analysis_final.json", "w"), indent=1)
print("\nwrote results_final.csv / analysis_final.json")
