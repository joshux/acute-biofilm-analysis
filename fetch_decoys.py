#!/usr/bin/env python3
"""Fetch one reference genome per co-infecting genus to act as mapping decoys.
Conserved genes (fusA, rpoB, tuf...) map near-identically across genera; without
a competitor in the index they collapse onto the target organism and inflate
whatever class they belong to. Adding these genomes lets minimap2's MAPQ
correctly penalise cross-genus hits."""
import json, os, urllib.request, urllib.parse, gzip, shutil, subprocess, sys

OUT = "/tmp/abf/decoy"
os.makedirs(OUT, exist_ok=True)

# genera to add as decoys (target genera Pseudomonas/Klebsiella excluded -- already in panel)
GENERA = [
    "Salmonella enterica", "Delftia acidovorans", "Sphingomonas paucimobilis",
    "Rhodococcus erythropolis", "Xanthomonas campestris", "Burkholderia cepacia",
    "Campylobacter jejuni", "Capnocytophaga ochracea", "Comamonas testosteroni",
    "Enterobacter cloacae", "Prevotella melaninogenica", "Mycoplasma pneumoniae",
    "Ralstonia pickettii", "Porphyromonas gingivalis", "Fusobacterium nucleatum",
    "Treponema denticola", "Neisseria meningitidis", "Bacteroides fragilis",
    "Tannerella forsythia", "Streptococcus mutans", "Veillonella parvula",
    "Cutibacterium acnes", "Moraxella catarrhalis", "Brucella melitensis",
    "Yersinia pestis", "Enterococcus faecalis", "Agrobacterium tumefaciens",
    "Rhizobium leguminosarum", "Herbaspirillum seropedicae", "Alcaligenes faecalis",
]


def get_ref(taxon):
    q = urllib.parse.quote(taxon)
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{q}"
           f"/dataset_report?filters.reference_only=true&page_size=1")
    r = json.load(urllib.request.urlopen(url, timeout=40))
    reps = r.get("reports", [])
    if not reps:
        # fall back to any complete assembly
        url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{q}"
               f"/dataset_report?filters.assembly_level=complete&page_size=1")
        r = json.load(urllib.request.urlopen(url, timeout=40))
        reps = r.get("reports", [])
    if not reps:
        return None
    rep = reps[0]
    return rep["accession"], rep["assembly_info"]["assembly_name"], rep["organism"]["organism_name"]


manifest = []
for taxon in GENERA:
    try:
        res = get_ref(taxon)
    except Exception as e:
        print(f"API FAIL {taxon}: {e}", flush=True); continue
    if not res:
        print(f"NO GENOME {taxon}", flush=True); continue
    acc, name, org = res
    # build FTP path
    parts = acc.split("_")[1]
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{parts[0:3]}/{parts[3:6]}/{parts[6:9]}/{acc}_{name}"
    url = f"{base}/{acc}_{name}_genomic.fna.gz"
    dest = f"{OUT}/{acc}.fna.gz"
    if not os.path.exists(dest) or os.path.getsize(dest) == 0:
        try:
            subprocess.run(["curl", "-sfL", "--retry", "2", "-o", dest, url], check=True, timeout=180)
        except Exception as e:
            print(f"DL FAIL {taxon} {url}: {e}", flush=True); continue
    manifest.append({"taxon": taxon, "acc": acc, "org": org, "path": dest})
    print(f"OK {taxon:32s} {acc:18s} {os.path.getsize(dest):>9,}b", flush=True)

json.dump(manifest, open(f"{OUT}/manifest.json", "w"), indent=1)
print(f"\n{len(manifest)}/{len(GENERA)} genomes fetched -> {OUT}/manifest.json")
