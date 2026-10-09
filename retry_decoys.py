import json, os, urllib.request, urllib.parse, subprocess, time
OUT="/tmp/abf/decoy"
FAILED=["Capnocytophaga ochracea","Comamonas testosteroni","Prevotella melaninogenica",
        "Moraxella catarrhalis","Brucella melitensis","Yersinia pestis","Enterococcus faecalis",
        "Agrobacterium tumefaciens","Rhizobium leguminosarum","Herbaspirillum seropedicae",
        "Streptococcus pneumoniae","Haemophilus influenzae","Escherichia coli","Acinetobacter baumannii"]
man=json.load(open(f"{OUT}/manifest.json"))
have={m["taxon"] for m in man}
def get_ref(taxon):
    q=urllib.parse.quote(taxon)
    for filt in ["filters.reference_only=true&","filters.assembly_level=complete&",""]:
        url=f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{q}/dataset_report?{filt}page_size=1"
        for attempt in range(4):
            try:
                r=json.load(urllib.request.urlopen(url,timeout=40))
                if r.get("reports"): 
                    rep=r["reports"][0]
                    return rep["accession"], rep["assembly_info"]["assembly_name"], rep["organism"]["organism_name"]
                break
            except Exception as e:
                time.sleep(6*(attempt+1))
    return None
for taxon in FAILED:
    if taxon in have: continue
    res=get_ref(taxon)
    if not res: print(f"STILL FAIL {taxon}",flush=True); continue
    acc,name,org=res
    parts=acc.split("_")[1]
    base=f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{parts[0:3]}/{parts[3:6]}/{parts[6:9]}/{acc}_{name}"
    url=f"{base}/{acc}_{name}_genomic.fna.gz"
    dest=f"{OUT}/{acc}.fna.gz"
    if not (os.path.exists(dest) and os.path.getsize(dest)>0):
        try: subprocess.run(["curl","-sfL","--retry","2","-o",dest,url],check=True,timeout=180)
        except Exception as e: print(f"DL FAIL {taxon}: {e}",flush=True); continue
    man.append({"taxon":taxon,"acc":acc,"org":org,"path":dest})
    print(f"OK {taxon:30s} {acc}",flush=True)
    time.sleep(3)
json.dump(man,open(f"{OUT}/manifest.json","w"),indent=1)
print(f"\ntotal {len(man)} genomes")
