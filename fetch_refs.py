import urllib.request, urllib.parse, json, os, time
QUERIES = {
 "acinetobacter": ["Acinetobacter baumannii ATCC 17978","Acinetobacter nosocomialis","Acinetobacter pittii",
                   "Acinetobacter johnsonii","Acinetobacter lwoffii","Acinetobacter baumannii AB5075"],
 "pseudomonas":   ["Pseudomonas aeruginosa PA14","Pseudomonas putida KT2440"],
 "stenotrophomonas":["Stenotrophomonas maltophilia K279a"],
 "klebsiella":    ["Klebsiella pneumoniae NTUH-K2044","Klebsiella pneumoniae HS11286","Klebsiella variicola"],
 "staph":         ["Staphylococcus aureus USA300 FPR3757"],
 "ecoli":         ["Escherichia coli O157:H7 EDL933"],
}
def get(u, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(u, timeout=60).read()
        except Exception as e:
            if '429' in str(e) and i<tries-1:
                time.sleep(3*(i+1)); continue
            raise
def esearch(q):
    return json.loads(get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={urllib.parse.quote(q)}&retmode=json"))['esearchresult'].get('idlist',[])
def esummary(i):
    return json.loads(get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={i}&retmode=json"))['result'][i]
os.makedirs("refs2",exist_ok=True)
manifest=[]
for genus,terms in QUERIES.items():
    for t in terms:
        try:
            ids=esearch(t); time.sleep(0.7)
            if not ids: print(f"  NO HIT: {t}"); continue
            d=esummary(ids[0]); time.sleep(0.7)
            ftp=d.get('ftppath_refseq') or d.get('ftppath_genbank')
            acc=d.get('assemblyaccession'); org=d.get('organism'); name=d.get('assemblyname')
            if not ftp: print(f"  no ftp: {t}"); continue
            base=ftp.replace("ftp://","https://").rstrip('/')
            fname=f"{base}/{base.split('/')[-1]}_genomic.fna.gz"
            out=f"refs2/{acc}.fna.gz"
            if not os.path.exists(out):
                with open(out,'wb') as fh: fh.write(get(fname))
            sz=os.path.getsize(out)
            manifest.append(dict(genus=genus,acc=acc,org=org,name=name,path=out,size=sz))
            print(f"  OK {genus:16} {acc:16} {org[:42]:42} {sz/1e6:.2f}MB")
        except Exception as e:
            print(f"  FAIL {t}: {e}")
        time.sleep(0.5)
json.dump(manifest,open("refs2/manifest.json","w"),indent=1)
print(f"\ntotal new genomes: {len(manifest)}")
