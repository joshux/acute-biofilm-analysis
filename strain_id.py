import os,sys,collections,json
os.chdir("/agent/workspace/abf")
MM2="bin/minimap2-2.28_x64-linux/minimap2"
NAMES={"GCF_026167805.1":"A.baumannii_17978","GCF_000241685.1":"A.baumannii_AB5075",
       "GCF_057224575.1":"A.nosocomialis","GCF_055475515.1":"A.pittii",
       "GCF_988236795.1":"A.johnsonii","GCF_050703235.1":"A.lwoffii"}
for s in sys.argv[1:]:
    fq=f"fastq/{s}.fastq.gz"
    if not os.path.exists(fq): continue
    c=collections.Counter()
    cmd=f"{MM2} -x sr -t 2 -a --secondary=no acine_multi.mmi {fq} 2>/dev/null"
    with os.popen(cmd) as fh:
        for line in fh:
            if line[0]=="@": continue
            f=line.split("\t",6)
            if len(f)<7 or int(f[1])&4: continue
            acc=f[2].split("|")[0]
            if acc.startswith("decoy"): c["<decoy>"]+=1; continue
            c[NAMES.get(acc,acc)]+=1
    tot=sum(v for k,v in c.items() if k!="<decoy>")
    top="  ".join(f"{k}:{v}" for k,v in c.most_common(4) if k!="<decoy>")
    print(f"{s} acine_reads={tot:6d}  {top}")
