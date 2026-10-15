#!/usr/bin/env bash
# Fetch all decoy genomes (fna only) from the local assembly summaries.
set -u
cd /agent/workspace/abf
mkdir -p decoy
n=0; miss=0
while read -r acc; do
  [ -z "$acc" ] && continue
  out="decoy/$acc.fna.gz"
  if [ -s "$out" ]; then n=$((n+1)); continue; fi
  line=$(grep -m1 -P "^$acc\t" genomes/refseq_bac.txt genomes/genbank_bac.txt 2>/dev/null | head -1)
  if [ -z "$line" ]; then echo "MISS $acc"; miss=$((miss+1)); continue; fi
  ftp=$(echo "$line" | awk -F'\t' '{print $20}')
  stem=$(basename "$ftp")
  curl -sfL --retry 3 -o "$out" "$ftp/${stem}_genomic.fna.gz" || { echo "FAILDL $acc"; miss=$((miss+1)); continue; }
  n=$((n+1))
done < recon/decoy_accs.txt
echo "have=$n miss=$miss"
ls decoy/*.fna.gz 2>/dev/null | wc -l
du -sh decoy
