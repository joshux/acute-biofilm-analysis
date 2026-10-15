#!/usr/bin/env bash
# Fetch the remaining target-organism genomes + GFFs.
set -u
cd /agent/workspace/abf
mkdir -p genomes
ACCS="GCF_000013425.1 GCF_026167805.1 GCF_000012185.1 GCF_000006885.1 GCF_000072485.1"
for acc in $ACCS; do
  line=$(grep -m1 -P "^$acc\t" genomes/refseq_bac.txt genomes/genbank_bac.txt 2>/dev/null | head -1)
  if [ -z "$line" ]; then echo "MISSING $acc"; continue; fi
  ftp=$(echo "$line" | awk -F'\t' '{print $20}')
  org=$(echo "$line" | awk -F'\t' '{print $8}')
  stem=$(basename "$ftp")
  for pair in "_genomic.fna.gz:fna.gz" "_genomic.gff.gz:gff.gz"; do
    suf="${pair%%:*}"; ext="${pair##*:}"
    out="genomes/$acc.$ext"
    [ -s "$out" ] || curl -sfL --retry 3 -o "$out" "$ftp/$stem$suf"
  done
  sz1=$(stat -c%s "genomes/$acc.fna.gz" 2>/dev/null || echo 0)
  sz2=$(stat -c%s "genomes/$acc.gff.gz" 2>/dev/null || echo 0)
  printf "%-18s %-34s fna=%10s gff=%10s\n" "$acc" "${org:0:34}" "$sz1" "$sz2"
done
echo DONE
