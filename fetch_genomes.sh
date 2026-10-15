#!/usr/bin/env bash
# Fetch genomes + GFFs by streaming the assembly_summary (no big in-memory dict).
set -u
cd /agent/workspace/abf
mkdir -p genomes
BASE="https://ftp.ncbi.nlm.nih.gov/genomes"

fetch_summary () {
  local name="$1" url="$2"
  if [ ! -s "genomes/$name" ]; then
    echo "downloading $name ..."
    curl -sfL --retry 3 -o "genomes/$name" "$url"
  fi
  echo "$name bytes: $(stat -c%s genomes/$name 2>/dev/null || echo 0)"
}

fetch_summary refseq_bac.txt "$BASE/refseq/bacteria/assembly_summary.txt"
fetch_summary genbank_bac.txt "$BASE/genbank/bacteria/assembly_summary.txt"

# accessions we need
ACCS="GCF_000240185.1 GCA_021989295.1 GCF_988229555.1 GCF_000016305.1
GCF_000006765.1 GCF_047118485.1 GCF_000005845.2
GCF_000006945.2 GCF_000018665.1 GCF_016027095.1 GCF_000696675.2 GCF_009586235.1 GCF_000009085.1"

for acc in $ACCS; do
  # find the line in either summary (strip leading # on header rows)
  line=$(grep -m1 -P "^$acc\t" genomes/refseq_bac.txt genomes/genbank_bac.txt 2>/dev/null | head -1)
  if [ -z "$line" ]; then
    echo "MISSING $acc"
    continue
  fi
  # ftp_path is column 20
  ftp=$(echo "$line" | awk -F'\t' '{print $20}')
  org=$(echo "$line" | awk -F'\t' '{print $8}')
  if [ -z "$ftp" ] || [ "$ftp" = "na" ]; then echo "no ftp_path $acc"; continue; fi
  stem=$(basename "$ftp")
  for pair in "_genomic.fna.gz:fna.gz" "_genomic.gff.gz:gff.gz"; do
    suf="${pair%%:*}"; ext="${pair##*:}"
    out="genomes/$acc.$ext"
    if [ ! -s "$out" ]; then
      curl -sfL --retry 3 -o "$out" "$ftp/$stem$suf"
    fi
  done
  sz1=$(stat -c%s "genomes/$acc.fna.gz" 2>/dev/null || echo 0)
  sz2=$(stat -c%s "genomes/$acc.gff.gz" 2>/dev/null || echo 0)
  printf "%-18s %-10s fna=%12s gff=%12s\n" "$acc" "${org:0:34}" "$sz1" "$sz2"
done
echo "DONE"
ls -la genomes/*.gz | wc -l
