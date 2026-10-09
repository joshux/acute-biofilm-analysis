#!/bin/bash
# Process every downloaded sample in a run list through quantify_v3.py.
# usage: process_set.sh <runs.txt> <pseu|kleb> [min_bytes]
cd /agent/workspace/abf
LIST=$1; ORG=$2
while read -r S; do
  [ -z "$S" ] && continue
  F=fastq/$S.fastq.gz
  [ -s "$F" ] || { echo "MISSING $S"; continue; }
  [ -s "out/${S}_${ORG}_v3.json" ] && { echo "SKIP $S"; continue; }
  python3 quantify_v3.py "$S" "$ORG" "$F" 20 > /dev/null 2>&1
  [ -s "out/${S}_${ORG}_v3.json" ] && echo "OK $S" || echo "FAIL $S"
done < "$LIST"
