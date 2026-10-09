#!/bin/bash
cd /agent/workspace/abf
while read -r s; do
  [ -z "$s" ] && continue
  [ -s "fastq/$s.fastq.gz" ] || continue
  [ -s "out/${s}_rp2.json" ] && continue
  python3 quantify_rp2.py "$s" >> /tmp/rp2run.log 2>&1
done < rp_all.txt
echo ALLDONE >> /tmp/rp2run.log
