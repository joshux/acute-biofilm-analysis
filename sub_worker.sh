#!/bin/bash
# Uniform-depth two-axis quantification for the decisive run list.
cd /tmp/abf
mkdir -p sub final_sub
while read S; do
  [ -z "$S" ] && continue
  [ -f "final_sub/${S}.json" ] && continue
  if [ ! -s "fastq/${S}_1.fastq.gz" ]; then
    bash dl_full.sh $S >> /tmp/dl2.log 2>&1
  fi
  bash sub500k.sh $S 500000 >> /tmp/sub.log 2>&1
  python3 quant_sub.py $S >> /tmp/quant_sub.log 2>&1
  echo "DONE $S" >> /tmp/quant_sub.log
done < /tmp/decide_srrs.txt
echo ALLDONE >> /tmp/quant_sub.log
