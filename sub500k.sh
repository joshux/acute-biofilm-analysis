#!/bin/bash
# Subsample both mates to NREADS reads each into sub/ ; usage: sub500k.sh SRR [nreads]
cd /tmp/abf
S=$1
N=${2:-500000}
mkdir -p sub
for M in 1 2; do
  SRC=fastq/${S}_${M}.fastq.gz
  DST=sub/${S}_${M}.fastq.gz
  [ -s "$DST" ] && continue
  [ -s "$SRC" ] || continue
  zcat $SRC 2>/dev/null | head -n $((4*N)) | gzip -c > $DST
done
ls -la sub/${S}_1.fastq.gz sub/${S}_2.fastq.gz 2>/dev/null | awk '{print $NF, $5}'
