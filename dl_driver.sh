#!/bin/bash
# Detached downloader: pulls each needed SRA run, converts to fastq.gz, cleans up.
cd /tmp/abf
FT=/tmp/sratoolkit.3.4.1-ubuntu64/bin/fasterq-dump
mkdir -p sra_tmp
LOG=/tmp/abf/dl_sra.log

# salvage any leftover uncompressed fastq from a previous run
for f in sra_tmp/*.fastq; do
  [ -s "$f" ] || continue
  b=$(basename "$f" .fastq)
  if [ ! -s "fastq/$b.fastq.gz" ]; then
    gzip -c "$f" > "fastq/$b.fastq.gz" && echo "SALVAGE $b" >> $LOG
  fi
  rm -f "$f"
done

while read -r S; do
  [ -z "$S" ] && continue
  if [ -s "fastq/$S.fastq.gz" ]; then echo "SKIP $S" >> $LOG; continue; fi
  rm -f sra_tmp/${S}*.fastq
  $FT --split-files -e 2 -O sra_tmp "$S" > /dev/null 2>&1
  f=$(ls sra_tmp/${S}*.fastq 2>/dev/null | head -1)
  if [ -s "$f" ]; then
    gzip -c "$f" > "fastq/$S.fastq.gz" && echo "OK $S $(du -h fastq/$S.fastq.gz | cut -f1)" >> $LOG
  else
    echo "FAIL $S" >> $LOG
  fi
  rm -f sra_tmp/${S}*.fastq
done < /tmp/abf/need_sra.txt
echo "DONE" >> $LOG
