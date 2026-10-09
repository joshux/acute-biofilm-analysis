#!/bin/bash
# Process one sample: quantify3.py -> <SRR>_<org>_q3.json
cd /tmp/abf
S=$1; O=$2
[ -s "${S}_${O}_q3.json" ] && { echo "SKIP $S"; exit 0; }
[ -s "fastq/${S}.fastq.gz" ] || { echo "NOFASTQ $S"; exit 0; }
NPERM=3000 python3 quantify3.py "$S" "$O" 20 > /dev/null 2>&1
[ -s "${S}_${O}_q3.json" ] && echo "OK $S $O" || echo "FAIL $S $O"
