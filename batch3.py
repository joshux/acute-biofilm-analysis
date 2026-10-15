#!/usr/bin/env python3
"""batch3.py — run dl_and_quant3 over a list of SRRs, one at a time, logging progress."""
import json, os, subprocess, sys, time

os.chdir("/agent/workspace/abf")
listfile = sys.argv[1]
mode = sys.argv[2] if len(sys.argv) > 2 else "500k"
srrs = [l.strip() for l in open(listfile) if l.strip()]

done = 0
for i, s in enumerate(srrs, 1):
    if os.path.exists(f"final3/{s}.json"):
        done += 1; continue
    t0 = time.time()
    subprocess.run(f"python3 -u dl_and_quant3.py {mode} {s}",
                   shell=True, timeout=7200)
    ok = os.path.exists(f"final3/{s}.json")
    print(f"[{i}/{len(srrs)}] {s} {'OK' if ok else 'FAIL'} {time.time()-t0:.0f}s",
          flush=True)
print(f"BATCH DONE {done}/{len(srrs)} pre-existing", flush=True)
