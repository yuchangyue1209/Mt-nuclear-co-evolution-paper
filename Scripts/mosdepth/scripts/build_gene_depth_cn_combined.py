#!/usr/bin/env python3
"""
Summarize mosdepth *.regions.bed.gz (combined confident+putative) into:
  1) gene × pool mean depth table
  2) gene × pool CN table (depth / per-pool median across genes)

Usage:
  python build_gene_depth_cn.py QC_DIR POOL_INFO PARTS4 OUT_DEPTH_TSV OUT_CN_TSV

- QC_DIR      : directory containing <POOL>.COMB.regions.bed.gz
- POOL_INFO   : pool_info.txt (header + SAMPLE ...)
- PARTS4      : 4-col BED used for --by (name like `gene|exonN` or `gene|span`)
- OUT_DEPTH_TSV: output gene depth table
- OUT_CN_TSV   : output gene CN table
"""
import os, sys, gzip, csv, statistics
from collections import defaultdict, OrderedDict

if len(sys.argv) != 6:
    sys.stderr.write("Usage: build_gene_depth_cn.py QC_DIR POOL_INFO PARTS4 OUT_DEPTH_TSV OUT_CN_TSV\n")
    sys.exit(1)

qcdir, poolinfo, parts4, out_dep, out_cn = sys.argv[1:]

# gene list from 4-col BED (col4 -> split by '|', take left part)
genes, seen = [], set()
with open(parts4) as f:
    for ln in f:
        if not ln.strip(): continue
        g = ln.split("\t", 3)[3].strip().split("|", 1)[0]
        if g and g not in seen:
            seen.add(g); genes.append(g)

# pools from POOL_INFO (skip header)
with open(poolinfo) as f:
    pools = [ln.strip().split()[0] for ln in f if ln.strip()]
pools = pools[1:]  # drop header row

def regions_path(pool: str) -> str:
    p1 = os.path.join(qcdir, f"{pool}.COMB.regions.bed.gz")
    if os.path.exists(p1): return p1
    p2 = os.path.join(qcdir, f"{pool}.regions.bed.gz")
    return p2

# gene×pool mean depth
D = { g: OrderedDict((p, None) for p in pools) for g in genes }

for p in pools:
    reg = regions_path(p)
    if not os.path.exists(reg):
        sys.stderr.write(f"[WARN] regions not found for {p}: {reg}\n")
        continue
    acc_sum = defaultdict(float)
    acc_n   = defaultdict(int)
    with gzip.open(reg, "rt") as fh:
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if len(c) < 5: continue
            name, depth_s = c[3], c[4]
            try:
                d = float(depth_s)
            except:
                continue
            g = name.split("|",1)[0]
            acc_sum[g] += d
            acc_n[g]   += 1
    for g in genes:
        n = acc_n.get(g, 0)
        if n > 0:
            D[g][p] = acc_sum[g] / n

# write depth table
with open(out_dep, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["gene"] + pools)
    for g in genes:
        w.writerow([g] + [("" if D[g][p] is None else f"{D[g][p]:.4f}") for p in pools])

# per-pool median across genes (exclude None/0)
pool_median = {}
for p in pools:
    vals = [D[g][p] for g in genes if D[g][p] is not None and D[g][p] > 0]
    pool_median[p] = (statistics.median(vals) if vals else None)

# write CN table
with open(out_cn, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["gene"] + [f"CN_{p}" for p in pools])
    for g in genes:
        row = [g]
        for p in pools:
            md = pool_median[p]
            x  = D[g][p]
            if md is None or x is None:
                row.append("")
            else:
                row.append(f"{x/md:.3f}")
        w.writerow(row)
