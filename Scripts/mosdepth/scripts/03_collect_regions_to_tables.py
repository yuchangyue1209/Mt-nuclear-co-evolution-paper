# 03_collect_regions_to_tables.py
#!/usr/bin/env python3
import argparse, glob, gzip, os, sys
from collections import defaultdict

p=argparse.ArgumentParser()
p.add_argument("--regions-glob", required=True, help="comma-separated globs to *.regions.bed.gz")
p.add_argument("--out-depth", required=True)
args=p.parse_args()

# gather files
files=[]
for patt in args.regions_glob.split(","):
    files.extend(glob.glob(patt))
files=sorted(set(f for f in files if f.endswith(".regions.bed.gz")))
if not files:
    sys.exit("no regions.bed.gz found")

# sample name from basename (remove .regions.bed.gz)
def samp(path):
    b=os.path.basename(path)
    return b.replace(".regions.bed.gz","")

# mean depth per gene per sample
genes=set()
table=defaultdict(lambda: defaultdict(lambda: [0.0,0]))  # gene -> sample -> [sum, n]

for fp in files:
    s = samp(fp)
    with gzip.open(fp,"rt") as f:
        for ln in f:
            if not ln.strip(): continue
            chrom,start,end,name,depth = ln.rstrip("\n").split("\t")
            gene = name.split("|",1)[0]
            try:
                d = float(depth)
            except:
                continue
            acc = table[gene][s]
            acc[0]+=d; acc[1]+=1
            genes.add(gene)

samples=sorted({samp(f) for f in files})
with open(args.out_depth,"w") as out:
    out.write("gene\t" + "\t".join(samples) + "\n")
    for g in sorted(genes):
        row=[g]
        for s in samples:
            sm, n = table[g][s]
            m = sm/n if n else ""
            row.append(f"{m:.6f}" if m!="" else "")
        out.write("\t".join(row)+"\n")
print(f"[DONE] depth -> {args.out_depth}")
