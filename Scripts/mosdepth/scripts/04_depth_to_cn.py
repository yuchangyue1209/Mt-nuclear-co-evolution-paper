# 04_depth_to_cn.py
#!/usr/bin/env python3
import argparse, csv, statistics as st

p=argparse.ArgumentParser()
p.add_argument("--depth", required=True, help="output of 03_collect_regions_to_tables.py")
p.add_argument("--pool-info", required=True, help="two columns: SAMPLE\tpool_size(or mapped_reads). Header allowed.")
p.add_argument("--out-cn", required=True)
args=p.parse_args()

# load pool size (mapped reads) as size factors; ignore header if present
sizes={}
with open(args.pool_info) as f:
    for ln in f:
        if not ln.strip(): continue
        a=ln.rstrip("\n").split()
        if not a: continue
        if a[0].lower().startswith("sample"): continue
        try:
            sizes[a[0]] = float("".join([c for c in a[1] if c.isdigit()]))
        except:
            continue

# read depth table
with open(args.depth) as f:
    header=f.readline().rstrip("\n").split("\t")
    samples=header[1:]
    rows=[ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
# scale by mapped reads (median-normalized)
vals=[sizes.get(s,1.0) for s in samples]
med = st.median([v for v in vals if v>0]) if vals else 1.0
sf  = {s:(sizes.get(s,1.0)/med if med>0 else 1.0) for s in samples}

def safe_float(x): 
    try: return float(x)
    except: return None

# normalized depth and CN per gene
with open(args.out_cn,"w") as out:
    cn_cols=[f"CN_{s}" for s in samples]
    out.write("gene\t" + "\t".join(cn_cols) + "\tCN_min\tCN_max\tN_outlier\toutlier_pools\n")
    for a in rows:
        gene=a[0]
        depths=[safe_float(a[i+1]) for i,_ in enumerate(samples)]
        # scale depths by size factor
        nd=[(depths[i]/sf[samples[i]] if depths[i] is not None and sf[samples[i]]>0 else None) for i in range(len(samples))]
        # gene-level median to anchor CN≈1
        med_g = st.median([x for x in nd if x is not None]) if any(x is not None for x in nd) else None
        cn=[]
        for i,x in enumerate(nd):
            if x is None or med_g in (None,0): cn.append(None)
            else: cn.append(x/med_g)
        # outliers
        outs=[]
        for i,c in enumerate(cn):
            if c is None: continue
            if c<0.5 or c>1.5: outs.append(samples[i])
        cn_min = min([c for c in cn if c is not None], default="")
        cn_max = max([c for c in cn if c is not None], default="")
        out.write(gene+"\t")
        out.write("\t".join("" if c is None else f"{c:.3f}" for c in cn))
        out.write(f"\t{'' if cn_min=='' else f'{cn_min:.3f}'}\t{'' if cn_max=='' else f'{cn_max:.3f}'}\t{len(outs)}\t{','.join(outs)}\n")
print(f"[DONE] CN -> {args.out_cn}")
