# 05_join_cnv_into_status.py
#!/usr/bin/env python3
import sys, csv

if len(sys.argv)<3:
    sys.exit("usage: 05_join_cnv_into_status.py STATUS_IN CN_TABLE > STATUS_OUT")

status_in, cn_tab = sys.argv[1], sys.argv[2]

# load CN summary per gene (CN_min/CN_max/N_outlier/outlier_pools)
cn = {}
with open(cn_tab) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        g=row['gene'].strip()
        if not g: continue
        cn[g.lower()] = {
            'CN_min': row.get('CN_min',''),
            'CN_max': row.get('CN_max',''),
            'N_outlier': row.get('N_outlier',''),
            'outlier_pools': row.get('outlier_pools','')
        }

# join into status; map semicolon alias list -> merge (min/min, max/max, sum outliers, union pools)
def merge_aliases(alias_str):
    if not alias_str: return ("","","","")
    aliases=[x.strip().lower() for x in alias_str.split(";") if x.strip()]
    mins=[]; maxs=[]; outs=0; pools=set()
    for a in aliases:
        rec = cn.get(a)
        if not rec: continue
        if rec['CN_min']!='': mins.append(float(rec['CN_min']))
        if rec['CN_max']!='': maxs.append(float(rec['CN_max']))
        try: outs += int(rec['N_outlier'] or 0)
        except: pass
        if rec['outlier_pools']:
            for p in rec['outlier_pools'].split(","):
                p=p.strip()
                if p: pools.add(p)
    CN_min = f"{min(mins):.3f}" if mins else ""
    CN_max = f"{max(maxs):.3f}" if maxs else ""
    N_out  = str(outs) if outs else ""
    Pools  = ",".join(sorted(pools)) if pools else ""
    return (CN_min, CN_max, N_out, Pools)

with open(status_in) as f:
    rdr = csv.DictReader(f, delimiter="\t")
    hdr = rdr.fieldnames + ["CN_min","CN_max","N_outlier","outlier_pools"]
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(hdr)
    for row in rdr:
        cn_min, cn_max, n_out, pools = merge_aliases(row.get("stickleback_name",""))
        row_out=[row.get(h,"") for h in rdr.fieldnames] + [cn_min,cn_max,n_out,pools]
        w.writerow(row_out)
