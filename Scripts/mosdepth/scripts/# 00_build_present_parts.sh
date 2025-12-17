# 00_build_present_parts.sh
#!/usr/bin/env bash
set -euo pipefail

# === paths ===
OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
BACK="$OUTDIR/OXPHOS_merged_status.backfilled.tsv"
PARTS6="$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.sorted.bed"
[ -s "$PARTS6" ] || PARTS6="$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.bed"

QCOUT="$OUTDIR/cnv_qc_backfilled_nuclear"
mkdir -p "$QCOUT"

# 0) checks
for f in "$BACK" "$PARTS6"; do
  [ -s "$f" ] || { echo "❌ missing: $f"; exit 1; }
done

# 1) fish names that are present=yes (confident + putative); lowercased unique
awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++){
    if($i=="present_in_fish") P=i
    else if($i=="stickleback_name") S=i
  } next
}
($P=="yes" && $S!=""){
  n=split($S,a,";")
  for(j=1;j<=n;j++){
    s=a[j]; gsub(/^ +| +$/,"",s)
    if(s!="") print tolower(s)
  }
}' "$BACK" | sort -u > "$QCOUT/fish_present_all.list"

# 2) pick all parts for those genes -> 4-col bed (chr start end name|exonN)
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6}' "$PARTS6" > "$QCOUT/parts.src.6col.bed"

awk 'NR==FNR{keep[$1]=1; next}
     {k=tolower($4); if(k in keep) print}' \
  "$QCOUT/fish_present_all.list" "$QCOUT/parts.src.6col.bed" \
> "$QCOUT/parts.present.6col.bed"

awk 'BEGIN{FS=OFS="\t"} {i[$4]++; print $1,$2,$3,$4"|exon"i[$4] }' \
  "$QCOUT/parts.present.6col.bed" \
| tr -d '\r' > "$QCOUT/parts.present.uid.4col.bed"

# sanity
awk 'BEGIN{FS=OFS="\t"} ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/){print "BADINT", NR, $0}' "$QCOUT/parts.present.uid.4col.bed" | sed -n '1,5p'
awk 'BEGIN{FS=OFS="\t"} ($3<=$2){print "BADLEN", NR, $0}' "$QCOUT/parts.present.uid.4col.bed" | sed -n '1,5p'
echo "✅ present parts -> $QCOUT/parts.present.uid.4col.bed"


## 01_add_putative_and_diamond_spans.sh
#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
QCOUT="$OUTDIR/cnv_qc_backfilled_nuclear"
SPAN="$OUTDIR/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus_tblastn.bed" # 同一张SPAN：注释 + putative + DIAMOND
STATUS="$OUTDIR/OXPHOS_merged_status.backfilled.tsv"

[ -s "$QCOUT/parts.present.uid.4col.bed" ] || { echo "❌ run 00_build_present_parts.sh first"; exit 1; }
[ -s "$SPAN" ] || { echo "❌ missing SPAN: $SPAN"; exit 1; }

# 1) aliases that still need SPAN rows (present=yes, have names, but CN 还没算过)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[tolower($i)]=i; next}
tolower($h["present_in_fish"])=="yes" && $h["stickleback_name"]!=""
{ n=split($h["stickleback_name"],a,";"); for(i=1;i<=n;i++){s=a[i]; gsub(/^ +| +$/,"",s); if(s!="") print tolower(s)} }' \
  "$STATUS" | sort -u > "$QCOUT/aliases_from_status.list"

# 2) find which aliases already exist in parts.present (strip |exonN)
cut -f4 "$QCOUT/parts.present.uid.4col.bed" \
| awk -F'|' '{print tolower($1)}' | sort -u > "$QCOUT/aliases_in_present.list"

# 3) need-to-add = aliases_from_status - aliases_in_present
grep -Fvwf "$QCOUT/aliases_in_present.list" "$QCOUT/aliases_from_status.list" \
  > "$QCOUT/aliases_need_span.list" || true

# 4) pull matching rows from SPAN and tag as |span
awk 'NR==FNR{need[$1]=1; next} BEGIN{FS=OFS="\t"} {k=tolower($4); if(k in need) print $1,$2,$3,$4"|span"}' \
  "$QCOUT/aliases_need_span.list" "$SPAN" \
> "$QCOUT/span_additions.from_DIAMOND.4col.bed"

# 5) combine parts + spans (present + putative + diamond-span)
cat "$QCOUT/parts.present.uid.4col.bed" "$QCOUT/span_additions.from_DIAMOND.4col.bed" \
| sort -k1,1 -k2,2n > "$QCOUT/parts.present+putative+spanDIAMOND.4col.bed"

echo "✅ COMB parts -> $QCOUT/parts.present+putative+spanDIAMOND.4col.bed"

# 02_run_mosdepth.sh
#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
QCOUT="$OUTDIR/cnv_qc_backfilled_nuclear"
COMB_PARTS="$QCOUT/parts.present+putative+spanDIAMOND.4col.bed"
POOL_INFO="/work/cyu/poolseq/pool_info.txt"
INPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
REFERENCE="/work/cyu/stickleback_nuclear_only.fa"

for f in "$COMB_PARTS" "$POOL_INFO" "$REFERENCE"; do
  [ -s "$f" ] || { echo "❌ missing: $f"; exit 1; }
done
command -v mosdepth >/dev/null || { echo "❌ mosdepth not found"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not found"; exit 1; }

tail -n +2 "$POOL_INFO" | awk '{print $1}' | while read -r SAMPLE; do
  BAM="${INPUT_DIR}/${SAMPLE}_rg.bam"
  [ -s "$BAM" ] || { echo "⚠️ skip $SAMPLE (no BAM)"; continue; }
  samtools quickcheck "$BAM" 2>/dev/null || samtools index -@8 "$BAM"
  echo "-> mosdepth on $SAMPLE"
  mosdepth -n -x -Q 1 --by "$COMB_PARTS" -f "$REFERENCE" "$QCOUT/${SAMPLE}.COMB" "$BAM"
done

echo "✅ mosdepth done."

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

# 06_make_safe_subunits.sh
#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
IN="$OUTDIR/OXPHOS_merged_status.backfilled.withCN.plus.tsv"
OUT_SAFE="$OUTDIR/subunits.safe_for_dn_ds.nuclear.txt"

awk -F'\t' '
BEGIN{OFS="\t"; 
  want_min=0.5; want_max=1.5
}
NR==1{
  for(i=1;i<=NF;i++){h[tolower($i)]=i}
  need="human_symbol role present_in_fish evidence cn_min cn_max"
  split(need,a," ")
  for(j in a){ if(!(a[j] in h)) { print "Missing col:", a[j] > "/dev/stderr"; exit 1 } }
  next
}
{
  role=$(h["role"]); pres=$(h["present_in_fish"]); ev=$(h["evidence"])
  cnmin=$(h["cn_min"]); cnmax=$(h["cn_max"])
  if(role!="subunit" || pres!="yes") next
  if(ev ~ /tBLASTn_putative/i) next        # 可改：若想包含 putative，注释掉这行
  if(cnmin=="" || cnmax=="") next
  if(cnmin+0.0 < want_min) next
  if(cnmax+0.0 > want_max) next
  print $h["human_symbol"]
}' "$IN" | sort -u > "$OUT_SAFE"

echo "✅ safe subunits -> $OUT_SAFE (n=$(wc -l < "$OUT_SAFE"))"

# RUN_ALL.sh  —— 一键顺序跑完
#!/usr/bin/env bash
set -euo pipefail
bash 00_build_present_parts.sh
bash 01_add_putative_and_diamond_spans.sh
bash 02_run_mosdepth.sh
python 03_collect_regions_to_tables.py --regions-glob "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/*.COMB.regions.bed.gz" \
  --out-depth "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.comb_gene_depth.tsv"
python 04_depth_to_cn.py \
  --depth "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.comb_gene_depth.tsv" \
  --pool-info "/work/cyu/poolseq/pool_info.txt" \
  --out-cn "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.comb_gene_cn.tsv"
python 05_join_cnv_into_status.py \
  "/work/cyu/oxphos_from_ref_no_biomart/06_igv/OXPHOS_merged_status.backfilled.tsv" \
  "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.comb_gene_cn.tsv" \
  > "/work/cyu/oxphos_from_ref_no_biomart/06_igv/OXPHOS_merged_status.backfilled.withCN.plus.tsv"
bash 06_make_safe_subunits.sh


#!/usr/bin/env bash
# refresh_cn_plus2.sh
# Purpose: For present==yes but missing CN, pull SPAN coordinates for aliases/IDs,
#          rebuild COMB parts, rerun mosdepth for all pools, recompute depth & CN,
#          and write back into OXPHOS_merged_status.backfilled.withCN.plus2.tsv

set -euo pipefail

# -------------------- paths (edit as needed) --------------------
OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
QCOUT="$OUTDIR/cnv_qc_backfilled_nuclear"
STATUS="$OUTDIR/OXPHOS_merged_status.backfilled.withCN.plus.tsv"      # 上一轮（plus）的主表
SPAN="$OUTDIR/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus_tblastn.bed"      # 含 annotation/DIAMOND/tBLASTn 的 SPAN
COMB_OLD="$QCOUT/parts.present+putative.4col.bed"                      # 之前的 present+putative 4列 BED
COMB_NEW="$QCOUT/parts.present+putative+spanDIAMOND.4col.bed"          # 本轮要生成的合并 BED
POOL_INFO="/work/cyu/poolseq/pool_info.txt"
INPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups"                     # ${SAMPLE}_rg.bam
REF="/work/cyu/stickleback_nuclear_only.fa"

mkdir -p "$QCOUT"

# -------------------- sanity checks --------------------
for f in "$STATUS" "$SPAN" "$COMB_OLD" "$POOL_INFO" "$REF"; do
  [ -s "$f" ] || { echo "❌ missing: $f"; exit 1; }
done
command -v mosdepth >/dev/null || { echo "❌ mosdepth not found"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not found"; exit 1; }

# -------------------- 1) pick aliases that still lack CN --------------------
# 条件：present=yes 且 stickleback_name 非空 且 CN_min 为空
awk -F'\t' 'NR==1{
  for(i=1;i<=NF;i++) h[tolower($i)]=i; next
}
tolower($(h["present_in_fish"]))=="yes" && $(h["stickleback_name"])!="" && $(h["cn_min"])==""
{
  n=split($(h["stickleback_name"]),a,";")
  for(i=1;i<=n;i++){
    s=a[i]; gsub(/^ +| +$/,"",s); if(s!="") print tolower(s)
  }
}' "$STATUS" | sort -u > "$QCOUT/aliases_need_span.list"

echo "need-span aliases: $(wc -l < "$QCOUT/aliases_need_span.list")"

# -------------------- 2) pull SPAN rows for those aliases -> 4-col with |span --------------------
# SPAN 第4列有基因名/ID（可能大小写不一），我们用小写匹配。
awk 'NR==FNR{need[$1]=1; next}
     BEGIN{FS=OFS="\t"}
{
  key=tolower($4)
  if(key in need){
    print $1,$2,$3,$4"|span"
  }
}' "$QCOUT/aliases_need_span.list" "$SPAN" \
| sort -k1,1 -k2,2n -u > "$QCOUT/span_additions.from_SPAN.4col.bed"

echo "span additions: $(wc -l < "$QCOUT/span_additions.from_SPAN.4col.bed")"

# -------------------- 3) merge with previous COMB parts --------------------
cat "$COMB_OLD" "$QCOUT/span_additions.from_SPAN.4col.bed" \
| sort -k1,1 -k2,2n -u > "$COMB_NEW"

echo "COMB_NEW -> $COMB_NEW"

# -------------------- 4) run mosdepth for all pools on COMB_NEW --------------------
# 产出：$QCOUT/<SAMPLE>.COMB2.regions.bed.gz
tail -n +2 "$POOL_INFO" | awk '{print $1}' | while read -r SAMPLE; do
  BAM="${INPUT_DIR}/${SAMPLE}_rg.bam"
  if [ ! -s "$BAM" ]; then
    echo "⚠️  skip $SAMPLE (missing $BAM)"; continue
  fi
  samtools quickcheck "$BAM" 2>/dev/null || samtools index -@8 "$BAM"
  echo "-> mosdepth (COMB2) on $SAMPLE"
  mosdepth -n -x -Q 1 --by "$COMB_NEW" -f "$REF" \
    "$QCOUT/${SAMPLE}.COMB2" "$BAM"
done

# -------------------- 5) summarize regions -> gene x sample depth table --------------------
# 写入一个小 Python 脚本，按 gene(去掉 |exonN/|span) 聚合均值
PY_SUM="$QCOUT/_summarize_regions_to_depth_plus2.py"
cat > "$PY_SUM" <<'PY'
import os, gzip, sys, glob, statistics
qc = sys.argv[1]
out_tab = sys.argv[2]
files = sorted(glob.glob(os.path.join(qc, "*.COMB2.regions.bed.gz")))
samples = [os.path.basename(f).split(".COMB2.")[0] for f in files]
per_gene = {}  # gene -> sample -> [depths]
for f,s in zip(files,samples):
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if not line or line[0]=="#": continue
            chrom,beg,end,name,depth = line.rstrip("\n").split("\t")
            gene = name.split("|",1)[0]
            per_gene.setdefault(gene,{})
            per_gene[gene].setdefault(s,[]).append(float(depth))
# mean over parts
genes = sorted(per_gene.keys())
with open(out_tab,"w") as out:
    out.write("gene\t" + "\t".join(samples) + "\n")
    for g in genes:
        vals=[]
        for s in samples:
            arr = per_gene[g].get(s, [])
            m = statistics.mean(arr) if arr else ""
            vals.append("" if arr==[] else f"{m:.4f}")
        out.write(g + "\t" + "\t".join(vals) + "\n")
PY

DEP_TAB="$QCOUT/oxphos_cnv_plus2.comb_gene_depth.tsv"
python "$PY_SUM" "$QCOUT" "$DEP_TAB"
echo "[DONE] depth -> $DEP_TAB"

# -------------------- 6) convert depth -> CN (per-gene median across samples = 1) --------------------
PY_CN="$QCOUT/_depth_to_cn_plus2.py"
cat > "$PY_CN" <<'PY'
import sys, csv, statistics
inp, outp = sys.argv[1], sys.argv[2]
rows = list(csv.reader(open(inp), delimiter="\t"))
samples = rows[0][1:]
table = {}
for r in rows[1:]:
    g=r[0]; vals=[]
    for v in r[1:]:
        vals.append(None if v=="" else float(v))
    table[g]=vals
# gene medians
gmed = {}
for g,vals in table.items():
    x=[v for v in vals if v is not None]
    gmed[g] = statistics.median(x) if x else None
# CN = depth / gene_median
def cn_row(vals, med):
    if med is None or med==0: return ["" if v is None else "" for v in vals]
    return ["" if v is None else f"{(v/med):.3f}" for v in vals]
# outliers
def out_stats(cns, thr_low=0.5, thr_high=1.5):
    xs=[(i,v) for i,v in enumerate(cns) if v not in ("",None)]
    mn= min(float(v) for i,v in xs) if xs else None
    mx= max(float(v) for i,v in xs) if xs else None
    outs=[samples[i] for i,v in xs if float(v)<thr_low or float(v)>thr_high]
    return mn, mx, outs
with open(outp,"w") as w:
    w.write("gene\t" + "\t".join([f"CN_{s}" for s in samples]) + "\tCN_min\tCN_max\tN_outlier\toutlier_pools\n")
    for g,vals in table.items():
        cns = cn_row(vals, gmed[g])
        mn,mx,outs = out_stats(cns)
        w.write(g + "\t" + "\t".join(cns) + "\t" +
                ("" if mn is None else f"{mn:.3f}") + "\t" +
                ("" if mx is None else f"{mx:.3f}") + "\t" +
                f"{len(outs)}\t" + (",".join(outs) if outs else "") + "\n")
PY

CN_TAB="$QCOUT/oxphos_cnv_plus2.comb_gene_cn.tsv"
python "$PY_CN" "$DEP_TAB" "$CN_TAB"
echo "[DONE]   CN  -> $CN_TAB"

# -------------------- 7) join CN back to the merged status -> plus2 --------------------
PY_JOIN="$QCOUT/_join_cnv_into_status_plus2.py"
cat > "$PY_JOIN" <<'PY'
import sys, csv
status, cn_tab = sys.argv[1], sys.argv[2]
# load CN
cn = {}
r = csv.DictReader(open(cn_tab), delimiter="\t")
for row in r:
    cn[row["gene"]] = {
        "CN_min": row.get("CN_min",""),
        "CN_max": row.get("CN_max",""),
        "N_outlier": row.get("N_outlier",""),
        "outlier_pools": row.get("outlier_pools",""),
    }
# status join: stickleback_name 里的别名逐个查，合并 min/max/outliers
def split_alias(s):
    return [x.strip() for x in s.split(";") if x.strip()]
def merge_cn(names):
    mins=[]; maxs=[]; outs=set()
    for n in names:
        rec = cn.get(n) or cn.get(n.lower()) or cn.get(n.upper())
        if not rec: continue
        if rec["CN_min"]: mins.append(float(rec["CN_min"]))
        if rec["CN_max"]: maxs.append(float(rec["CN_max"]))
        if rec["N_outlier"] and rec["outlier_pools"]:
            outs.update(rec["outlier_pools"].split(","))
    if not mins and not maxs and not outs: return "","","",""
    mn = f"{min(mins):.3f}" if mins else ""
    mx = f"{max(maxs):.3f}" if maxs else ""
    outlist = sorted([x for x in outs if x])
    return mn, mx, str(len(outlist)), ",".join(outlist)
rd = csv.reader(open(status), delimiter="\t")
hdr = next(rd)
# ensure CN columns
for col in ("CN_min","CN_max","N_outlier","outlier_pools"):
    if col not in hdr: hdr.append(col)
out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
out.writerow(hdr)
hmap = {c:i for i,c in enumerate(hdr)}
for row in rd:
    # position of columns (safe whether exist previously or not)
    try: sidx = hmap["stickleback_name"]
    except KeyError: sidx = None
    names = split_alias(row[sidx]) if sidx is not None and sidx < len(row) and row[sidx] else []
    mn,mx,n,outs = merge_cn(names)
    # append or replace
    def setcol(name,val):
        if name in hmap: row[hmap[name]] = val
        else: row.append(val)
    setcol("CN_min", mn)
    setcol("CN_max", mx)
    setcol("N_outlier", n)
    setcol("outlier_pools", outs)
    out.writerow(row)
PY

STATUS2="$OUTDIR/OXPHOS_merged_status.backfilled.withCN.plus2.tsv"
python "$PY_JOIN" "$STATUS" "$CN_TAB" > "$STATUS2"
echo "[DONE] merged status -> $STATUS2"

# small preview
sed -n '1,20p' "$STATUS2" | column -t
