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
