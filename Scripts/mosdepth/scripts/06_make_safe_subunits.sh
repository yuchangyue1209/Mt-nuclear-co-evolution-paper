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
