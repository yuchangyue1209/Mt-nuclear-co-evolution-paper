# 01_add_putative_and_diamond_spans.sh
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
