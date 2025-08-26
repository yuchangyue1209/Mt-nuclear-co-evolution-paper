#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_env.sh"
source "$(dirname "$0")/utils.sh"

SPAN_IN="${OUTDIR}/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus_tblastn.bed"
MISS7="${OUTDIR}/nuOXPHOS_vs_human.v3map_mix_DIA_TBL.nuclear_only.OXPHOS_subunits.missing.txt"
HUMMISS="${OUTDIR}/human_missing_subunits.fa"

# 1) subset FASTA to the current 7
SUB7="${OUTDIR}/human_missing_subunits.7.fa"
: > "$SUB7"
while read -r g; do
  [ -z "$g" ] && continue
  awk -v G="$g" 'BEGIN{IGNORECASE=1; keep=0}
       /^>/{keep=(tolower(substr($0,2))==tolower(G)); if(keep) print ">" G; next}
       {if(keep) print}' "$HUMMISS" >> "$SUB7"
done < "$MISS7"

# 2) relaxed tBLASTn
TBL2="${DBROOT}/missing7.tblastn.tsv"
tblastn -query "$SUB7" -db "$DBROOT/stk" -out "$TBL2" \
  -evalue 1 -word_size 2 -seg yes -comp_based_stats F \
  -max_target_seqs 400 -max_hsps 400 \
  -outfmt '6 qseqid sseqid pident length qlen sstart send evalue bitscore'

# 3) cluster with relaxed cutoffs (qcov>=0.30, pident>=0.15)
SPAN_TBL2="${OUTDIR}/tblastn_putative7.SPAN.bed"
MAP_TBL2="${OUTDIR}/tblastn_putative7.map.tsv"
python "${PYSRC}/tblastn_to_spans.py" \
  --tblastn "$TBL2" \
  --span_out "$SPAN_TBL2" \
  --map_out "$MAP_TBL2" \
  --gap 3000 --pad 200 --min_qcov 0.30 --min_pident 0.15

# 4) merge & diff
SPAN_OUT="${OUTDIR}/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus_tblastn_x2.bed"
sort -k1,1 -k2,2n "$SPAN_IN" "$SPAN_TBL2" > "$SPAN_OUT"

MAP_MERGED="${OUTDIR}/ortholog_map.subunits.tryDIAMOND.tbl2.tsv"
cat "$OUTDIR/ortholog_map.subunits.tryDIAMOND.tbl1.tsv" "$MAP_TBL2" \
| awk 'BEGIN{FS=OFS="\t"} NR==1{print; next}
       {key=toupper($1) FS toupper($2); if(!seen[key]++){print}}' > "$MAP_MERGED"

OUT="${OUTDIR}/nuOXPHOS_vs_human.v3map_mix_DIA_TBLx2.nuclear_only"
python /work/cyu/oxphos_from_ref_no_biomart/diff_oxphos_sets_mapped.py \
  --span-bed "$SPAN_OUT" \
  --ortholog-map "$MAP_MERGED" \
  --out-prefix "$OUT" \
  --ignore-mt --keep-ab

summ_line "[tBLASTn round2 relaxed]" "${OUT}.summary.tsv"  # → 85
