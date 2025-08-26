#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_env.sh"
source "$(dirname "$0")/utils.sh"

SPAN="${OUTDIR}/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus_tblastn_x2.bed"
MAP_ALL="${OUTDIR}/ortholog_map.subunits.ALL.tsv"

# merge *all* maps you created along the way (add files if you have more)
cat \
  "$MAP_TRY1" \
  "${OUTDIR}/ortholog_map.from_diamond.v3.tsv" \
  "${OUTDIR}/tblastn_putative.map.tsv" \
  "${OUTDIR}/tblastn_putative7.map.tsv" \
| awk 'BEGIN{FS=OFS="\t"} NR==1{print; next}
       {key=toupper($1) FS toupper($2); if(!seen[key]++){print}}' > "$MAP_ALL"

OUT="${OUTDIR}/nuOXPHOS_vs_human.vFINAL_COMBINED.nuclear_only"
python /work/cyu/oxphos_from_ref_no_biomart/diff_oxphos_sets_mapped.py \
  --span-bed "$SPAN" \
  --ortholog-map "$MAP_ALL" \
  --out-prefix "$OUT" \
  --ignore-mt --keep-ab

summ_line "[FINAL]" "${OUT}.summary.tsv"
echo "Missing list:"
cat "${OUT}.OXPHOS_subunits.missing.txt"
