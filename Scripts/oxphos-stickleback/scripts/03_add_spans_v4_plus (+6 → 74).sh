#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_env.sh"
source "$(dirname "$0")/utils.sh"

# five supplemental loci you validated earlier (example; adjust if needed)
cat > "${OUTDIR}/add5.spans.bed" <<EOF
chrV    14308310 14313226 uqcrc2b               . -
chrVI   17225319 17227730 ensgacg00000011687    . +
chrXI   14075448 14076086 ensgacg00000013583    . +
chrXVI   9920476  9922109 ensgacg00000004591    . +
chrXVII  8839771  8844687 si:dkey-31b16.7       . -
EOF

SPAN_PLUS="${OUTDIR}/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus.bed"
sort -k1,1 -k2,2n "$SPAN_BASE" "${OUTDIR}/add5.spans.bed" > "$SPAN_PLUS"

# merge maps: try1 + v3 + manual mapping corresponding to the 5 entries (if any)
MAP_MERGED="${OUTDIR}/ortholog_map.subunits.tryDIAMOND.v4_plus.tsv"
cat "$MAP_TRY1" "${OUTDIR}/ortholog_map.from_diamond.v3.tsv" 2>/dev/null \
| awk 'BEGIN{FS=OFS="\t"} NR==1{print; next}
       {key=toupper($1) FS toupper($2); if(!seen[key]++){print}}' > "$MAP_MERGED"

OUT="${OUTDIR}/nuOXPHOS_vs_human.v3map_mix_DIAv4_plus.nuclear_only"
python /work/cyu/oxphos_from_ref_no_biomart/diff_oxphos_sets_mapped.py \
  --span-bed "$SPAN_PLUS" \
  --ortholog-map "$MAP_MERGED" \
  --out-prefix "$OUT" \
  --ignore-mt --keep-ab

summ_line "[add 5 spans (v4_plus)]" "${OUT}.summary.tsv"  # → 74
