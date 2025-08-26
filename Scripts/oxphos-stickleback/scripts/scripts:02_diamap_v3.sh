#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_env.sh"
source "$(dirname "$0")/utils.sh"

DIAM="${OUTDIR}/missing_subunits.diamond.tsv"     # your existing file
MAP_D3="${OUTDIR}/ortholog_map.from_diamond.v3.tsv"

python "${PYSRC}/diamond_to_map_resolve_gene_v4.py" \
  --diamond "$DIAM" \
  --gff "$GFF" \
  --span "$SPAN_BASE" \
  --out-map "$MAP_D3"

# merge with try1
MAP_MERGED="${OUTDIR}/ortholog_map.subunits.tryDIAMOND.v3.tsv"
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR{if(NR==1){print; next} key=toupper($1) FS toupper($2); if(!seen[key]++){print}; next}
     FNR==1{next} {key=toupper($1) FS toupper($2); if(!seen[key]++){print}}' \
     "$MAP_TRY1" "$MAP_D3" > "$MAP_MERGED"

OUT="${OUTDIR}/nuOXPHOS_vs_human.v3map_mix_DIAv3.nuclear_only"
python /work/cyu/oxphos_from_ref_no_biomart/diff_oxphos_sets_mapped.py \
  --span-bed "$SPAN_BASE" \
  --ortholog-map "$MAP_MERGED" \
  --out-prefix "$OUT" \
  --ignore-mt --keep-ab

summ_line "[DIAmapping v3]" "${OUT}.summary.tsv"  # → 68
