#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_env.sh"
source "$(dirname "$0")/utils.sh"

OUT="${OUTDIR}/nuOXPHOS_vs_human.v3map.nuclear_only"

python /work/cyu/oxphos_from_ref_no_biomart/diff_oxphos_sets_mapped.py \
  --span-bed "$SPAN_BASE" \
  --ortholog-map "$MAP_TRY1" \
  --out-prefix "$OUT" \
  --ignore-mt --keep-ab

summ_line "[baseline]" "${OUT}.summary.tsv"   # → OXPHOS subunits … 67 …
