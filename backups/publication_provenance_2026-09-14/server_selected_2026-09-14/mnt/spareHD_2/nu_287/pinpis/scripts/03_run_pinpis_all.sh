#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/spareHD_2/nu_287/pinpis"
POOL="${ROOT}/meta/pool_sizes.tsv"

POP_NU="${ROOT}/meta/pop_order.nuclear.txt"
POP_MT="${ROOT}/meta/pop_order.mt.txt"

OUT_NU="${ROOT}/results/pinpis_nuclear.tsv"
OUT_MT="${ROOT}/results/pinpis_mt.tsv"
OUT_ALL="${ROOT}/results/pinpis_all_300.tsv"

# 287 nuclear
python3 "${ROOT}/scripts/02_pinpis_from_sync.py" \
  --sync_glob "/mnt/spareHD_2/nu_287/sync/*.sync" \
  --params_suffix ".params" \
  --pool_sizes "$POOL" \
  --pop_order "$POP_NU" \
  --genetic_code nuclear \
  --out_tsv "$OUT_NU"

# 13 mt PCG（你这个目录里已经是 per-gene sync）
python3 "${ROOT}/scripts/02_pinpis_from_sync.py" \
  --sync_glob "/work/cyu/poolseq/PPalign_output/ann_mt_pergene_sync/{ND1,ND2,ND3,ND4,ND4L,ND5,ND6,COX1,COX2,COX3,ATP6,ATP8,CYTB}.sync" \
  --pool_sizes "$POOL" \
  --pop_order "$POP_MT" \
  --genetic_code mt \
  --out_tsv "$OUT_MT"

# merge
{ head -n 1 "$OUT_NU"; tail -n +2 "$OUT_NU"; tail -n +2 "$OUT_MT"; } > "$OUT_ALL"

echo "[OK] nuclear: $OUT_NU"
echo "[OK] mt     : $OUT_MT"
echo "[OK] merged : $OUT_ALL"
echo "[OK] merged lines: $(wc -l < "$OUT_ALL")"
