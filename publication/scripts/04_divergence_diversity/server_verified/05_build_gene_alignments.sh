#!/usr/bin/env bash
set -euo pipefail

PY="/path/to/workspace/Kuster2026_genomewide_reanalysis/05_build_gene_alignments.py"

echo "===== Step 05 started ====="
date

python3 "$PY"

echo "===== Step 05 finished ====="
date
