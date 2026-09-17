#!/usr/bin/env bash
set -euo pipefail

# Run this script on the Mac from this repository. It downloads the exact
# finalized stage 00–07 scripts currently present on the stickleback server.

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
REMOTE="user@server"
REMOTE_ROOT="/path/to/workspace/Kuster2026_genomewide_reanalysis"

mkdir -p "$HERE/server_verified"

FILES=(
  00_prepare_genomewide_targets.sh
  01_call_variants_genomewide_hap1.sh
  03_make_masks_and_filtered_vcf.sh
  04_build_masked_consensus_cds.sh
  05_build_gene_alignments.py
  05_build_gene_alignments.sh
  06_prepare_codeml.py
  06_run_codeml_pilot.sh
  07_codeml_worker.sh
  07_run_codeml_genomewide.sh
)

for file in "${FILES[@]}"; do
  echo "[sync] $file"
  scp "$REMOTE:$REMOTE_ROOT/$file" "$HERE/server_verified/$file"
done

echo "[OK] downloaded to $HERE/server_verified"

