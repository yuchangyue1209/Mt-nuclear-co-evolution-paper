#!/usr/bin/env bash
set -euo pipefail

source "${1:?Provide config/paths.sh}"
CONSENSUS_DIR="$DATA_ROOT/mt_consensus"
OUT_DIR="$RESULTS_ROOT/mitochondrial_phylogeny"
mkdir -p "$OUT_DIR"
mafft --auto "$CONSENSUS_DIR/all_populations.fasta" > "$OUT_DIR/mtDNA.aligned.fasta"
iqtree2 -s "$OUT_DIR/mtDNA.aligned.fasta" -m GTR+G -bb 1000 -alrt 1000 \
  -nt AUTO -pre "$OUT_DIR/mtDNA"
