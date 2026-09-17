#!/usr/bin/env bash
set -euo pipefail

CODE_ROOT="/path/to/workspace/Kuster2026_genomewide_reanalysis"
DATA_ROOT="/path/to/data/genomewide_codeml_kuster"

# Existing genome-wide nuclear M0 outputs (15,241 variable genes).
NU_RESULTS="$DATA_ROOT/07_codeml_genomewide/results"
NU_GENES="$DATA_ROOT/07_codeml_genomewide/codeml_variable_genes.txt"

# Required: rerun the 13 mtPCGs under codeml M0 on the exact same 27-tip PAML
# topology used for the nuclear genes. Each directory must contain an mlc file:
#   <MT_RESULTS>/ATP6/mlc ... <MT_RESULTS>/ND6/mlc
MT_RESULTS="$DATA_ROOT/10_erc_mt13_codeml/results"

METADATA="$CODE_ROOT/stickleback_mapping/stickleback_gene_classification.final.tsv"
SCRIPT="$CODE_ROOT/genomewide_erc/run_genomewide_erc.py"
OUT="$DATA_ROOT/10_erc_results/mtPCG_composite_t_spearman"

python3 "$SCRIPT" \
  --nuclear-results "$NU_RESULTS" \
  --nuclear-genes "$NU_GENES" \
  --mt-results "$MT_RESULTS" \
  --metadata "$METADATA" \
  --out-dir "$OUT" \
  --metric t \
  --method spearman \
  --expected-nuclear 15241

echo "[OK] ERC outputs: $OUT"
