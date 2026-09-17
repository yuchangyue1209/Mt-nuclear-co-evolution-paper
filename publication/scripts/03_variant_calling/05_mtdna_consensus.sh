#!/usr/bin/env bash
set -euo pipefail

source "${1:?Provide config/paths.sh}"
VCF_DIR="$DATA_ROOT/mt_variants/overlap_filtered"
OUT_DIR="$DATA_ROOT/mt_consensus"
mkdir -p "$OUT_DIR"
shopt -s nullglob
vcfs=("$VCF_DIR"/*.overlap.vcf.gz)
(( ${#vcfs[@]} > 0 )) || { echo "No overlap VCFs found" >&2; exit 1; }
: > "$OUT_DIR/all_populations.fasta"

for vcf in "${vcfs[@]}"; do
  sample="$(basename "$vcf" .overlap.vcf.gz)"
  fasta="$OUT_DIR/${sample}.fasta"
  bcftools consensus -f "$MITO_REFERENCE_FASTA" -M N -s "$sample" "$vcf" |
    awk -v sample="$sample" '/^>/{print ">" sample; next} {print}' > "$fasta"
  cat "$fasta" >> "$OUT_DIR/all_populations.fasta"
done
