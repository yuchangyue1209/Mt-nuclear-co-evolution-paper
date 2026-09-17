#!/usr/bin/env bash
set -euo pipefail

# Independent haploid bcftools branch used for cross-caller validation.
source "${1:?Provide config/paths.sh}"
BAM_DIR="$DATA_ROOT/mt_bam"
OUT_DIR="$DATA_ROOT/mt_variants/bcftools"
mkdir -p "$OUT_DIR"
shopt -s nullglob
bams=("$BAM_DIR"/*.mt.sorted.bam)
(( ${#bams[@]} > 0 )) || { echo "No mitochondrial BAMs found" >&2; exit 1; }

for bam in "${bams[@]}"; do
  sample="$(basename "$bam" .mt.sorted.bam)"
  raw="$OUT_DIR/${sample}.raw.vcf.gz"
  final="$OUT_DIR/${sample}.vcf.gz"
  bcftools mpileup -Ou -f "$MITO_REFERENCE_FASTA" -d 5000 \
    -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP "$bam" |
    bcftools call --ploidy 1 -mv -Oz -o "$raw"
  bcftools norm -m -any -f "$MITO_REFERENCE_FASTA" -Oz -o "$final" "$raw"
  bcftools index -f "$final"
  rm -f "$raw"
done
