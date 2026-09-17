#!/usr/bin/env bash
set -euo pipefail

# Historical mtDNA GATK branch. Pool size is used as HaplotypeCaller ploidy.
# Calls are retained only when supported by the independent bcftools branch.
source "${1:?Provide config/paths.sh}"
OUT_DIR="$DATA_ROOT/mt_variants/gatk_gvcf"
BAM_DIR="$DATA_ROOT/mt_bam"
mkdir -p "$OUT_DIR"

awk -F '\t' 'NR>1 && $1!="" && $2~/^[0-9]+$/ {print $1, $2}' "$POOL_INFO" |
while read -r sample pool_size; do
  bam="$BAM_DIR/${sample}.mt.sorted.bam"
  [[ -f "$bam" ]] || { echo "Missing BAM: $bam" >&2; continue; }
  gatk HaplotypeCaller -R "$MITO_REFERENCE_FASTA" -I "$bam" \
    -O "$OUT_DIR/${sample}.g.vcf.gz" -ERC GVCF --ploidy "$pool_size" \
    --minimum-mapping-quality 20 -mbq 13 \
    --indel-size-to-eliminate-in-ref-model 12 \
    -G AS_StandardAnnotation -G StandardAnnotation --sample-name "$sample"
done
