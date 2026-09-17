#!/usr/bin/env bash
set -euo pipefail

# Jointly genotype the GATK branch, then emit normalized one-sample VCFs for
# comparison with the independent bcftools calls.
source "${1:?Provide config/paths.sh}"
GVCF_DIR="$DATA_ROOT/mt_variants/gatk_gvcf"
DB_DIR="$DATA_ROOT/mt_variants/gatk_genomicsdb"
JOINT_DIR="$DATA_ROOT/mt_variants/gatk_joint"
SINGLE_DIR="$DATA_ROOT/mt_variants/gatk_vcf"
SAMPLE_MAP="$GVCF_DIR/sample_map.tsv"
mkdir -p "$JOINT_DIR" "$SINGLE_DIR"

shopt -s nullglob
gvcfs=("$GVCF_DIR"/*.g.vcf.gz)
(( ${#gvcfs[@]} > 0 )) || { echo "No GVCFs found" >&2; exit 1; }
: > "$SAMPLE_MAP"
for gvcf in "${gvcfs[@]}"; do
  sample="$(basename "$gvcf" .g.vcf.gz)"
  printf '%s\t%s\n' "$sample" "$gvcf" >> "$SAMPLE_MAP"
done
[[ ! -e "$DB_DIR" ]] || { echo "Remove existing GenomicsDB first: $DB_DIR" >&2; exit 1; }

gatk GenomicsDBImport --sample-name-map "$SAMPLE_MAP" \
  --genomicsdb-workspace-path "$DB_DIR" -L "$MITO_CONTIG" \
  --genomicsdb-vcf-buffer-size 1048576 --reader-threads 4
gatk GenotypeGVCFs -R "$MITO_REFERENCE_FASTA" -V "gendb://$DB_DIR" \
  -O "$JOINT_DIR/all_samples.vcf.gz"

while IFS=$'\t' read -r sample _; do
  bcftools view -s "$sample" -Ou "$JOINT_DIR/all_samples.vcf.gz" |
    bcftools norm -m -any -f "$MITO_REFERENCE_FASTA" \
      -Oz -o "$SINGLE_DIR/${sample}.vcf.gz"
  bcftools index -f "$SINGLE_DIR/${sample}.vcf.gz"
done < "$SAMPLE_MAP"
