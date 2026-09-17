#!/usr/bin/env bash
set -euo pipefail

# Retain variants called by both branches and apply historical filters.
# Input GATK files must be normalized per-sample VCFs, not GVCFs.
source "${1:?Provide config/paths.sh}"
GATK_DIR="$DATA_ROOT/mt_variants/gatk_vcf"
BCF_DIR="$DATA_ROOT/mt_variants/bcftools"
ISEC_DIR="$DATA_ROOT/mt_variants/intersections"
OUT_DIR="$DATA_ROOT/mt_variants/overlap_filtered"
MIN_DP="${MIN_DP:-10}"; MAX_DP="${MAX_DP:-5000}"; MIN_QUAL="${MIN_QUAL:-30}"
mkdir -p "$ISEC_DIR" "$OUT_DIR"
shopt -s nullglob
vcfs=("$BCF_DIR"/*.vcf.gz)
(( ${#vcfs[@]} > 0 )) || { echo "No bcftools VCFs found" >&2; exit 1; }

for bcf_vcf in "${vcfs[@]}"; do
  sample="$(basename "$bcf_vcf" .vcf.gz)"
  gatk_vcf="$GATK_DIR/${sample}.vcf.gz"
  [[ -f "$gatk_vcf" ]] || { echo "Missing normalized GATK VCF: $gatk_vcf" >&2; continue; }
  bcftools index -f "$gatk_vcf"
  sample_isec="$ISEC_DIR/$sample"
  [[ ! -e "$sample_isec" ]] || { echo "Remove existing directory first: $sample_isec" >&2; exit 1; }
  bcftools isec "$bcf_vcf" "$gatk_vcf" -p "$sample_isec" -Oz
  bcftools view -i "QUAL>=${MIN_QUAL} && FORMAT/DP>=${MIN_DP} && FORMAT/DP<=${MAX_DP}" \
    -Oz -o "$OUT_DIR/${sample}.overlap.vcf.gz" "$sample_isec/0002.vcf.gz"
  bcftools index -f "$OUT_DIR/${sample}.overlap.vcf.gz"
done
