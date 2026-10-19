#!/usr/bin/env bash
set -euo pipefail
# Reuses existing normalized caller outputs; preserves historical intersections.
source "${1:?Provide config/paths.sh}"
MIN_DP="${MIN_DP:-10}"; MAX_DP="${MAX_DP:-5000}"; MIN_QUAL="${MIN_QUAL:-30}"
OUT="$DATA_ROOT/mt_variants/genotype_checked"
mkdir -p "$OUT"
shopt -s nullglob
vcfs=("$DATA_ROOT/mt_variants/bcftools/"*.vcf.gz)
(( ${#vcfs[@]} )) || { echo 'No bcftools VCFs found' >&2; exit 1; }
for bcf in "${vcfs[@]}"; do
  sample=$(basename "$bcf" .vcf.gz)
  gatk="$DATA_ROOT/mt_variants/gatk_vcf/${sample}.vcf.gz"
  [[ -s "$gatk" ]] || { echo "Missing GATK VCF: $gatk" >&2; exit 1; }
  dir="$OUT/$sample"
  mkdir -p "$dir"
  # Explicitly select the population before checking its genotype. With split
  # biallelic records, GT=alt must support this record's ALT, not another ALT.
  for caller in bcf gatk; do
    input="$bcf"; [[ "$caller" != gatk ]] || input="$gatk"
    bcftools view -s "$sample" -Ou "$input" |
      bcftools norm -f "$MITO_REFERENCE_FASTA" -m -any -Ou |
      bcftools view -i 'GT="alt"' -Ou |
      bcftools sort -Oz -o "$dir/${caller}.alt.vcf.gz"
    bcftools index -f "$dir/${caller}.alt.vcf.gz"
  done
  # -c none requires exact REF/ALT matches; retain the BCFtools genotype.
  bcftools isec -c none -n=2 -w1 -Oz \
    -o "$dir/shared.vcf.gz" "$dir/bcf.alt.vcf.gz" "$dir/gatk.alt.vcf.gz"
  bcftools view -v snps -m2 -M2 \
    -i "QUAL>=${MIN_QUAL} && FORMAT/DP>=${MIN_DP} && FORMAT/DP<=${MAX_DP}" \
    -Oz -o "$dir/pass_snps.vcf.gz" "$dir/shared.vcf.gz"
  bcftools index -f "$dir/pass_snps.vcf.gz"
done
