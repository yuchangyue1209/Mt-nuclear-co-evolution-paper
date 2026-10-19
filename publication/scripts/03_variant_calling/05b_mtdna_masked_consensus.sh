#!/usr/bin/env bash
set -euo pipefail
source "${1:?Provide config/paths.sh}"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MIN_DP="${MIN_DP:-10}"; MAX_DP="${MAX_DP:-5000}"
MIN_MAPQ="${MIN_MAPQ:-20}"; MIN_BASEQ="${MIN_BASEQ:-13}"
MIN_MAJOR_AF="${MIN_MAJOR_AF:-0.80}"; INDEL_PAD="${INDEL_PAD:-5}"
PILEUP_CAP="${PILEUP_CAP:-100000}"
(( PILEUP_CAP > MAX_DP )) || { echo 'PILEUP_CAP must exceed MAX_DP' >&2; exit 1; }
OUT="$DATA_ROOT/mt_consensus_masked"
mkdir -p "$OUT"
[[ -s "${MITO_REFERENCE_FASTA}.fai" ]] || samtools faidx "$MITO_REFERENCE_FASTA"
shopt -s nullglob
dirs=("$DATA_ROOT/mt_variants/genotype_checked/"*/)
(( ${#dirs[@]} )) || { echo 'Run 04b first' >&2; exit 1; }
aggregate="$OUT/all_populations.fasta.tmp"
: > "$aggregate"
for dir in "${dirs[@]}"; do
  sample=$(basename "$dir")
  bam="$DATA_ROOT/mt_bam/${sample}.mt.sorted.bam"
  [[ -s "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }
  [[ -s "${bam}.bai" || -s "${bam%.bam}.bai" ]] || samtools index "$bam"
  work="$OUT/$sample"
  mkdir -p "$work"
  [[ $(bcftools query -l "$dir/pass_snps.vcf.gz") == "$sample" ]] || exit 1
  # No -v/call step: retain reference positions as well as candidate variants.
  # Unreported positions are initialized as masked by the Python helper.
  bcftools mpileup -Ou -f "$MITO_REFERENCE_FASTA" -r "$MITO_CONTIG" \
    -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -d "$PILEUP_CAP" \
    -a FORMAT/DP,FORMAT/AD "$bam" |
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' \
      > "$work/pileup.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    "$dir/pass_snps.vcf.gz" > "$work/pass.tsv"
  : > "$work/candidates.tsv"
  for caller in bcf gatk; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
      "$dir/${caller}.alt.vcf.gz" >> "$work/candidates.tsv"
  done
  python3 "$SCRIPT_DIR/mtdna_mask.py" \
    --fai "${MITO_REFERENCE_FASTA}.fai" --contig "$MITO_CONTIG" \
    --pileup "$work/pileup.tsv" --passed "$work/pass.tsv" \
    --candidates "$work/candidates.tsv" --min-dp "$MIN_DP" \
    --max-dp "$MAX_DP" --major-af "$MIN_MAJOR_AF" --indel-pad "$INDEL_PAD" \
    --bed "$work/mask.bed" --qa "$work/mask_qa.json"
  samtools faidx "$MITO_REFERENCE_FASTA" "$MITO_CONTIG" |
    bcftools consensus \
    -m "$work/mask.bed" -M N -s "$sample" "$dir/pass_snps.vcf.gz" |
    awk -v sample="$sample" '/^>/{print ">" sample; next} {print}' \
    > "$work/consensus.fasta.tmp"
  mv "$work/consensus.fasta.tmp" "$OUT/${sample}.fasta"
  cat "$OUT/${sample}.fasta" >> "$aggregate"
  printf 'min_mapq=%s\nmin_baseq=%s\nmin_dp=%s\nmax_dp=%s\nmajor_af=%s\nindel_pad=%s\npileup_cap=%s\n' \
    "$MIN_MAPQ" "$MIN_BASEQ" "$MIN_DP" "$MAX_DP" "$MIN_MAJOR_AF" "$INDEL_PAD" "$PILEUP_CAP" \
    > "$work/parameters.txt"
done
mv "$aggregate" "$OUT/all_populations.fasta"
echo "Masked sequences: $OUT (historical mt_consensus is unchanged)"
