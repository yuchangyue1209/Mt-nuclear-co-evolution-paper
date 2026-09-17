#!/usr/bin/env bash
set -euo pipefail

REF="/path/to/workspace/stickleback_nuclear_only.fa"
BED="/path/to/data/genomewide_codeml_kuster/00_targets_codeml_ready/stickleback_20348_canonical_cds_merged.bed"
BAM_DIR="/path/to/data/nuclear_with_readgroups"

ROOT="/path/to/data/genomewide_codeml_kuster"
OUT="$ROOT/01_vcf_hap1"
LOG="$ROOT/logs/01_call_variants"

THREADS="${THREADS:-4}"

mkdir -p "$OUT" "$LOG"

if [[ ! -s "$REF" ]]; then
    echo "ERROR: reference not found: $REF" >&2
    exit 1
fi

if [[ ! -s "$BED" ]]; then
    echo "ERROR: target BED not found: $BED" >&2
    exit 1
fi

[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

echo "[01] Reference: $REF"
echo "[01] Target BED: $BED"
echo "[01] Target intervals: $(wc -l < "$BED")"
echo "[01] BAM directory: $BAM_DIR"
echo "[01] Output: $OUT"

echo "[01] Disk space"
df -h "$ROOT"

SAMPLE_LIST="$OUT/samples.list"

find "$BAM_DIR" -maxdepth 1 -type f \
    \( -name '*_rg.bam' -o -name '*_rg.sorted.bam' \) \
    -printf '%f\n' |
sed -E '
    s/_rg\.sorted\.bam$//
    s/_rg\.bam$//
' |
sort -u > "$SAMPLE_LIST"

NSAMPLES=$(wc -l < "$SAMPLE_LIST")
echo "[01] Samples detected: $NSAMPLES"

if [[ "$NSAMPLES" -eq 0 ]]; then
    echo "ERROR: no BAM samples detected" >&2
    exit 1
fi

while read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    if [[ -f "$BAM_DIR/${SAMPLE}_rg.sorted.bam" ]]; then
        BAM="$BAM_DIR/${SAMPLE}_rg.sorted.bam"
    elif [[ -f "$BAM_DIR/${SAMPLE}_rg.bam" ]]; then
        BAM="$BAM_DIR/${SAMPLE}_rg.bam"
    else
        echo "ERROR: BAM missing for $SAMPLE" >&2
        exit 1
    fi

    VCF="$OUT/${SAMPLE}.genomewide_cds.hap1.vcf.gz"
    TMP="$OUT/.${SAMPLE}.tmp.vcf.gz"
    SAMPLE_LOG="$LOG/${SAMPLE}.log"

    if [[ -s "$VCF" && -s "${VCF}.csi" ]]; then
        if bcftools index -n "$VCF" >/dev/null 2>&1; then
            echo "[01] SKIP complete: $SAMPLE"
            continue
        fi
    fi

    echo "[01] START: $SAMPLE"
    echo "[01] BAM: $BAM"

    if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
        echo "[01] Indexing BAM: $SAMPLE"
        samtools index -@ "$THREADS" "$BAM"
    fi

    {
        echo "sample=$SAMPLE"
        echo "bam=$BAM"
        date

        bcftools mpileup \
            --threads "$THREADS" \
            -Ou \
            -f "$REF" \
            -R "$BED" \
            -q 30 \
            -Q 25 \
            -a FORMAT/AD,FORMAT/DP \
            "$BAM" |
        bcftools call \
            --threads "$THREADS" \
            -mv \
            --ploidy 1 \
            -Ou |
        bcftools norm \
            --threads "$THREADS" \
            -f "$REF" \
            -m -any \
            -Ou |
        bcftools sort \
            -Oz \
            -o "$TMP"

        bcftools index -f "$TMP"

        mv "$TMP" "$VCF"
        mv "${TMP}.csi" "${VCF}.csi"

        echo "variants=$(bcftools index -n "$VCF")"
        date
    } > "$SAMPLE_LOG" 2>&1

    if ! bcftools index -n "$VCF" >/dev/null 2>&1; then
        echo "ERROR: invalid output VCF for $SAMPLE" >&2
        exit 1
    fi

    NVAR=$(bcftools index -n "$VCF")
    echo "[01] DONE: $SAMPLE variants=$NVAR"

done < "$SAMPLE_LIST"

echo "[01] Final audit"

printf "Expected samples\t%s\n" "$NSAMPLES"

printf "Completed VCFs\t"
find "$OUT" -maxdepth 1 \
  -name '*.genomewide_cds.hap1.vcf.gz' |
wc -l

printf "Completed indexes\t"
find "$OUT" -maxdepth 1 \
  -name '*.genomewide_cds.hap1.vcf.gz.csi' |
wc -l

echo "[01] Per-sample variant counts"

for VCF in "$OUT"/*.genomewide_cds.hap1.vcf.gz; do
    SAMPLE=$(basename "$VCF" .genomewide_cds.hap1.vcf.gz)
    printf "%s\t%s\n" \
      "$SAMPLE" \
      "$(bcftools index -n "$VCF")"
done |
sort -k2,2n > "$OUT/variant_counts.tsv"

cat "$OUT/variant_counts.tsv"

echo "[01] All samples completed"
