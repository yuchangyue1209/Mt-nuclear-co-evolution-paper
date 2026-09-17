#!/usr/bin/env bash
set -euo pipefail

ROOT="/path/to/data/genomewide_codeml_kuster"
VCF_DIR="$ROOT/01_vcf_hap1"
CALL_ROOT="$ROOT/02_callability"
DEPTH_DIR="$CALL_ROOT/depth"
THRESHOLDS="$CALL_ROOT/depth_thresholds.tsv"

OUT="$ROOT/03_filtered_variants"
DEPTH_MASK="$OUT/depth_masks"
VARIANT_MASK="$OUT/uncertain_variant_masks"
COMBINED_MASK="$OUT/combined_masks"
PASS_VCF="$OUT/pass_vcf"
TMP_DIR="$OUT/tmp"
QA="$OUT/filter_QA.tsv"

ONLY_SAMPLE="${ONLY_SAMPLE:-}"

mkdir -p \
    "$DEPTH_MASK" \
    "$VARIANT_MASK" \
    "$COMBINED_MASK" \
    "$PASS_VCF" \
    "$TMP_DIR"

printf "sample\tmax_dp\tdepth_mask_bp\tuncertain_variant_mask_bp\tcombined_mask_bp\tpass_snps\n" \
> "$QA"

for DEPTH in "$DEPTH_DIR"/*.cds.depth.tsv.gz; do
    SAMPLE=$(basename "$DEPTH" .cds.depth.tsv.gz)

    if [[ -n "$ONLY_SAMPLE" && "$SAMPLE" != "$ONLY_SAMPLE" ]]; then
        continue
    fi

    echo "[03] START $SAMPLE"

    MAX_DP=$(
        awk -F'\t' -v sample="$SAMPLE" \
          'NR>1 && $1==sample {print $10}' \
          "$THRESHOLDS"
    )

    if [[ -z "$MAX_DP" ]]; then
        echo "ERROR: missing MAX_DP for $SAMPLE" >&2
        exit 1
    fi

    mapfile -t VCF_MATCHES < <(
        find "$VCF_DIR" -maxdepth 1 -type f \
          -name "${SAMPLE}*.vcf.gz" |
        sort
    )

    if [[ ${#VCF_MATCHES[@]} -ne 1 ]]; then
        echo "ERROR: $SAMPLE matched ${#VCF_MATCHES[@]} VCFs" >&2
        printf '%s\n' "${VCF_MATCHES[@]}" >&2
        exit 1
    fi

    VCF="${VCF_MATCHES[0]}"

    D_MASK="$DEPTH_MASK/${SAMPLE}.depth_mask.bed"
    V_MASK="$VARIANT_MASK/${SAMPLE}.uncertain_variant_mask.bed"
    C_MASK="$COMBINED_MASK/${SAMPLE}.combined_mask.bed"
    P_VCF="$PASS_VCF/${SAMPLE}.pass_snps.vcf.gz"

    PASS_KEYS="$TMP_DIR/${SAMPLE}.pass.keys.tsv"
    ALL_VARIANTS="$TMP_DIR/${SAMPLE}.all_variants.tsv"
    RAW_DEPTH_MASK="$TMP_DIR/${SAMPLE}.depth_mask.raw.bed"
    RAW_VARIANT_MASK="$TMP_DIR/${SAMPLE}.variant_mask.raw.bed"
    RAW_COMBINED="$TMP_DIR/${SAMPLE}.combined.raw.bed"

    echo "[03] $SAMPLE MAX_DP=$MAX_DP"

    # 1. High-confidence biallelic SNP VCF
    bcftools view \
      -v snps \
      -m2 -M2 \
      -i 'QUAL>=30 && FORMAT/DP>=10 && FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1])>=0.80' \
      -Oz \
      -o "${P_VCF}.tmp" \
      "$VCF"

    bcftools index -f "${P_VCF}.tmp"
    mv "${P_VCF}.tmp" "$P_VCF"
    mv "${P_VCF}.tmp.csi" "${P_VCF}.csi"

    # 2. Depth mask: DP<10 or DP>sample-specific MAX_DP
    zcat "$DEPTH" |
    awk -v maxdp="$MAX_DP" 'BEGIN{OFS="\t"} $3<10 || $3>maxdp {print $1,$2-1,$2}' \
    > "$RAW_DEPTH_MASK"

    sort -k1,1 -k2,2n -k3,3n "$RAW_DEPTH_MASK" |
    bedtools merge -i - \
    > "$D_MASK"

    # 3. Identify every original variant that failed the strict filter
    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\n' \
      "$P_VCF" \
      > "$PASS_KEYS"

    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\n' \
      "$VCF" \
      > "$ALL_VARIANTS"

    awk -F'\t' '
    BEGIN {OFS="\t"}
    NR==FNR {
        key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4
        pass[key]=1
        next
    }
    {
        key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4

        if (!(key in pass)) {
            start=$2-1
            end=start+length($3)

            if (end<=start)
                end=start+1

            print $1,start,end
        }
    }
    ' "$PASS_KEYS" "$ALL_VARIANTS" \
    > "$RAW_VARIANT_MASK"

    if [[ -s "$RAW_VARIANT_MASK" ]]; then
        sort -k1,1 -k2,2n -k3,3n "$RAW_VARIANT_MASK" |
        bedtools merge -i - \
        > "$V_MASK"
    else
        : > "$V_MASK"
    fi

    # 4. Union of depth and uncertain-variant masks
    {
        cat "$D_MASK"
        cat "$V_MASK"
    } |
    sort -k1,1 -k2,2n -k3,3n \
    > "$RAW_COMBINED"

    bedtools merge -i "$RAW_COMBINED" \
    > "$C_MASK"

    DEPTH_BP=$(awk '{sum+=$3-$2} END {print sum+0}' "$D_MASK")
    VARIANT_BP=$(awk '{sum+=$3-$2} END {print sum+0}' "$V_MASK")
    COMBINED_BP=$(awk '{sum+=$3-$2} END {print sum+0}' "$C_MASK")
    PASS_N=$(bcftools index -n "$P_VCF")

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$SAMPLE" "$MAX_DP" "$DEPTH_BP" "$VARIANT_BP" \
      "$COMBINED_BP" "$PASS_N" \
      >> "$QA"

    rm -f \
      "$PASS_KEYS" \
      "$ALL_VARIANTS" \
      "$RAW_DEPTH_MASK" \
      "$RAW_VARIANT_MASK" \
      "$RAW_COMBINED"

    echo "[03] DONE $SAMPLE pass_snps=$PASS_N combined_mask_bp=$COMBINED_BP"
done

echo "[03] COMPLETE"
column -t -s $'\t' "$QA"
