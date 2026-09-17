#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/spareHD_2/genomewide_codeml_kuster"

REF="/work/cyu/stickleback_nuclear_only.fa"
GTF="$ROOT/00_targets_codeml_ready/stickleback_20347_unique_codeml.gtf"

VCF_DIR="$ROOT/03_filtered_variants/pass_vcf"
MASK_DIR="$ROOT/03_filtered_variants/combined_masks"

OUT="$ROOT/04_consensus_cds"
GENOME_OUT="$OUT/masked_genomes"
CDS_OUT="$OUT/canonical_cds"
QA_OUT="$OUT/qa"
LOG_DIR="$ROOT/logs/04_consensus_cds"

EXPECTED_GENES=20347

mkdir -p \
    "$GENOME_OUT" \
    "$CDS_OUT" \
    "$QA_OUT" \
    "$LOG_DIR"

echo "===== 04 masked consensus and canonical CDS ====="
date
echo "[input] Reference: $REF"
echo "[input] GTF: $GTF"
echo "[input] VCF directory: $VCF_DIR"
echo "[input] Mask directory: $MASK_DIR"
echo "[output] $OUT"
echo "[expected] Canonical CDS per sample: $EXPECTED_GENES"

for PROGRAM in bcftools samtools gffread seqkit; do
    if ! command -v "$PROGRAM" >/dev/null 2>&1; then
        echo "ERROR: required program not found: $PROGRAM" >&2
        exit 1
    fi

    printf "[program] %-10s %s\n" \
        "$PROGRAM" \
        "$(command -v "$PROGRAM")"
done

for FILE in "$REF" "$GTF"; do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: missing or empty input: $FILE" >&2
        exit 1
    fi
done

if [[ ! -s "${REF}.fai" ]]; then
    echo "[index] Reference FASTA"
    samtools faidx "$REF"
fi

mapfile -t VCFS < <(
    find "$VCF_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*.pass_snps.vcf.gz' |
    sort
)

N_SAMPLES=${#VCFS[@]}

echo "[input] Filtered VCFs detected: $N_SAMPLES"

if [[ "$N_SAMPLES" -ne 27 ]]; then
    echo "ERROR: expected 27 filtered VCFs, found $N_SAMPLES" >&2
    exit 1
fi

SUMMARY="$QA_OUT/consensus_cds_summary.tsv"

printf "sample\tconsensus_contigs\tcds_records\tunique_cds_ids\tnontriplet_cds\tN_bases\tstatus\n" \
> "$SUMMARY"

for VCF in "${VCFS[@]}"; do
    BASE=$(basename "$VCF")
    SAMPLE=${BASE%.pass_snps.vcf.gz}

    MASK="$MASK_DIR/${SAMPLE}.combined_mask.bed"

    GENOME="$GENOME_OUT/${SAMPLE}.masked_consensus.fa"
    CDS="$CDS_OUT/${SAMPLE}.canonical_cds.fa"

    TMP_GENOME="${GENOME}.tmp"
    TMP_CDS="${CDS}.tmp"

    SAMPLE_LOG="$LOG_DIR/${SAMPLE}.log"
    LENGTH_QA="$QA_OUT/${SAMPLE}.cds_lengths.tsv"
    NONTRIPLET_QA="$QA_OUT/${SAMPLE}.nontriplet_cds.tsv"
    DUPLICATE_QA="$QA_OUT/${SAMPLE}.duplicate_cds_ids.tsv"

    echo
    echo "[sample] $SAMPLE"

    if [[ ! -s "$MASK" ]]; then
        echo "ERROR: missing or empty mask: $MASK" >&2
        exit 1
    fi

    VCF_SAMPLE=$(bcftools query -l "$VCF")

    if [[ "$VCF_SAMPLE" != "$SAMPLE" ]]; then
        echo "ERROR: VCF sample mismatch" >&2
        echo "File sample: $SAMPLE" >&2
        echo "VCF sample:  $VCF_SAMPLE" >&2
        exit 1
    fi

    # A completed CDS can be safely skipped after rechecking its dimensions.
    if [[ -s "$GENOME" && -s "$CDS" ]]; then
        EXISTING_CDS=$(grep -c '^>' "$CDS" || true)

        EXISTING_UNIQUE=$(
            grep '^>' "$CDS" |
            sed 's/^>//; s/[[:space:]].*$//' |
            sort -u |
            wc -l
        )

        EXISTING_NONTRIPLET=$(
            seqkit fx2tab -n -l "$CDS" |
            awk '$2 % 3 != 0 {n++} END {print n+0}'
        )

        if [[ "$EXISTING_CDS" -eq "$EXPECTED_GENES" &&
              "$EXISTING_UNIQUE" -eq "$EXPECTED_GENES" &&
              "$EXISTING_NONTRIPLET" -eq 0 ]]; then

            CONSENSUS_CONTIGS=$(grep -c '^>' "$GENOME" || true)

            N_BASES=$(
                grep -v '^>' "$GENOME" |
                tr -cd 'Nn' |
                wc -c
            )

            printf "%s\t%s\t%s\t%s\t%s\t%s\tcomplete_existing\n" \
                "$SAMPLE" \
                "$CONSENSUS_CONTIGS" \
                "$EXISTING_CDS" \
                "$EXISTING_UNIQUE" \
                "$EXISTING_NONTRIPLET" \
                "$N_BASES" \
                >> "$SUMMARY"

            echo "[skip complete] $SAMPLE"
            continue
        fi

        echo "[redo] Existing output failed completeness checks"
    fi

    rm -f \
        "$TMP_GENOME" \
        "${TMP_GENOME}.fai" \
        "$TMP_CDS"

    {
        echo "sample=$SAMPLE"
        echo "vcf=$VCF"
        echo "mask=$MASK"
        echo "reference=$REF"
        echo "gtf=$GTF"
        date

        echo "===== Build masked consensus ====="

        bcftools consensus \
            -f "$REF" \
            -m "$MASK" \
            -s "$SAMPLE" \
            "$VCF" \
            > "$TMP_GENOME"

        if [[ ! -s "$TMP_GENOME" ]]; then
            echo "ERROR: consensus FASTA is empty" >&2
            exit 1
        fi

        CONSENSUS_CONTIGS=$(grep -c '^>' "$TMP_GENOME" || true)

        if [[ "$CONSENSUS_CONTIGS" -eq 0 ]]; then
            echo "ERROR: consensus FASTA has no headers" >&2
            exit 1
        fi

        echo "Consensus contigs: $CONSENSUS_CONTIGS"

        samtools faidx "$TMP_GENOME"

        echo "===== Extract canonical CDS ====="

        gffread \
            "$GTF" \
            -g "$TMP_GENOME" \
            -x "$TMP_CDS"

        if [[ ! -s "$TMP_CDS" ]]; then
            echo "ERROR: extracted CDS FASTA is empty" >&2
            exit 1
        fi

        seqkit fx2tab -n -l "$TMP_CDS" \
            > "$LENGTH_QA"

        CDS_RECORDS=$(wc -l < "$LENGTH_QA")

        UNIQUE_CDS_IDS=$(
            cut -f1 "$LENGTH_QA" |
            sort -u |
            wc -l
        )

        awk '$2 % 3 != 0' "$LENGTH_QA" \
            > "$NONTRIPLET_QA"

        NONTRIPLET_CDS=$(wc -l < "$NONTRIPLET_QA")

        cut -f1 "$LENGTH_QA" |
        sort |
        uniq -d \
            > "$DUPLICATE_QA"

        DUPLICATE_CDS_IDS=$(wc -l < "$DUPLICATE_QA")

        echo "Extracted CDS: $CDS_RECORDS"
        echo "Unique CDS IDs: $UNIQUE_CDS_IDS"
        echo "Duplicate CDS IDs: $DUPLICATE_CDS_IDS"
        echo "Non-triplet CDS: $NONTRIPLET_CDS"

        if [[ "$CDS_RECORDS" -ne "$EXPECTED_GENES" ]]; then
            echo "ERROR: expected $EXPECTED_GENES CDS records; found $CDS_RECORDS" >&2
            exit 1
        fi

        if [[ "$UNIQUE_CDS_IDS" -ne "$EXPECTED_GENES" ]]; then
            echo "ERROR: expected $EXPECTED_GENES unique CDS IDs; found $UNIQUE_CDS_IDS" >&2
            exit 1
        fi

        if [[ "$DUPLICATE_CDS_IDS" -ne 0 ]]; then
            echo "ERROR: duplicated CDS IDs detected" >&2
            exit 1
        fi

        if [[ "$NONTRIPLET_CDS" -ne 0 ]]; then
            echo "ERROR: non-triplet CDS records detected" >&2
            exit 1
        fi

        N_BASES=$(
            grep -v '^>' "$TMP_GENOME" |
            tr -cd 'Nn' |
            wc -c
        )

        echo "Masked N bases: $N_BASES"

        mv "$TMP_GENOME" "$GENOME"
        mv "${TMP_GENOME}.fai" "${GENOME}.fai"
        mv "$TMP_CDS" "$CDS"

        printf "%s\t%s\t%s\t%s\t%s\t%s\tcomplete\n" \
            "$SAMPLE" \
            "$CONSENSUS_CONTIGS" \
            "$CDS_RECORDS" \
            "$UNIQUE_CDS_IDS" \
            "$NONTRIPLET_CDS" \
            "$N_BASES" \
            >> "$SUMMARY"

        echo "[done] $SAMPLE"
        date

    } > "$SAMPLE_LOG" 2>&1

    echo "[done] $SAMPLE"
done

sort -t $'\t' -k1,1 "$SUMMARY" \
    > "${SUMMARY}.sorted"

{
    head -1 "$SUMMARY"
    tail -n +2 "${SUMMARY}.sorted"
} > "${SUMMARY}.final"

mv "${SUMMARY}.final" "$SUMMARY"
rm -f "${SUMMARY}.sorted"

echo
echo "===== Final audit ====="

MASKED_GENOMES=$(
    find "$GENOME_OUT" \
        -maxdepth 1 \
        -type f \
        -name '*.masked_consensus.fa' |
    wc -l
)

CDS_FILES=$(
    find "$CDS_OUT" \
        -maxdepth 1 \
        -type f \
        -name '*.canonical_cds.fa' |
    wc -l
)

SUMMARY_ROWS=$(
    awk 'END {print NR-1}' "$SUMMARY"
)

echo "Expected samples: 27"
echo "Masked genomes: $MASKED_GENOMES"
echo "Canonical CDS files: $CDS_FILES"
echo "Summary rows: $SUMMARY_ROWS"

if [[ "$MASKED_GENOMES" -ne 27 ||
      "$CDS_FILES" -ne 27 ||
      "$SUMMARY_ROWS" -ne 27 ]]; then
    echo "ERROR: final sample count is incomplete" >&2
    exit 1
fi

echo
column -t -s $'\t' "$SUMMARY"

echo
echo "[04] All 27 samples completed"
echo "[output] $OUT"
date
