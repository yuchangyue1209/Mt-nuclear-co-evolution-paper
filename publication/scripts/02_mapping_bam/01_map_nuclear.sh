#!/usr/bin/env bash
set -euo pipefail

# Maps trimmed Pool-seq reads to the nuclear-only reference. The historical
# workflow removed duplicate fragments before adding one read group per pool.
source "${1:?Provide config/paths.sh}"
THREADS="${2:-8}"
INDEX_PREFIX="$DATA_ROOT/reference_indices/nuclear"
SORTED_DIR="$DATA_ROOT/nuclear_bam/sorted"
DEDUP_DIR="$DATA_ROOT/nuclear_bam/deduplicated"
FINAL_DIR="$DATA_ROOT/nuclear_bam/final"
mkdir -p "$(dirname "$INDEX_PREFIX")" "$SORTED_DIR" "$DEDUP_DIR" "$FINAL_DIR"

[[ -f "${INDEX_PREFIX}.1.bt2" || -f "${INDEX_PREFIX}.1.bt2l" ]] || \
  bowtie2-build "$NUCLEAR_REFERENCE_FASTA" "$INDEX_PREFIX"

shopt -s nullglob
r1_files=("$TRIMMED_READS_DIR"/trimmed_R1_*.fastq.gz)
(( ${#r1_files[@]} > 0 )) || { echo "No trimmed reads found" >&2; exit 1; }

for r1 in "${r1_files[@]}"; do
  sample="$(basename "$r1" .fastq.gz)"; sample="${sample#trimmed_R1_}"
  r2="$TRIMMED_READS_DIR/trimmed_R2_${sample}.fastq.gz"
  [[ -f "$r2" ]] || { echo "Missing mate for $sample" >&2; continue; }
  sorted="$SORTED_DIR/${sample}.sorted.bam"
  dedup="$DEDUP_DIR/${sample}.dedup.bam"
  final="$FINAL_DIR/${sample}.bam"

  bowtie2 -x "$INDEX_PREFIX" -1 "$r1" -2 "$r2" -p "$THREADS" \
    --no-mixed --no-discordant -X 2000 2> "$SORTED_DIR/${sample}.bowtie2.log" |
    samtools view -@ "$THREADS" -b -q 20 - |
    samtools sort -@ "$THREADS" -o "$sorted" -

  picard MarkDuplicates I="$sorted" O="$dedup" \
    M="$DEDUP_DIR/${sample}.metrics.txt" REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

  picard AddOrReplaceReadGroups I="$dedup" O="$final" \
    RGID="$sample" RGLB="$sample" RGPL=ILLUMINA RGPU="${sample}.unit1" \
    RGSM="$sample" VALIDATION_STRINGENCY=LENIENT
  samtools index "$final"
done
