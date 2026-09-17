#!/usr/bin/env bash
set -euo pipefail

# Independent mapping to the mitochondrial reference. Duplicate removal was
# not applied because mtDNA occurs at high biological copy number.
source "${1:?Provide config/paths.sh}"
THREADS="${2:-8}"
INDEX_PREFIX="$DATA_ROOT/reference_indices/mitochondrial"
OUT_DIR="$DATA_ROOT/mt_bam"
mkdir -p "$(dirname "$INDEX_PREFIX")" "$OUT_DIR"

[[ -f "${INDEX_PREFIX}.1.bt2" || -f "${INDEX_PREFIX}.1.bt2l" ]] || \
  bowtie2-build "$MITO_REFERENCE_FASTA" "$INDEX_PREFIX"

shopt -s nullglob
r1_files=("$TRIMMED_READS_DIR"/trimmed_R1_*.fastq.gz)
(( ${#r1_files[@]} > 0 )) || { echo "No trimmed reads found" >&2; exit 1; }

for r1 in "${r1_files[@]}"; do
  sample="$(basename "$r1" .fastq.gz)"; sample="${sample#trimmed_R1_}"
  r2="$TRIMMED_READS_DIR/trimmed_R2_${sample}.fastq.gz"
  [[ -f "$r2" ]] || { echo "Missing mate for $sample" >&2; continue; }
  bam="$OUT_DIR/${sample}.mt.sorted.bam"

  bowtie2 -x "$INDEX_PREFIX" -1 "$r1" -2 "$r2" -p "$THREADS" \
    --very-sensitive-local --no-mixed --no-discordant -X 2000 \
    2> "$OUT_DIR/${sample}.bowtie2.log" |
    samtools view -@ "$THREADS" -b -q 20 - |
    samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
  samtools depth "$bam" | gzip -c > "$OUT_DIR/${sample}.depth.tsv.gz"
done
