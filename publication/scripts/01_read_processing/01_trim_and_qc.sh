#!/usr/bin/env bash
set -euo pipefail

# Usage: bash 01_trim_and_qc.sh /path/to/config/paths.sh [threads]
source "${1:?Provide config/paths.sh}"
THREADS="${2:-8}"
QC_DIR="$DATA_ROOT/fastqc_trimmed"
mkdir -p "$TRIMMED_READS_DIR" "$QC_DIR"

shopt -s nullglob
r1_files=("$RAW_READS_DIR"/*_R1_001.fastq.gz)
(( ${#r1_files[@]} > 0 )) || { echo "No paired FASTQ files found" >&2; exit 1; }

for r1 in "${r1_files[@]}"; do
  r2="${r1/_R1_/_R2_}"
  [[ -f "$r2" ]] || { echo "Missing mate for $r1" >&2; continue; }
  sample="$(basename "$r1" _R1_001.fastq.gz)"
  out1="$TRIMMED_READS_DIR/trimmed_R1_${sample}.fastq.gz"
  out2="$TRIMMED_READS_DIR/trimmed_R2_${sample}.fastq.gz"

  bbduk.sh in1="$r1" in2="$r2" out1="$out1" out2="$out2" \
    ref="$ADAPTER_FASTA" ktrim=rl trimq=20 minlength=25 ftl=10 \
    tossbrokenreads=t threads="$THREADS" \
    > "$TRIMMED_READS_DIR/${sample}.bbduk.log" 2>&1
  fastqc "$out1" "$out2" --outdir="$QC_DIR" --threads="$THREADS"
done
