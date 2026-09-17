#!/usr/bin/env bash
set -euo pipefail

# Converts all mitochondrial BAMs to a synchronized Pool-seq allele-count file.
# Requires Popoolation2 mpileup2sync.pl on PATH or supplied as argument 3.
source "${1:?Provide config/paths.sh}"
THREADS="${2:-8}"
MPILEUP2SYNC="${3:-mpileup2sync.pl}"
BAM_DIR="$DATA_ROOT/mt_bam"
OUT_DIR="$DATA_ROOT/mt_sync"
mkdir -p "$OUT_DIR"

find "$BAM_DIR" -maxdepth 1 -name '*.mt.sorted.bam' -print | sort > "$OUT_DIR/bam.list"
[[ -s "$OUT_DIR/bam.list" ]] || { echo "No mitochondrial BAMs found" >&2; exit 1; }

samtools mpileup -B -f "$MITO_REFERENCE_FASTA" -b "$OUT_DIR/bam.list" \
  -q 30 -Q 30 -d 5000 -o "$OUT_DIR/all_populations.mpileup"
perl "$MPILEUP2SYNC" --input "$OUT_DIR/all_populations.mpileup" \
  --output "$OUT_DIR/all_populations.sync" --fastq-type sanger --min-qual 20
gzip -f "$OUT_DIR/all_populations.mpileup"
