#!/usr/bin/env bash
set -euo pipefail

# --- Paths you already use ---
export OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
export GENOME="/work/cyu/stickleback_nuclear_only.fa"
export GFF="/work/cyu/stickleback_v5.gff3.gz"
export FAI="${GENOME}.fai"                 # samtools faidx GENOME if missing

# span/parts from your v3 “CDS + exon fallback” build
export SPAN_BASE="${OUTDIR}/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.bed"

# (optional) more complete span that already includes previous tBLASTn adds
# export SPAN_PLUS="${OUTDIR}/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus_tblastn.bed"

# human vs fish set definitions & maps
export SET_DIFF="${OUTDIR}/nuOXPHOS_vs_human"
export MAP_TRY1="${OUTDIR}/ortholog_map.subunits.try1.tsv"

# DIAMOND / BLAST work on the big drive
export DBROOT="/mnt/spareHD_2/blastdb_stickleback"
mkdir -p "$DBROOT" "$DBROOT/tmp"
export TMPDIR="$DBROOT/tmp"

# helper scripts
export PYSRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/../python" && pwd)"
