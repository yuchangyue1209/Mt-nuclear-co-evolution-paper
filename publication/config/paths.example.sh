#!/usr/bin/env bash

# Copy this file to config/paths.sh and replace the example paths.
# config/paths.sh is local configuration and should not be committed.

export PROJECT_ROOT="/path/to/stickleback-mitonuclear-project"
export DATA_ROOT="/path/to/large-analysis-data"
export RAW_READS_DIR="$DATA_ROOT/raw_reads"
export TRIMMED_READS_DIR="$DATA_ROOT/trimmed_reads"
export REFERENCE_FASTA="/path/to/GAculeatus_UGA_version5.fa"
export NUCLEAR_REFERENCE_FASTA="/path/to/GAculeatus_UGA_version5_nuclear.fa"
export MITO_REFERENCE_FASTA="/path/to/stickleback_mitochondrial_reference.fa"
export MITO_CONTIG="MH205729.1"
export REFERENCE_GFF="/path/to/GAculeatus_UGA_version5.gff3.gz"
export ADAPTER_FASTA="/path/to/bbmap/resources/adapters.fa"
export POPULATION_METADATA="$PROJECT_ROOT/metadata/populations.tsv"
export POOL_INFO="$PROJECT_ROOT/metadata/pool_info.tsv"
export RESULTS_ROOT="$PROJECT_ROOT/results"
export SIFT4G_BIN="sift4g"
