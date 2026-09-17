#!/usr/bin/env bash

# Actual server locations. The historical directory name is retained on the
# server for compatibility; the local archive uses analysis-based naming.
export ANALYSIS_CODE_ROOT="/work/cyu/Kuster2026_genomewide_reanalysis"
export DNDS_DATA_ROOT="/mnt/spareHD_2/genomewide_codeml_kuster"
export PINPIS_DATA_ROOT="$ANALYSIS_CODE_ROOT/genomewide_pinpis"

export GENE_CLASSIFICATION="$ANALYSIS_CODE_ROOT/stickleback_mapping/stickleback_gene_classification.final.tsv"
export FUNCTIONAL_LABELS="$ANALYSIS_CODE_ROOT/old287_vs_genomewide/old287_three_labels_codeml20347.tsv"
export CANONICAL_GTF="$DNDS_DATA_ROOT/00_targets_codeml_ready/stickleback_20347_unique_codeml.gtf"
export REFERENCE_FASTA="/work/cyu/stickleback_nuclear_only.fa"
export GENOMEWIDE_SYNC="/work/cyu/nuclear_withNorway.merged.sync"

export DNDS_MASTER="$DNDS_DATA_ROOT/07_codeml_genomewide/codeml_master_analysis.tsv"
export PINPIS_POPULATION="$PINPIS_DATA_ROOT/results/genomewide_20347_pinpis_population.tsv"
export PINPIS_GENE="$PINPIS_DATA_ROOT/results/genomewide_20347_pinpis_gene_level.tsv"
export PINPIS_MT13="$PINPIS_DATA_ROOT/results/mt13_pinpis_population.corrected.tsv"
export PINPIS_COMBINED="$PINPIS_DATA_ROOT/results/genomewide_nuclear_mt13_pinpis_combined.tsv"

export EXPECT_CLASSIFIED=20426
export EXPECT_CODON_READY=20347
export EXPECT_DNDS_PRIMARY=17965
export EXPECT_PINPIS_PRIMARY=18696
export EXPECT_MATCHED=17918
export EXPECT_MT_PCG=13
