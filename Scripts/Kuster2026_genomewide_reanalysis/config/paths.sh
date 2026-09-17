#!/usr/bin/env bash

# Central server paths. Source this file; do not execute it directly.

export KUSTER_CODE_ROOT="/work/cyu/Kuster2026_genomewide_reanalysis"
export KUSTER_DATA_ROOT="/mnt/spareHD_2/genomewide_codeml_kuster"

export KUSTER_CLASSIFICATION="$KUSTER_CODE_ROOT/stickleback_mapping/stickleback_gene_classification.final.tsv"
export KUSTER_OLD_LABELS="$KUSTER_CODE_ROOT/old287_vs_genomewide/old287_three_labels_codeml20347.tsv"

export KUSTER_GTF="$KUSTER_DATA_ROOT/00_targets_codeml_ready/stickleback_20347_unique_codeml.gtf"
export KUSTER_ALIGNMENT_QA="$KUSTER_DATA_ROOT/05_gene_alignments/qa/gene_alignment_summary.with_class.tsv"
export KUSTER_CODEML_RESULTS="$KUSTER_DATA_ROOT/07_codeml_genomewide/codeml_results.tsv"
export KUSTER_CODEML_SITES="$KUSTER_DATA_ROOT/07_codeml_genomewide/codeml_site_information.tsv"
export KUSTER_MASTER="$KUSTER_DATA_ROOT/07_codeml_genomewide/codeml_master_analysis.tsv"

export KUSTER_OLD_MT_RESULTS="/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/codeml_sites_summary.merged.tsv"
export KUSTER_FIG2_OUT="$KUSTER_DATA_ROOT/09_figures/Figure2_complete_updated"

# Genome-wide standing-variation analysis.
export KUSTER_PINPIS_ROOT="$KUSTER_CODE_ROOT/genomewide_pinpis"
export KUSTER_SYNC="/work/cyu/nuclear_withNorway.merged.sync"
export KUSTER_REFERENCE="/work/cyu/stickleback_nuclear_only.fa"
export KUSTER_PINPIS_COLUMN_MAP="$KUSTER_PINPIS_ROOT/meta/sync_column_map.tsv"
export KUSTER_PINPIS_POPULATION="$KUSTER_PINPIS_ROOT/results/genomewide_20347_pinpis_population.tsv"
export KUSTER_PINPIS_GENE="$KUSTER_PINPIS_ROOT/results/genomewide_20347_pinpis_gene_level.tsv"
export KUSTER_PINPIS_MT13="$KUSTER_PINPIS_ROOT/results/mt13_pinpis_population.corrected.tsv"
export KUSTER_PINPIS_COMBINED="$KUSTER_PINPIS_ROOT/results/genomewide_nuclear_mt13_pinpis_combined.tsv"
export KUSTER_PINPIS_TESTS="$KUSTER_PINPIS_ROOT/results/pinpis_group_tests.tsv"

export KUSTER_EXPECT_CLASSIFIED=20426
export KUSTER_EXPECT_CODON_READY=20347
export KUSTER_EXPECT_PRIMARY=17965
export KUSTER_EXPECT_VARIABLE=15241
export KUSTER_EXPECT_INVARIANT=2724
export KUSTER_EXPECT_PINPIS_PRIMARY=18696
export KUSTER_EXPECT_PINPIS_DNDS_MATCHED=17918
export KUSTER_EXPECT_MT_PCG=13
