#!/usr/bin/env bash

# Central server paths. Source this file; do not execute it directly.

export GENOMEWIDE_CODE_ROOT="/path/to/workspace/genomewide_reanalysis"
export GENOMEWIDE_DATA_ROOT="/path/to/data/genomewide_codeml"

export GENOMEWIDE_CLASSIFICATION="$GENOMEWIDE_CODE_ROOT/stickleback_mapping/stickleback_gene_classification.final.tsv"
export GENOMEWIDE_OLD_LABELS="$GENOMEWIDE_CODE_ROOT/old287_vs_genomewide/old287_three_labels_codeml20347.tsv"

export GENOMEWIDE_GTF="$GENOMEWIDE_DATA_ROOT/00_targets_codeml_ready/stickleback_20347_unique_codeml.gtf"
export GENOMEWIDE_ALIGNMENT_QA="$GENOMEWIDE_DATA_ROOT/05_gene_alignments/qa/gene_alignment_summary.with_class.tsv"
export GENOMEWIDE_CODEML_RESULTS="$GENOMEWIDE_DATA_ROOT/07_codeml_genomewide/codeml_results.tsv"
export GENOMEWIDE_CODEML_SITES="$GENOMEWIDE_DATA_ROOT/07_codeml_genomewide/codeml_site_information.tsv"
export GENOMEWIDE_MASTER="$GENOMEWIDE_DATA_ROOT/07_codeml_genomewide/codeml_master_analysis.tsv"

export GENOMEWIDE_OLD_MT_RESULTS="/path/to/data/oxphos_codeml_ready/09_codeml_sites_models/codeml_sites_summary.merged.tsv"
export GENOMEWIDE_FIG2_OUT="$GENOMEWIDE_DATA_ROOT/09_figures/Figure2_complete_updated"

# Genome-wide standing-variation analysis.
export GENOMEWIDE_PINPIS_ROOT="$GENOMEWIDE_CODE_ROOT/genomewide_pinpis"
export GENOMEWIDE_SYNC="/path/to/workspace/nuclear_withNorway.merged.sync"
export GENOMEWIDE_REFERENCE="/path/to/workspace/stickleback_nuclear_only.fa"
export GENOMEWIDE_PINPIS_COLUMN_MAP="$GENOMEWIDE_PINPIS_ROOT/meta/sync_column_map.tsv"
export GENOMEWIDE_PINPIS_POPULATION="$GENOMEWIDE_PINPIS_ROOT/results/genomewide_20347_pinpis_population.tsv"
export GENOMEWIDE_PINPIS_GENE="$GENOMEWIDE_PINPIS_ROOT/results/genomewide_20347_pinpis_gene_level.tsv"
export GENOMEWIDE_PINPIS_MT13="$GENOMEWIDE_PINPIS_ROOT/results/mt13_pinpis_population.corrected.tsv"
export GENOMEWIDE_PINPIS_COMBINED="$GENOMEWIDE_PINPIS_ROOT/results/genomewide_nuclear_mt13_pinpis_combined.tsv"
export GENOMEWIDE_PINPIS_TESTS="$GENOMEWIDE_PINPIS_ROOT/results/pinpis_group_tests.tsv"

export GENOMEWIDE_EXPECT_CLASSIFIED=20426
export GENOMEWIDE_EXPECT_CODON_READY=20347
export GENOMEWIDE_EXPECT_PRIMARY=17965
export GENOMEWIDE_EXPECT_VARIABLE=15241
export GENOMEWIDE_EXPECT_INVARIANT=2724
export GENOMEWIDE_EXPECT_PINPIS_PRIMARY=18696
export GENOMEWIDE_EXPECT_PINPIS_DNDS_MATCHED=17918
export GENOMEWIDE_EXPECT_MT_PCG=13
