#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$HERE/config/paths.sh"

test "$(awk 'END{print NR-1}' "$KUSTER_ALIGNMENT_QA")" -eq "$KUSTER_EXPECT_CODON_READY"
test "$(awk 'END{print NR-1}' "$KUSTER_CODEML_RESULTS")" -eq "$KUSTER_EXPECT_VARIABLE"

Rscript "$HERE/R/09_build_codeml_master_table.R"
test "$(awk 'END{print NR-1}' "$KUSTER_MASTER")" -eq "$KUSTER_EXPECT_PRIMARY"

Rscript "$HERE/R/10_primary_codeml_tests.R"
Rscript "$HERE/R/11_plot_Figure2_complete.R"

