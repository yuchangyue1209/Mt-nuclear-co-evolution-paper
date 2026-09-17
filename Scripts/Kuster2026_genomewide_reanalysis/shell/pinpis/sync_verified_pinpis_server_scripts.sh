#!/usr/bin/env bash
set -euo pipefail

# Run on the Mac. Override SERVER if needed.
SERVER="${SERVER:-cyu@144.92.58.154}"
REMOTE="/work/cyu/Kuster2026_genomewide_reanalysis/genomewide_pinpis"
LOCAL="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

mkdir -p "$LOCAL/python/pinpis" "$LOCAL/R/pinpis"

scp "$SERVER:$REMOTE/03_build_gene_level_pinpis.R" \
    "$LOCAL/R/pinpis/03_build_gene_level_pinpis.R"
scp "$SERVER:$REMOTE/05_build_combined_pinpis_table.R" \
    "$LOCAL/R/pinpis/05_build_combined_pinpis_table.R"
scp "$SERVER:$REMOTE/06_test_pinpis_groups.R" \
    "$LOCAL/R/pinpis/06_test_pinpis_groups.R"
scp "$SERVER:$REMOTE/07_plot_pinpis_Figure2_large_fonts.R" \
    "$LOCAL/R/pinpis/07_plot_pinpis_Figure2_large_fonts.R"

echo "Synced verified piN/piS scripts into: $LOCAL"
