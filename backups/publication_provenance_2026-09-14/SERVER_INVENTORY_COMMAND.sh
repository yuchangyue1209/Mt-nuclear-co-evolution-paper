#!/usr/bin/env bash
set -euo pipefail

OUT="${1:-$PWD/mitonuclear_server_inventory.tsv}"

find /work/cyu /mnt/spareHD_2 -type f \
  \( -iname '*.R' -o -iname '*.Rmd' -o -iname '*.py' \
     -o -iname '*.sh' -o -iname '*.slurm' -o -iname '*.ini' \
     -o -iname '*.tsv' -o -iname '*.csv' -o -iname '*.txt' \
     -o -iname '*.yml' -o -iname '*.yaml' -o -iname '*.log' \
     -o -iname '*.pdf' -o -iname '*.png' \) \
  \( -ipath '*mito*' -o -ipath '*oxphos*' -o -ipath '*deltaAF*' \
     -o -ipath '*codeml*' -o -ipath '*pbs*' -o -ipath '*pinpis*' \
     -o -ipath '*erc*' -o -ipath '*cox4i1*' -o -ipath '*ndufs2*' \
     -o -ipath '*nu_287*' -o -ipath '*Kuster2026*' \) \
  -printf '%s\t%TY-%Tm-%TdT%TH:%TM:%TS\t%p\n' 2>/dev/null \
  | sort -k3,3 > "$OUT"

printf 'Wrote %s records to %s\n' "$(wc -l < "$OUT")" "$OUT"
