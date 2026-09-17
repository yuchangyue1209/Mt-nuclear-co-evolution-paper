#!/usr/bin/env bash
set -euo pipefail

ROOT="/path/to/data/genomewide_codeml_kuster"
PILOT="$ROOT/06_codeml/pilot"
LOG="$ROOT/logs/06_codeml_pilot.log"

if ! command -v codeml >/dev/null 2>&1; then
    echo "ERROR: codeml not found in PATH" >&2
    exit 1
fi

printf "gene_id\texit_status\tmlc_present\ttime_present\n" \
> "$PILOT/pilot_run_status.tsv"

tail -n +2 "$PILOT/pilot_manifest.tsv" |
while IFS=$'\t' read -r GENE CLASS VARIABLE LENGTH; do
    DIR="$PILOT/$GENE"

    echo "[pilot] $GENE class=$CLASS"

    rm -f \
        "$DIR/mlc" \
        "$DIR/rst" \
        "$DIR/rst1" \
        "$DIR/2NG.dN" \
        "$DIR/2NG.dS" \
        "$DIR/2NG.t" \
        "$DIR/lnf" \
        "$DIR/rub"

    set +e

    (
        cd "$DIR"
        codeml codeml.ctl
    ) > "$DIR/codeml.stdout.log" 2>&1

    STATUS=$?

    set -e

    MLC_PRESENT="no"
    TIME_PRESENT="no"

    [[ -s "$DIR/mlc" ]] && MLC_PRESENT="yes"

    if [[ -s "$DIR/mlc" ]] &&
       grep -q 'Time used' "$DIR/mlc"; then
        TIME_PRESENT="yes"
    fi

    printf "%s\t%s\t%s\t%s\n" \
        "$GENE" \
        "$STATUS" \
        "$MLC_PRESENT" \
        "$TIME_PRESENT" \
        >> "$PILOT/pilot_run_status.tsv"

    echo "[pilot done] $GENE status=$STATUS"
done

echo
echo "===== Pilot summary ====="

column -t -s $'\t' \
"$PILOT/pilot_run_status.tsv"

echo "[06 pilot] COMPLETE"
