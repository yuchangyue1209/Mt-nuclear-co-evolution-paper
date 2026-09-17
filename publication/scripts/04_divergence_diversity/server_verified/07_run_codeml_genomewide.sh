#!/usr/bin/env bash
set -euo pipefail

ROOT="/path/to/data/genomewide_codeml_kuster"

QA="$ROOT/06_codeml/qa/codeml_input_QA.tsv"
OUT="$ROOT/07_codeml_genomewide"
RESULTS="$OUT/results"
RUN_LIST="$OUT/codeml_variable_genes.txt"
INVARIANT="$OUT/codeml_invariant_genes.tsv"

WORKER="/path/to/workspace/Kuster2026_genomewide_reanalysis/07_codeml_worker.sh"

JOBS=${JOBS:-8}

mkdir -p "$OUT" "$RESULTS"

if ! command -v codeml >/dev/null 2>&1; then
    echo "ERROR: codeml not found" >&2
    exit 1
fi

if [[ ! -s "$QA" ]]; then
    echo "ERROR: missing Step 06 QA: $QA" >&2
    exit 1
fi

awk -F'\t' '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="gene_id") gene=i
        if ($i=="codeml_input_status") status=i
    }
    next
}
$status=="variable" {print $gene}
' "$QA" > "$RUN_LIST"

awk -F'\t' '
NR==1 || $NF=="invariant"
' "$QA" > "$INVARIANT"

TOTAL=$(wc -l < "$RUN_LIST")
INVARIANT_N=$(awk 'END {print NR-1}' "$INVARIANT")

echo "===== Genome-wide codeml ====="
date
echo "Parallel jobs: $JOBS"
echo "Variable genes to run: $TOTAL"
echo "Invariant genes excluded: $INVARIANT_N"
echo "Results: $RESULTS"

df -h "$ROOT"

START=$(date +%s)

set +e

xargs -P "$JOBS" -n 1 \
"$WORKER" \
< "$RUN_LIST"

XARGS_STATUS=$?

set -e

END=$(date +%s)
ELAPSED=$((END - START))

SUCCESS=$(
    find "$RESULTS" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '.complete' |
    wc -l
)

FAILED=$(
    find "$RESULTS" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name 'failed.tsv' |
    wc -l
)

echo
echo "===== Final run audit ====="
echo "Target variable genes: $TOTAL"
echo "Completed genes: $SUCCESS"
echo "Failed records: $FAILED"
echo "Elapsed seconds: $ELAPSED"
echo "xargs status: $XARGS_STATUS"

SUMMARY="$OUT/codeml_results.tsv"
FAILURE_SUMMARY="$OUT/codeml_failures.tsv"

printf "gene_id\texit_status\telapsed_seconds\tlnL\tkappa\tomega\ttree_length_dN\ttree_length_dS\n" \
> "$SUMMARY"

while read -r GENE; do
    RESULT="$RESULTS/$GENE/result.tsv"

    if [[ -s "$RESULT" ]]; then
        tail -n 1 "$RESULT" >> "$SUMMARY"
    fi
done < "$RUN_LIST"

printf "gene_id\texit_status\telapsed_seconds\tstatus\n" \
> "$FAILURE_SUMMARY"

while read -r GENE; do
    FAILURE="$RESULTS/$GENE/failed.tsv"

    if [[ -s "$FAILURE" ]]; then
        tail -n 1 "$FAILURE" >> "$FAILURE_SUMMARY"
    fi
done < "$RUN_LIST"

RESULT_ROWS=$(awk 'END {print NR-1}' "$SUMMARY")
FAILURE_ROWS=$(awk 'END {print NR-1}' "$FAILURE_SUMMARY")

echo "Combined result rows: $RESULT_ROWS"
echo "Combined failure rows: $FAILURE_ROWS"
echo "Summary: $SUMMARY"
echo "Failures: $FAILURE_SUMMARY"

if [[ "$RESULT_ROWS" -eq "$TOTAL" ]]; then
    echo "[07] ALL VARIABLE GENES COMPLETE"
else
    echo "[07] INCOMPLETE: rerun the same command to resume"
fi

date
