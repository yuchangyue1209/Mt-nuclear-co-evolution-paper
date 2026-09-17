#!/usr/bin/env bash
set -euo pipefail

LIST=${1:-SERVER_FILES_TO_RETRIEVE.txt}
OUT=${2:-/work/cyu/mitonuclear_publication_code_candidates.tar.gz}
TMP_LIST=$(mktemp)
trap 'rm -f "$TMP_LIST"' EXIT

while IFS= read -r path; do
    [[ -z "$path" || "$path" == \#* ]] && continue
    if [[ -f "$path" ]]; then
        printf '%s\n' "${path#/}" >> "$TMP_LIST"
    else
        printf '[MISSING] %s\n' "$path" >&2
    fi
done < "$LIST"

if [[ ! -s "$TMP_LIST" ]]; then
    printf '[STOP] No listed files were found.\n' >&2
    exit 1
fi

tar -C / -czf "$OUT" -T "$TMP_LIST"
printf '[OK] Bundle: %s\n' "$OUT"
printf '[OK] Files: %s\n' "$(wc -l < "$TMP_LIST" | tr -d ' ')"
sha256sum "$OUT"
