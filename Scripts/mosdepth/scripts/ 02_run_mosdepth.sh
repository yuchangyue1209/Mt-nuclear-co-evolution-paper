# 02_run_mosdepth.sh
#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
QCOUT="$OUTDIR/cnv_qc_backfilled_nuclear"
COMB_PARTS="$QCOUT/parts.present+putative+spanDIAMOND.4col.bed"
POOL_INFO="/work/cyu/poolseq/pool_info.txt"
INPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
REFERENCE="/work/cyu/stickleback_nuclear_only.fa"

for f in "$COMB_PARTS" "$POOL_INFO" "$REFERENCE"; do
  [ -s "$f" ] || { echo "❌ missing: $f"; exit 1; }
done
command -v mosdepth >/dev/null || { echo "❌ mosdepth not found"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not found"; exit 1; }

tail -n +2 "$POOL_INFO" | awk '{print $1}' | while read -r SAMPLE; do
  BAM="${INPUT_DIR}/${SAMPLE}_rg.bam"
  [ -s "$BAM" ] || { echo "⚠️ skip $SAMPLE (no BAM)"; continue; }
  samtools quickcheck "$BAM" 2>/dev/null || samtools index -@8 "$BAM"
  echo "-> mosdepth on $SAMPLE"
  mosdepth -n -x -Q 1 --by "$COMB_PARTS" -f "$REFERENCE" "$QCOUT/${SAMPLE}.COMB" "$BAM"
done

echo "✅ mosdepth done."
