#!/usr/bin/env bash
set -euo pipefail

# ---- paths (edit as needed) ----
OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
QCOUT="$OUTDIR/cnv_qc_backfilled_nuclear"
COMB_PARTS="$QCOUT/parts.present+putative.4col.bed"   # 已合好的 4-col BED（confident + putative）
POOL_INFO="/work/cyu/poolseq/pool_info.txt"           # 第一列 SAMPLE，含表头
INPUT_DIR="/mnt/spareHD_2/nuclear_with_readgroups"    # ${SAMPLE}_rg.bam
REFERENCE="/work/cyu/stickleback_nuclear_only.fa"     # nuclear-only FASTA
THREADS=8

mkdir -p "$QCOUT" "$QCOUT/tmp"
export TMPDIR="$QCOUT/tmp"

# ---- sanity ----
for f in "$COMB_PARTS" "$POOL_INFO" "$REFERENCE"; do
  [ -s "$f" ] || { echo "❌ missing: $f"; exit 1; }
done
command -v mosdepth >/dev/null || { echo "❌ mosdepth 未安装"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools 未安装"; exit 1; }

echo "== SETTINGS =="
echo "COMB_PARTS : $COMB_PARTS"
echo "POOL_INFO  : $POOL_INFO"
echo "INPUT_DIR  : $INPUT_DIR"
echo "REFERENCE  : $REFERENCE"
echo "THREADS    : $THREADS"
echo

# ---- 1) run mosdepth for every pool on the combined parts ----
tail -n +2 "$POOL_INFO" | awk '{print $1}' | while read -r SAMPLE; do
  [ -n "${SAMPLE:-}" ] || continue
  BAM="${INPUT_DIR}/${SAMPLE}_rg.bam"
  if [ ! -s "$BAM" ]; then
    echo "⚠️  跳过 $SAMPLE（找不到 $BAM）"
    continue
  fi
  # 确保有 index
  samtools quickcheck "$BAM" 2>/dev/null || samtools index -@8 "$BAM" || true

  echo "-> mosdepth (COMB) on $SAMPLE"
  mosdepth \
    -t "$THREADS" \
    -n -x -Q 1 \
    --by "$COMB_PARTS" \
    -f "$REFERENCE" \
    "$QCOUT/${SAMPLE}.COMB" \
    "$BAM"

  ls -lh "$QCOUT/${SAMPLE}.COMB.regions.bed.gz" || true
done

echo "✅ mosdepth 完成：pattern = $QCOUT/<POOL>.COMB.regions.bed.gz"
echo "下一步：用 build_gene_depth_cn.py 汇总 depth → CN"
