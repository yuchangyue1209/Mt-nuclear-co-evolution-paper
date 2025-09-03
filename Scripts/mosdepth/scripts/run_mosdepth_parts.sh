#!/usr/bin/env bash
set -euo pipefail

# Run mosdepth on present OXPHOS genes (nuclear), using exon/CDS parts,
# producing per-pool *.regions.bed.gz, and the 4-col BED used for --by.
#
# Requirements: bash, awk, coreutils, samtools, mosdepth
# Author: you
# License: MIT

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Required:
  -o, --outdir DIR           Working/output directory (has backfilled TSV) [default: /work/cyu/oxphos_from_ref_no_biomart/06_igv]
  -B, --backfilled FILE      Backfilled status table (TSV, has 'present_in_fish' & 'stickleback_name')
  -P, --parts FILE           6-column parts BED (prefer sorted); will fallback if missing
  -I, --pool-info FILE       pool_info.txt (header + SAMPLE [SIZE ...])
  -D, --bam-dir DIR          Directory of BAMs: \${SAMPLE}_rg.bam
  -r, --reference FASTA      Nuclear-only reference FASTA (for mosdepth -f)

Optional:
  -t, --threads INT          Threads for mosdepth (default: 8)
  -s, --test-sample ID       One sample to quick test mosdepth first (default: 1_FG)
  --qcsub DIR                QC subdir under outdir (default: cnv_qc_backfilled_nuclear)
  -h, --help                 Show this help

Outputs (under qcsub):
  - parts.present.uid.4col.bed     4-column BED used by mosdepth --by
  - \${SAMPLE}.regions.bed.gz      regions depth per pool/sample
EOF
}

# -------- defaults (match your current layout) --------
OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv"
BACK=""
PARTS6="$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.sorted.bed"
POOL_INFO="/work/cyu/poolseq/pool_info.txt"
BAM_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
REF_FASTA="/work/cyu/stickleback_nuclear_only.fa"
THREADS=8
TEST_SAMPLE="1_FG"
QCSUB="cnv_qc_backfilled_nuclear"

# -------- parse args --------
while (( "$#" )); do
  case "$1" in
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -B|--backfilled) BACK="$2"; shift 2;;
    -P|--parts) PARTS6="$2"; shift 2;;
    -I|--pool-info) POOL_INFO="$2"; shift 2;;
    -D|--bam-dir) BAM_DIR="$2"; shift 2;;
    -r|--reference) REF_FASTA="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -s|--test-sample) TEST_SAMPLE="$2"; shift 2;;
    --qcsub) QCSUB="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

# -------- resolve paths --------
mkdir -p "$OUTDIR"
QCOUT="$OUTDIR/$QCSUB"
mkdir -p "$QCOUT"

# fallback to unsorted parts if sorted is missing
if [ ! -s "$PARTS6" ]; then
  if [ -s "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.bed" ]; then
    PARTS6="$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.bed"
  fi
fi

# -------- checks --------
command -v mosdepth >/dev/null || { echo "❌ mosdepth not found"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not found"; exit 1; }

: "${BACK:?"❌ Please provide --backfilled TSV"}"
for f in "$BACK" "$PARTS6" "$POOL_INFO" "$REF_FASTA"; do
  [ -s "$f" ] || { echo "❌ Missing file: $f"; exit 1; }
done

echo "== SETTINGS =="
echo "OUTDIR     : $OUTDIR"
echo "QCOUT      : $QCOUT"
echo "BACKFILLED : $BACK"
echo "PARTS6     : $PARTS6"
echo "POOL_INFO  : $POOL_INFO"
echo "BAM_DIR    : $BAM_DIR"
echo "REF_FASTA  : $REF_FASTA"
echo "THREADS    : $THREADS"
echo "TEST_SAMPLE: $TEST_SAMPLE"
echo

# -------- 1) present gene list from backfilled table --------
PRESENT_LIST="$QCOUT/fish_present_all.list"
awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++){
    k=$i
    if(k=="present_in_fish") P=i
    else if(k=="stickleback_name") S=i
  }
  if(!P || !S){ print "ERR: columns present_in_fish/stickleback_name not found" > "/dev/stderr"; exit 2 }
  next
}
($P=="yes" && $S!=""){
  n=split($S,a,";")
  for(j=1;j<=n;j++){
    s=a[j]; gsub(/^ +| +$/,"",s)
    if(s!="") print tolower(s)
  }
}' "$BACK" | sort -u > "$PRESENT_LIST"

echo "present(all): $(wc -l < "$PRESENT_LIST")"
[ -s "$PRESENT_LIST" ] || { echo "❌ present list empty"; exit 1; }

# -------- 2) filter parts & build 4-col BED (name = gene|exonN) --------
SRC6="$QCOUT/parts.src.6col.bed"
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6}' "$PARTS6" > "$SRC6"

FIL6="$QCOUT/parts.present.6col.bed"
awk 'NR==FNR{keep[$1]=1; next}
     {k=tolower($4); if(k in keep) print}' \
  "$PRESENT_LIST" "$SRC6" > "$FIL6"

BED4="$QCOUT/parts.present.uid.4col.bed"
awk 'BEGIN{FS=OFS="\t"} {i[$4]++; print $1,$2,$3,$4"|exon"i[$4] }' "$FIL6" \
  | tr -d '\r' > "$BED4"

echo "4-col BED -> $BED4  (lines: $(wc -l < "$BED4"))"

# -------- 3) sanity checks --------
echo "== BED 4-col sanity (expect no output) =="
awk 'BEGIN{FS=OFS="\t"} ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/){print "BADINT", NR, $0}' "$BED4" | sed -n '1,10p' || true
awk 'BEGIN{FS=OFS="\t"} ($3<=$2){print "BADLEN", NR, $0}' "$BED4" | sed -n '1,10p' || true
sed -n '1,5p' "$BED4" | cat -A

# -------- 4) test one sample --------
BAM_TEST="$BAM_DIR/${TEST_SAMPLE}_rg.bam"
if [ -s "$BAM_TEST" ]; then
  samtools quickcheck "$BAM_TEST" 2>/dev/null || samtools index -@8 "$BAM_TEST"
  echo "-> mosdepth (test) on ${TEST_SAMPLE}"
  mosdepth -t "$THREADS" -n --by "$BED4" -f "$REF_FASTA" "$QCOUT/${TEST_SAMPLE}.TEST" "$BAM_TEST"
  ls -lh "$QCOUT/${TEST_SAMPLE}.TEST.regions.bed.gz" || true
else
  echo "⚠️ TEST sample $TEST_SAMPLE not found at $BAM_TEST, skip test"
fi

# -------- 5) run all samples --------
tail -n +2 "$POOL_INFO" | awk '{print $1}' | while read -r SAMPLE; do
  [ -n "${SAMPLE:-}" ] || continue
  BAM="$BAM_DIR/${SAMPLE}_rg.bam"
  if [ ! -s "$BAM" ]; then
    echo "⚠️ skip $SAMPLE: missing $BAM"
    continue
  fi
  samtools quickcheck "$BAM" 2>/dev/null || samtools index -@8 "$BAM"
  echo "-> mosdepth on ${SAMPLE}"
  mosdepth -t "$THREADS" -n --by "$BED4" -f "$REF_FASTA" "$QCOUT/${SAMPLE}" "$BAM"
done

echo "✅ mosdepth done. 4-col BED: $BED4"
