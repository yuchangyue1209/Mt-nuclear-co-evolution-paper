#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_env.sh"
source "$(dirname "$0")/utils.sh"

SPAN_IN="${OUTDIR}/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus.bed"
HUMMISS="${OUTDIR}/human_missing_subunits.fa"   # you already built this FASTA

# 1) make BLAST DB once
makeblastdb -in "$GENOME" -dbtype nucl -parse_seqids -out "$DBROOT/stk"

# 2) tBLASTn, outfmt 6 (qseqid sseqid pident length qlen sstart send evalue bitscore)
TBL="${DBROOT}/missing_subunits.tblastn.tsv"
tblastn -query "$HUMMISS" -db "$DBROOT/stk" -out "$TBL" \
  -evalue 1e-3 -word_size 2 -seg yes -comp_based_stats F \
  -max_target_seqs 200 -max_hsps 200 \
  -outfmt '6 qseqid sseqid pident length qlen sstart send evalue bitscore'

# 3) cluster to SPANs (strict pass: qcov>=0.45, pident>=0.20)
SPAN_TBL="${OUTDIR}/tblastn_putative.SPAN.bed"
MAP_TBL="${OUTDIR}/tblastn_putative.map.tsv"
python "${PYSRC}/tblastn_to_spans.py" \
  --tblastn "$TBL" \
  --span_out "$SPAN_TBL" \
  --map_out "$MAP_TBL" \
  --gap 3000 --pad 200 --min_qcov 0.45 --min_pident 0.20

# 4) merge spans + maps, re-diff
SPAN_OUT="${OUTDIR}/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.plus_tblastn.bed"
sort -k1,1 -k2,2n "$SPAN_IN" "$SPAN_TBL" > "$SPAN_OUT"

MAP_MERGED="${OUTDIR}/ortholog_map.subunits.tryDIAMOND.tbl1.tsv"
cat "$MAP_TRY1" "${OUTDIR}/ortholog_map.from_diamond.v3.tsv" "$MAP_TBL" \
| awk 'BEGIN{FS=OFS="\t"} NR==1{print; next}
       {key=toupper($1) FS toupper($2); if(!seen[key]++){print}}' > "$MAP_MERGED"

OUT="${OUTDIR}/nuOXPHOS_vs_human.v3map_mix_DIA_TBL.nuclear_only"
python /work/cyu/oxphos_from_ref_no_biomart/diff_oxphos_sets_mapped.py \
  --span-bed "$SPAN_OUT" \
  --ortholog-map "$MAP_MERGED" \
  --out-prefix "$OUT" \
  --ignore-mt --keep-ab

summ_line "[tBLASTn round1]" "${OUT}.summary.tsv"  # → 82
