#!/usr/bin/env bash
set -euo pipefail

REF="/work/cyu/stickleback_nuclear_only.fa"
BED="/work/cyu/nuOXPHOS_genes_with_complex_core.v2.bed"   # 用新的 v2 BED
BAM_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
OUT_VCF="/mnt/spareHD_2/oxphos_codeml_ready/01_vcf_hap1"
mkdir -p "$OUT_VCF"

# 参考索引
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# 按参考顺序排序/合并 BED
BED_SORTED="${OUT_VCF}/nuOXPHOS_genes_with_complex_core.v2.sorted.bed"
bedtools sort -faidx "${REF}.fai" -i "$BED" | bedtools merge > "$BED_SORTED"

for BAM in "$BAM_DIR"/*_rg.bam; do
  SAMPLE=$(basename "$BAM" _rg.bam)
  echo "[00] $SAMPLE → VCF (diploid call; later mask→hap1)"

  # 若无 .bai，生成索引
  [ -f "${BAM}.bai" ] || samtools index -@ 4 "$BAM"

  # mpileup → call(默认二倍体) → norm → sort → index
  bcftools mpileup -Ou -f "$REF" -T "$BED_SORTED" -q 30 -Q 25 -a FORMAT/AD,FORMAT/DP "$BAM" \
  | bcftools call -mv -Ou \
  | bcftools norm -f "$REF" -m -any -Ou \
  | bcftools sort -Oz -o "$OUT_VCF/${SAMPLE}.oxphos.hap1.vcf.gz"

  bcftools index -f "$OUT_VCF/${SAMPLE}.oxphos.hap1.vcf.gz"
done

echo "[00] all done -> $OUT_VCF"

#01b_fix_overlaps.sh
#!/usr/bin/env bash
set -euo pipefail

REF="/work/cyu/stickleback_nuclear_only.fa"
IN_DIR="/mnt/spareHD_2/oxphos_codeml_ready/01_vcf_hap1"
OUT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/01_vcf_hap1_fixed"
mkdir -p "$OUT_DIR"

[ -f "${REF}.fai" ] || samtools faidx "$REF"

for VCF in "$IN_DIR"/*.oxphos.hap1.vcf.gz; do
  S=$(basename "$VCF" .oxphos.hap1.vcf.gz)
  echo "[fix] $S"
  TMP=$(mktemp -d)

  # 规范化：拆多等位、排序、去重复
  bcftools norm -f "$REF" -m -any -Ou "$VCF" \
  | bcftools sort -Ou \
  | bcftools norm -d all -Ob -o "$TMP/normalized.bcf"
  bcftools index -f "$TMP/normalized.bcf"

  # 拆成 indel / snp
  bcftools view -v indels "$TMP/normalized.bcf" -Oz -o "$TMP/indels.vcf.gz"
  bcftools view -v snps  "$TMP/normalized.bcf" -Oz -o "$TMP/snps.vcf.gz"
  bcftools index -f "$TMP/indels.vcf.gz" "$TMP/snps.vcf.gz"

  # 由 indel 推出影响区间（最长长度）
  bcftools query -f '%CHROM\t%POS0\t%REF\t%ALT\n' "$TMP/indels.vcf.gz" \
  | awk 'BEGIN{OFS="\t"}{
        n=split($4,a,","); L=1;
        for(i=1;i<=n;i++){d=length(a[i])-length($3); if(d<0)d=-d; if(d>L)L=d}
        print $1,$2,$2+L
      }' \
  | bedtools sort -faidx "${REF}.fai" | bedtools merge > "$TMP/indels.bed"

  # 从 SNP 中剔除与 indel 区重叠的位点（保持 header）
  bcftools view -T ^"$TMP/indels.bed" "$TMP/snps.vcf.gz" -Oz -o "$TMP/snps.clean.vcf.gz"
  bcftools index -f "$TMP/snps.clean.vcf.gz"

  # 合并回去并排序输出
  bcftools concat -a "$TMP/indels.vcf.gz" "$TMP/snps.clean.vcf.gz" -Ou \
  | bcftools sort -Oz -o "$OUT_DIR/${S}.oxphos.hap1.fixed.vcf.gz"

  bcftools index -f "$OUT_DIR/${S}.oxphos.hap1.fixed.vcf.gz"
  rm -rf "$TMP"
done

echo "[fix] all -> $OUT_DIR"

#01_build_masks.sh
#!/usr/bin/env bash
set -euo pipefail

# ---- paths ----
REF="/work/cyu/stickleback_nuclear_only.fa"
BED="/work/cyu/nuOXPHOS_genes_with_complex_core.v2.bed"   # 你现在用的 v2
BAM_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
VCF_DIR="/mnt/spareHD_2/oxphos_codeml_ready/01_vcf_hap1_fixed"  # ← 01b 的输出
MASK_DIR="/mnt/spareHD_2/oxphos_codeml_ready/02_masks"

# ---- thresholds (可按需微调) ----
MIN_DP=6          # depth <6 -> mask
AF_LO=0.3       
AF_HI=0.7
INDEL_PAD=3       # ±3 bp 围绕 indel 的掩蔽

mkdir -p "$MASK_DIR"

# 参考索引
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# 按参考顺序排序/合并 BED（稳妥起见）
BED_SORTED="${MASK_DIR}/nuOXPHOS_genes_with_complex_core.v2.sorted.bed"
bedtools sort -faidx "${REF}.fai" -i "$BED" | bedtools merge > "$BED_SORTED"

# 清理旧的 mask（可注释掉）
rm -f "$MASK_DIR"/*.hap1.mask.bed 2>/dev/null || true

for VCF in "$VCF_DIR"/*.vcf.gz; do
  [ -e "$VCF" ] || continue
  SAMPLE=$(basename "$VCF" .vcf.gz)
  SAMPLE=${SAMPLE%.oxphos.hap1}
  SAMPLE=${SAMPLE%.oxphos.hap1.fixed}

  # 选择匹配的 BAM（有 .sorted 优先）
  if [ -f "${BAM_DIR}/${SAMPLE}_rg.sorted.bam" ]; then
    BAM="${BAM_DIR}/${SAMPLE}_rg.sorted.bam"
  else
    BAM="${BAM_DIR}/${SAMPLE}_rg.bam"
  fi

  echo "[01] ${SAMPLE} → building hap1 mask"
  TMP_AF="$(mktemp)"; TMP_DP="$(mktemp)"; TMP_INDEL="$(mktemp)"
  trap 'rm -f "$TMP_AF" "$TMP_DP" "$TMP_INDEL"' EXIT

  # 1) VCF 中的等位深度 → 过滤 (AF在0.2~0.8 或 DP<MIN_DP)
  # 01b 后 VCF 仍带 FORMAT/AD（若没有可改为用 INFO/DP4 或直接依赖 samtools depth）
  bcftools query -f '%CHROM\t%POS0\t%POS\t[%AD]\n' "$VCF" \
  | awk -v OFS='\t' -v mdp="$MIN_DP" -v lo="$AF_LO" -v hi="$AF_HI" '
      {
        split($4,a,","); dp=a[1]+a[2];
        if (dp==0) next;
        if (dp<mdp) {print $1,$2,$3; next}
        af=a[2]/dp;
        if (af>lo && af<hi) print $1,$2,$3
      }' > "$TMP_AF"

  # 2) 目标区间内的测序深度（即使不是变异位点）
  samtools depth -b "$BED_SORTED" "$BAM" \
  | awk -v mdp="$MIN_DP" 'BEGIN{OFS="\t"} $3<mdp {print $1,$2-1,$2}' > "$TMP_DP"

  # 3) indel 周边 ±INDEL_PAD
  bcftools view -v indels -H "$VCF" \
  | awk -v p="$INDEL_PAD" 'BEGIN{OFS="\t"}{ s=$2-1-p; if(s<0)s=0; e=$2+p; print $1,s,e }' > "$TMP_INDEL"

  # 4) 合并三类掩蔽区
  cat "$TMP_AF" "$TMP_DP" "$TMP_INDEL" \
  | bedtools sort -faidx "${REF}.fai" \
  | bedtools merge > "$MASK_DIR/${SAMPLE}.hap1.mask.bed"

  rm -f "$TMP_AF" "$TMP_DP" "$TMP_INDEL"
  trap - EXIT
  echo "[01] ${SAMPLE} ✓ mask -> $MASK_DIR/${SAMPLE}.hap1.mask.bed"
done

echo "[01] All masks done -> $MASK_DIR"


# 02_make_consensus_hap1.sh  —— 用 01b/fixed 的 VCF 生成共识（含轻清洗）
#!/usr/bin/env bash
set -euo pipefail
set -o pipefail

REF="/work/cyu/stickleback_nuclear_only.fa"
VCF_DIR="/mnt/spareHD_2/oxphos_codeml_ready/01_vcf_hap1_fixed"
MASK_DIR="/mnt/spareHD_2/oxphos_codeml_ready/02_masks"
CONS_DIR="/mnt/spareHD_2/oxphos_codeml_ready/03_consensus_hap1_snpOnly"
LOG_DIR="${CONS_DIR}/_logs"
TMP_DIR="${CONS_DIR}/_tmp"
mkdir -p "$CONS_DIR" "$LOG_DIR" "$TMP_DIR"
[ -f "${REF}.fai" ] || samtools faidx "$REF"

shopt -s nullglob
for VCF in "$VCF_DIR"/*.vcf.gz; do
  base=$(basename "$VCF" .vcf.gz)
  SAMPLE="${base%.oxphos.hap1.fixed}"; SAMPLE="${SAMPLE%.oxphos.hap1}"
  [[ -z "$SAMPLE" ]] && { echo "[02] bad sample from $VCF"; continue; }

  OUT="$CONS_DIR/${SAMPLE}.hap1.consensus.fa"
  LOG="$LOG_DIR/${SAMPLE}.log"
  MASK="$MASK_DIR/${SAMPLE}.hap1.mask.bed"
  CLEAN_BCF="$TMP_DIR/${SAMPLE}.clean.bcf"
  SNP_BCF="$TMP_DIR/${SAMPLE}.snps.bcf"     # ← 新增：把 SNP-only 落盘
  echo "[02] $SAMPLE → clean & SNP-only consensus"

  {
    set -x
    # 轻清洗到 BCF
    bcftools norm -f "$REF" -m -any -Ou "$VCF" \
    | bcftools norm -f "$REF" -c x -Ou \
    | bcftools norm -d all -Ou \
    | bcftools sort -Ob -o "$CLEAN_BCF"
    bcftools index -f "$CLEAN_BCF"

    # 只保留 SNP，写成 BCF 文件并建索引（避免用管道传给 consensus）
    bcftools view -v snps -Ob -o "$SNP_BCF" "$CLEAN_BCF"
    bcftools index -f "$SNP_BCF"

    # 跑 consensus：输入用文件，不用“-”
    if [[ -s "$MASK" ]]; then
      bcftools consensus -f "$REF" -m "$MASK" -M N "$SNP_BCF" > "$OUT"
    else
      echo "[02] WARN: mask missing -> $MASK (proceed without -m)"
      bcftools consensus -f "$REF" -M N "$SNP_BCF" > "$OUT"
    fi
    set +x

    # 简单 sanity check：非空且首行是 '>'
    if [[ ! -s "$OUT" ]] || ! head -n1 "$OUT" | grep -q '^>'; then
      echo "[02] ERROR: empty/invalid consensus: $OUT"
      exit 2
    fi

    samtools faidx "$OUT"
    echo "[02] ${SAMPLE} ✓ -> $OUT"
  } >"$LOG" 2>&1 || {
    echo "[02] ${SAMPLE} FAILED, see $LOG"
    rm -f "$OUT" "$OUT.fai"
    continue
  }

done
shopt -u nullglob


#!/usr/bin/env bash
set -euo pipefail

# ======= inputs (按需调整 TXLIST 到“本次 287 个 tx 的清单”) =======
GFF3="/work/cyu/stickleback_v5.gff3"    # 可为 .gff3 或 .gff3.gz
TXLIST="/work/cyu/oxphos_from_ref_no_biomart/06_igv/dnds_annotations_final/chosen.tx.final.list"  # 需为 287 行

# ======= outputs =======
ROOT="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep"
OUT="${ROOT}/gff3_fix_from_tx"
SUB_GFF3="${OUT}/nuOXPHOS_subset.gff3"
CANON_GTF="${ROOT}/nuOXPHOS_canonical.gtf"
ID2NAME="${ROOT}/id2name.tsv"
TXMAP="${ROOT}/tx_gene_symbol.tsv"
mkdir -p "$OUT"

# auto reader for (gz|plain)
if [[ "${GFF3##*.}" == "gz" ]]; then
  READER="zcat"
else
  READER="cat"
fi

echo "[03a] input GFF3: $GFF3"
echo "[03a] tx list    : $TXLIST"
echo "[03a] out dir    : $OUT"

# 0) 读入目标转录本集合（应为 287；若不是，仅警告继续）
TX_N=$(grep -vc '^\s*$' "$TXLIST" || true)
echo "[03a] transcripts in list: ${TX_N}"
if [[ "$TX_N" -ne 287 ]]; then
  echo "[03a][WARN] TXLIST 行数 ${TX_N} != 287，请确认你提供的是当前 287 的清单。"
fi

# 1) 建 transcript->gene 映射（mRNA/transcript 的 ID 与 Parent）
T2G="${OUT}/t2g.tsv"
$READER "$GFF3" \
| awk -F'\t' '
  $0!~/^#/ && ($3=="mRNA"||$3=="transcript"){
    tid=""; gid="";
    if (match($9,/ID=([^;]+)/,I))   tid=I[1];
    if (match($9,/Parent=([^;]+)/,P)) gid=P[1];
    if (tid!="" && gid!="") print tid"\t"gid
  }' > "$T2G"
echo "[03a] built t2g: $(wc -l < "$T2G") rows"

# 2) 用 tx 清单推断需要的 gene 列表
KEEP_G="${OUT}/keep_gene_ids.txt"
awk 'NR==FNR{want[$1]=1; next} ($1 in want){print $2}' "$TXLIST" "$T2G" \
  | sort -u > "$KEEP_G"
echo "[03a] genes inferred from txlist: $(wc -l < "$KEEP_G")"

# 3) 子集 GFF3：保留 (gene in keep) + (transcript in list) + (对应 exon/CDS/UTR)
$READER "$GFF3" \
| awk -F'\t' -v LTX="$TXLIST" -v KG="$KEEP_G" '
  BEGIN{
    while((getline t<LTX)>0){if(t!="") keepTx[t]=1}
    while((getline g<KG)>0){keepG[g]=1}
  }
  $0~/^#/ {next}
  {
    feat=$3; attr=$9; id=""; par="";
    if (match(attr,/ID=([^;]+)/,I)) id=I[1];
    if (match(attr,/Parent=([^;]+)/,P)) par=P[1];

    keep=0;
    if (feat=="gene" && keepG[id]) keep=1;
    else if ((feat=="mRNA"||feat=="transcript") && keepTx[id]) keep=1;
    else if ((feat=="exon"||feat=="CDS"||feat=="five_prime_UTR"||feat=="three_prime_UTR") && keepTx[par]) keep=1;

    if (keep) print
  }' > "$SUB_GFF3"
echo "[03a] subset gff3 lines: $(wc -l < "$SUB_GFF3") -> $SUB_GFF3"

# 4) GFF3 -> GTF（gffread 会写出 gene_id / transcript_id）
gffread -T -o "$CANON_GTF" "$SUB_GFF3"
echo "[03a] wrote canonical GTF: $CANON_GTF"

# 5) 导出 gene_id ↔ gene_symbol（从 gene 的 Name=；小写；若缺则用 gene_id 代替）
$READER "$GFF3" \
| awk -F"\t" '
  $0!~/^#/ && $3=="gene" {
    gid=""; nm="";
    if (match($9,/ID=([^;]+)/,I))   gid=I[1];
    if (match($9,/Name=([^;]+)/,N)) nm=tolower(N[1]);
    if (gid!="") { if(nm=="") nm=gid; print gid"\t"nm }
  }' > "$ID2NAME"
echo "[03a] wrote: $ID2NAME"

# 6) 生成 transcript_id ↔ gene_id ↔ symbol 三列表（供 03b/04 使用）
awk -F'\t' '$3=="transcript"{
  if (match($9,/transcript_id "([^"]+)"/,t) && match($9,/gene_id "([^"]+)"/,g))
    print t[1]"\t"g[1]
}' "$CANON_GTF" \
| sort -u \
| awk 'NR==FNR{s[$1]=$2; next}{gid=$2; sym=(gid in s?s[gid]:gid); print $1"\t"$2"\t"sym}' \
     "$ID2NAME" - > "$TXMAP"
echo "[03a] wrote: $TXMAP"

# 7) 快速 QC
echo -n "[QC] transcripts kept (unique): "
awk -F'\t' '$3=="transcript"{if (match($9,/transcript_id "([^"]+)"/,m)) print m[1] }' "$CANON_GTF" | sort -u | wc -l
echo -n "[QC] CDS lines: "
awk -F'\t' '$3=="CDS"{c++} END{print c+0}' "$CANON_GTF"

echo "[03a] DONE (for 287 genes)."


#03b 287gene
#!/usr/bin/env bash
set -euo pipefail

# === Inputs ===
GTF="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/nuOXPHOS_canonical.gtf"
# 若你这次的共识在 _snpOnly 目录，请保持这一行；否则改成 03_consensus_hap1
CONS_DIR="/mnt/spareHD_2/oxphos_codeml_ready/03_consensus_hap1_snpOnly"

# === Outputs ===
CDS_DIR="/mnt/spareHD_2/oxphos_codeml_ready/04_cds_by_sample"
PEP_DIR="/mnt/spareHD_2/oxphos_codeml_ready/05_pep_by_sample"
QC_DIR="/mnt/spareHD_2/oxphos_codeml_ready/00_qc"
mkdir -p "$CDS_DIR" "$PEP_DIR" "$QC_DIR"

[ -s "$GTF" ] || { echo "ERROR: GTF not found -> $GTF"; exit 1; }
[ -d "$CONS_DIR" ] || { echo "ERROR: consensus dir not found -> $CONS_DIR"; exit 1; }

echo "[03b] extracting using GTF: $GTF"
echo "[03b] from consensus dir  : $CONS_DIR"

# 结果汇总（覆盖旧文件）
PER_SAMPLE="$QC_DIR/03b_per_sample_qc.tsv"
P_STOP="$QC_DIR/pep_internal_stop.tsv"   # cols: sample  tx_id
P_LEN3="$QC_DIR/cds_len_mod3.tsv"        # cols: sample  tx_id
: > "$PER_SAMPLE"; : > "$P_STOP"; : > "$P_LEN3"

shopt -s nullglob
for F in "$CONS_DIR"/*.hap1.consensus.fa; do
  S=$(basename "$F" .hap1.consensus.fa)
  echo "[03b] $S → CDS/PEP"

  gffread -g "$F" -S "$GTF" \
          -x "$CDS_DIR/${S}.hap1.cds.fa" \
          -y "$PEP_DIR/${S}.hap1.pep.fa"

  samtools faidx "$CDS_DIR/${S}.hap1.cds.fa" >/dev/null 2>&1 || true
  samtools faidx "$PEP_DIR/${S}.hap1.pep.fa" >/dev/null 2>&1 || true

  # 逐转录本检查并记录问题清单
  # 1) CDS 长度%3!=0 的 tx
  awk -v samp="$S" 'BEGIN{RS=">";FS="\n"}
    NR>1{
      tx=$1; gsub(/\r/,"",tx);
      seq=""; for(i=2;i<=NF;i++) seq=seq $i;
      if (length(seq)%3!=0) print samp"\t"tx
    }' "$CDS_DIR/${S}.hap1.cds.fa" >> "$P_LEN3"

  # 2) 蛋白内部 *（忽略末尾终止*）
  awk -v samp="$S" 'BEGIN{RS=">";FS="\n"}
    NR>1{
      tx=$1; gsub(/\r/,"",tx);
      p=""; for(i=2;i<=NF;i++) p=p $i;
      sub(/\*$/,"",p);
      if (index(p,"*")) print samp"\t"tx
    }' "$PEP_DIR/${S}.hap1.pep.fa" >> "$P_STOP"

  # 统计数量
  cds_bad=$(awk 'END{print NR+0}' "$P_LEN3") || cds_bad=0
  pep_bad=$(awk 'END{print NR+0}' "$P_STOP") || pep_bad=0

  # 仅统计本样本的条数
  cds_bad_s=$(awk -v s="$S" '$1==s' "$P_LEN3" | wc -l)
  pep_bad_s=$(awk -v s="$S" '$1==s' "$P_STOP" | wc -l)

  printf "%s\t%s\t%s\n" "$S" "$cds_bad_s" "$pep_bad_s" >> "$PER_SAMPLE"
  if (( cds_bad_s>0 || pep_bad_s>0 )); then
    echo "[03b][WARN] $S : CDS_len%3!=0=${cds_bad_s}  internal*=${pep_bad_s}"
  fi
done
shopt -u nullglob

# 加表头（如无）
if [[ -s "$PER_SAMPLE" && "$(head -n1 "$PER_SAMPLE" | cut -f1)" != "sample" ]]; then
  sed -i '1isample\tcds_len_mod3_neq0\tpep_contains_internal_stop' "$PER_SAMPLE"
fi
if [[ -s "$P_STOP" && "$(head -n1 "$P_STOP" | awk "{print \$1}")" != "sample" ]]; then
  sed -i '1isample\ttx_id' "$P_STOP"
fi
if [[ -s "$P_LEN3" && "$(head -n1 "$P_LEN3" | awk "{print \$1}")" != "sample" ]]; then
  sed -i '1isample\ttx_id' "$P_LEN3"
fi

echo "[03b] done."
echo "  CDS: $CDS_DIR"
echo "  PEP: $PEP_DIR"
echo "  QC : $PER_SAMPLE"
echo "  QC : $P_STOP"
echo "  QC : $P_LEN3"


#04
#!/usr/bin/env bash
set -euo pipefail

# === Inputs ===
GTF="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/nuOXPHOS_canonical.gtf"
TXMAP="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/tx_gene_symbol.tsv"   # 3列: tx_id  gene_id  gene_symbol(小写)
CDS_DIR="/mnt/spareHD_2/oxphos_codeml_ready/04_cds_by_sample"
PEP_DIR="/mnt/spareHD_2/oxphos_codeml_ready/05_pep_by_sample"

# === Outputs ===
ALIGN_DIR="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
GENE_LIST="$ALIGN_DIR/genes72.list"   # 以此清单为准（应为 287）
mkdir -p "$ALIGN_DIR"

# 若你已有 287 的 list，可直接复用；否则从 TXMAP 导出（symbol 唯一）
if [[ ! -s "$GENE_LIST" ]]; then
  cut -f3 "$TXMAP" | sort -u > "$GENE_LIST"
fi

# 小工具：精确提取指定 header 的序列（FASTA）
extract_one () {
  local fa="$1" id="$2"
  awk -v id="$id" 'BEGIN{RS=">";FS="\n"}
    NR>1{hdr=$1; if(hdr==id){printf(">%s\n",hdr); for(i=2;i<=NF;i++)print $i; exit}}' "$fa"
}

# 过滤函数：返回 0 表示通过；1 表示该样本对该转录本应跳过（frame 或内部*）
ok_this_sample_tx () {
  local cds="$1" pep="$2" tx="$3"
  # CDS %3
  awk -v id="$tx" 'BEGIN{RS=">";FS="\n"}
    NR>1{
      h=$1; s="";
      for(i=2;i<=NF;i++) s=s $i;
      if(h==id){ exit (length(s)%3!=0 ? 1 : 0) }
    }' "$cds" || return 1
  # PEP 内部 *
  awk -v id="$tx" 'BEGIN{RS=">";FS="\n"}
    NR>1{
      h=$1; p="";
      for(i=2;i<=NF;i++) p=p $i;
      sub(/\*$/,"",p);
      if(h==id){ exit (index(p,"*") ? 1 : 0) }
    }' "$pep" || return 1
  return 0
}

echo "[04] building per-gene alignments from list: $GENE_LIST"
total=$(wc -l < "$GENE_LIST" 2>/dev/null || echo 0)
echo "[04] genes planned: $total"

while read -r G; do
  [[ -z "$G" ]] && continue

  # 找该基因的代表转录本（TXMAP：第3列=gene_symbol）
  TX=$(awk -v g="$G" '$3==g{print $1; exit}' "$TXMAP") || true
  if [[ -z "$TX" ]]; then
    echo "[skip] $G : no tx in TXMAP"
    continue
  fi

  outd="$ALIGN_DIR/$G"
  mkdir -p "$outd"
  pep_all="$outd/$G.pep.faa"
  cds_all="$outd/$G.cds.fna"
  : > "$pep_all"; : > "$cds_all"

  ok=0
  shopt -s nullglob
  for PEP in "$PEP_DIR"/*.hap1.pep.fa; do
    S=$(basename "$PEP" .hap1.pep.fa)
    CDS="$CDS_DIR/${S}.hap1.cds.fa"
    [[ -s "$PEP" && -s "$CDS" ]] || continue

    pseq=$(extract_one "$PEP" "$TX") || true
    cseq=$(extract_one "$CDS" "$TX") || true
    [[ -z "$pseq" || -z "$cseq" ]] && continue

    # 过滤掉 frame/内部* 的样本条目
    if ! ok_this_sample_tx "$CDS" "$PEP" "$TX"; then
      # echo "[skip-one] $G $S (frame/stop)"
      continue
    fi

    # header 改成样本名
    echo "$pseq" | sed "1 s/^>.*/>${S}/" >> "$pep_all"
    echo "$cseq" | sed "1 s/^>.*/>${S}/" >> "$cds_all"
    ok=$((ok+1))
  done
  shopt -u nullglob

  # codeml 至少需要 3 个 taxa 才稳妥；少于3就跳过
  if (( ok < 3 )); then
    echo "[skip] $G : usable_seqs=$ok (<3)"
    rm -f "$pep_all" "$cds_all"
    rmdir "$outd" 2>/dev/null || true
    continue
  fi

  # 蛋白对齐 → PAL2NAL
  mafft --maxiterate 1000 --localpair "$pep_all" > "$outd/$G.pep.aln.faa"
  pal2nal.pl "$outd/$G.pep.aln.faa" "$cds_all" -output fasta -nogap > "$outd/$G.codon.fas"

  # 快检
  nseq=$(grep -c '^>' "$outd/$G.codon.fas" || echo 0)
  if (( nseq < 3 )); then
    echo "[warn] $G : pal2nal output seqs=$nseq (<3) → remove"
    rm -f "$outd/$G.codon.fas" "$outd/$G.pep.aln.faa" "$pep_all" "$cds_all"
    rmdir "$outd" 2>/dev/null || true
    continue
  fi

  echo "[04] $G ✓ ($ok seqs)"
done < "$GENE_LIST"

echo "[04] all done -> $ALIGN_DIR"



#codmel
#!/usr/bin/env bash
set -euo pipefail

# ================= 路径 =================
TREEFILE="/mnt/spareHD_2/oxphos_gene_tree/unrooted_NJ_tree_pruned.standard.tree"
NU_GENE_DIR="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"

# 输出根目录（与 mt 版本保持相同布局）
OUTROOT="/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/nu"
mkdir -p "$OUTROOT/gene"
LOGERR="$OUTROOT/error_runs.log"
: > "$LOGERR"

# 汇总表（追加写入）
SUMMARY="$OUTROOT/codeml_sites_summary.tsv"
echo -e "gene\tmodel\tlnL\tkappa\tomega" > "$SUMMARY"

# ================= 模型列表（固定顺序） =================
declare -A NSMAP=([M0]=0 [M1a]=1 [M2a]=2 [M7]=7 [M8]=8)
MODELS=(M0 M1a M2a M7 M8)

# 依赖检查
command -v nw_prune >/dev/null 2>&1 || { echo "[FATAL] nw_prune not found in PATH"; exit 1; }
command -v codeml   >/dev/null 2>&1 || { echo "[FATAL] codeml not found in PATH"; exit 1; }

# ================= 小工具 =================
names_from_fa(){ grep '^>' "$1" | sed 's/^>//' | cut -d' ' -f1 | tr -d '\r'; }

write_ctl(){  # $1 seqfile $2 treefile $3 outprefix $4 NSsites
  cat > "$3.ctl" <<EOF
seqfile = $1
treefile = $2
outfile = $3.txt

runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 0
NSsites = $4
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
EOF
}

extract_stats(){ # $1 result.txt -> print "lnL\tkappa\tomega"
  local f="$1" lnL kappa omega
  lnL=$(awk   '/lnL/ && /time/ {for(i=1;i<=NF;i++) if($i=="lnL"){print $(i+2); exit}}' "$f" || true)
  kappa=$(awk '/kappa \(ts\/tv\)/ {print $4; exit}' "$f" || true)
  omega=$(awk '/omega \(dN\/dS\)/ {print $4; exit}' "$f" || true)
  echo -e "${lnL:-}\t${kappa:-}\t${omega:-}"
}

# ================= 运行函数 =================
run_codeml() {
  local align=$1     # 多序列密码子比对（绝对或相对路径都可）
  local tag=$2       # 基因名
  local outdir=$3    # 输出目录: $OUTROOT/gene/<gene>
  mkdir -p "$outdir"

  # 为该基因裁剪子树（只用 nw_prune；失败即跳过）
  names_from_fa "$align" > "$outdir/taxa.keep"
  local ntaxa; ntaxa=$(wc -l < "$outdir/taxa.keep" 2>/dev/null || echo 0)
  if (( ntaxa < 3 )); then
    echo "[skip] $tag : ntaxa=$ntaxa (<3)" | tee -a "$LOGERR"
    for m in "${MODELS[@]}"; do echo -e "$tag\t$m\t\t\t" >> "$SUMMARY"; done
    return 0
  fi
  if ! nw_prune -v "$TREEFILE" $(tr '\n' ' ' < "$outdir/taxa.keep") > "$outdir/tree.nwk" 2>/dev/null; then
    echo "[WARN] prune failed for $tag" | tee -a "$LOGERR"
    for m in "${MODELS[@]}"; do echo -e "$tag\t$m\t\t\t" >> "$SUMMARY"; done
    return 0
  fi

  # 模型循环
  for model in "${MODELS[@]}"; do
    local model_dir="$outdir/$model"
    mkdir -p "$model_dir"

    # 最简 ctl（兼容 conda 版 codeml）
    write_ctl "$align" "$outdir/tree.nwk" "$model_dir/result" "${NSMAP[$model]}"

    # 跑 codeml（清理旧产物；错误不终止全局）
    (
      cd "$model_dir"
      rm -f result.txt mlc lnf rst rst1 2NG.* 4fold.* run.log
      codeml result.ctl > run.log 2>&1
    ) || true

    if [[ -s "$model_dir/result.txt" ]]; then
      # 屏幕提示 lnL
      awk '/lnL/ && /time/ {for(i=1;i<=NF;i++) if($i=="lnL"){printf("[ok] %-12s %-3s : lnL=%s\n","'"$tag"'", "'"$model"'", $(i+2)); exit}}' \
        "$model_dir/result.txt"
      # 汇总写入
      read -r lnL kappa omega < <(extract_stats "$model_dir/result.txt")
      echo -e "$tag\t$model\t$lnL\t$kappa\t$omega" >> "$SUMMARY"
    else
      echo "[WARN] no output: $tag $model" | tee -a "$LOGERR"
      tail -n 8 "$model_dir/run.log" 2>/dev/null || true
      echo -e "$tag\t$model\t\t\t" >> "$SUMMARY"
    fi
  done
}

# ================= gene 级别 =================
echo "[plan] Searching nuclear gene alignments..."
mapfile -t NU_GENES < <(find "$NU_GENE_DIR" -type f -name "*.codon.fas" | sort)
echo "[plan] Found ${#NU_GENES[@]} gene alignments."

for fas in "${NU_GENES[@]}"; do
  gene=$(basename "$fas" .codon.fas)
  run_codeml "$fas" "$gene" "$OUTROOT/gene/$gene"
done

echo "[done] All nuclear codeml runs finished."
echo "Summary TSV: $SUMMARY"
echo "Warn log   : $LOGERR"






#mt
#!/usr/bin/env bash
# === mt dN/dS (gene-level) — drop-in replacement ===

TREEFILE="/mnt/spareHD_2/oxphos_gene_tree/unrooted_NJ_tree_pruned.standard.tree"
MT_GENE_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13"

OUTROOT="/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/mt"
mkdir -p "$OUTROOT/gene"

SUMMARY="$OUTROOT/codeml_sites_summary.tsv"
: > "$SUMMARY"
echo -e "gene\tmodel\tlnL\tkappa\tomega" >> "$SUMMARY"

# 固定迭代顺序（关联数组遍历无序，避免结果顺序漂移）
declare -A NSMAP=([M0]=0 [M1a]=1 [M2a]=2 [M7]=7 [M8]=8)
MODELS=(M0 M1a M2a M7 M8)

grab_stats(){ # $1=result.txt -> echo "lnL kappa omega"
  local f="$1" lnL kappa omega
  lnL=$(awk   '/lnL/ && /time/ {for(i=1;i<=NF;i++) if($i=="lnL"){print $(i+2); exit}}' "$f" 2>/dev/null)
  kappa=$(awk '/kappa \(ts\/tv\)/ {print $4; exit}' "$f" 2>/dev/null)
  omega=$(awk '/omega \(dN\/dS\)/ {print $4; exit}' "$f" 2>/dev/null)
  echo "$lnL $kappa $omega"
}

run_codeml() {
    local align=$1
    local tag=$2
    local outdir=$3
    mkdir -p "$outdir"

    # 至少 3 条序列才有意义
    local ntaxa
    ntaxa=$(grep -c '^>' "$align" 2>/dev/null || echo 0)
    if (( ntaxa < 3 )); then
      echo "[skip] $tag : ntaxa=$ntaxa (<3)"
      for m in "${MODELS[@]}"; do echo -e "$tag\t$m\t\t\t" >> "$SUMMARY"; done
      return 0
    fi

    # 按对齐样本裁剪树（若无 nw_prune 则用原树）
    grep '^>' "$align" | sed 's/^>//; s/ .*//' > "$outdir/keep.txt"
    if command -v nw_prune >/dev/null 2>&1; then
      if ! nw_prune -v "$TREEFILE" $(tr '\n' ' ' < "$outdir/keep.txt") > "$outdir/tree.nwk" 2>/dev/null; then
        echo "[WARN] prune failed for $tag — using unpruned tree"
        cp -f "$TREEFILE" "$outdir/tree.nwk"
      fi
    else
      cp -f "$TREEFILE" "$outdir/tree.nwk"
    fi

    for model in "${MODELS[@]}"; do
        local model_dir="$outdir/$model"
        mkdir -p "$model_dir"

        cat > "$model_dir/codeml.ctl" <<EOF
seqfile = $align
treefile = $outdir/tree.nwk
outfile = $model_dir/result.txt

noisy = 3
verbose = 1
runmode = 0
seqtype = 1       * codon sequences
CodonFreq = 2
clock = 0
aaDist = 0
model = 0         * same omega distribution for all branches
NSsites = ${NSMAP[$model]}
icode = 1        * Vertebrate mitochondrial code (关键)
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
EOF
        ( cd "$model_dir" && rm -f result.txt mlc lnf rst rst1 2NG.* 4fold.* run.log
          codeml codeml.ctl > run.log 2>&1 || true )

        if [[ -s "$model_dir/result.txt" ]]; then
          read -r lnL kappa omega < <(grab_stats "$model_dir/result.txt")
          echo -e "$tag\t$model\t$lnL\t$kappa\t$omega" >> "$SUMMARY"
          printf "[ok] %-12s %-3s lnL=%s\n" "$tag" "$model" "$lnL"
        else
          echo -e "$tag\t$model\t\t\t" >> "$SUMMARY"
          echo "[WARN] no result: $tag $model"
        fi
    done
}

echo "[plan] Searching mtDNA gene alignments..."
mapfile -t MT_GENES < <(find "$MT_GENE_DIR" -type f -name "*.codon.fas" | sort)
echo "[plan] Found ${#MT_GENES[@]} gene alignments."

for fas in "${MT_GENES[@]}"; do
    gene=$(basename "$fas" .codon.fas)
    run_codeml "$fas" "$gene" "$OUTROOT/gene/$gene"
done

echo "Summary → $SUMMARY"





#merge into table
#!/usr/bin/env bash
set -euo pipefail

# ===== 路径：按需改 =====
ROOT="/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/nu/gene"
OUT="$ROOT/../codeml_sites_summary.tsv"   # 输出放在 nu/ 下面，避免和 gene/ 混一起

# 头
echo -e "gene\tmodel\tlnL\tkappa\tomega\tdN\tdS" > "$OUT"

# 遍历所有 result.txt
# 目录结构假定为: .../nu/gene/<gene>/<model>/result.txt
find "$ROOT" -type f -name "result.txt" | sort | while read -r F; do
  gene=$(basename "$(dirname "$(dirname "$F")")")
  model=$(basename "$(dirname "$F")")

  # 逐字段从 result.txt 中抓数；抓不到就 NA
  # lnL
  lnl=$(awk '/lnL\(/ {print $5; exit}' "$F")
  [[ -z "${lnl:-}" ]] && lnl="NA"

  # kappa：兼容 "kappa (ts/tv) =" 或 "kappa ="
  kappa=$(awk '
    $0 ~ /kappa/ {
      for(i=1;i<=NF;i++){
        if($i=="kappa" || index($i,"kappa")>0){
          # 找等号后的第一个数
          for(j=i;j<=NF;j++){
            if($(j)=="=" && (j+1)<=NF){print $(j+1); exit}
            if(match($(j),"=([0-9.+-eE]+)",m)){print m[1]; exit}
          }
        }
      }
    }' "$F" | head -n1)
  [[ -z "${kappa:-}" ]] && kappa="NA"

  # omega：兼容 "omega (dN/dS) =" 或 "w (dN/dS) ="，取第一个出现的
  omega=$(awk '
    /(^|[ \t])(omega|w)[ \t]*\(/ || /(^|[ \t])(omega|w)[ \t]*=/ {
      for(i=1;i<=NF;i++){
        if($i=="omega" || $i=="w" || index($i,"omega")>0){
          for(j=i;j<=NF;j++){
            if($(j)=="=" && (j+1)<=NF){print $(j+1); exit}
            if(match($(j),"=([0-9.+-eE]+)",m)){print m[1]; exit}
          }
        }
      }
    }' "$F" | head -n1)
  [[ -z "${omega:-}" ]] && omega="NA"

  # dN / dS 的 tree length
  dN=$(awk '/tree length for dN:/ {print $(NF); exit}' "$F")
  dS=$(awk '/tree length for dS:/ {print $(NF); exit}' "$F")
  [[ -z "${dN:-}" ]] && dN="NA"
  [[ -z "${dS:-}" ]] && dS="NA"

  echo -e "${gene}\t${model}\t${lnl}\t${kappa}\t${omega}\t${dN}\t${dS}" >> "$OUT"
done

# 简短反馈
echo "[done] wrote: $OUT"
echo "[hint] 预览前几行："
head -n 10 "$OUT"





#nu branch model1
#!/usr/bin/env bash
set -euo pipefail

# ====== 输入路径（保持与你现有布局一致）======
TREEFILE="/mnt/spareHD_2/oxphos_gene_tree/unrooted_NJ_tree_pruned.standard.tree"
NU_GENE_DIR="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"

# ====== 全新输出路径（不会覆盖你之前的 site models）======
OUTROOT="/mnt/spareHD_2/oxphos_codeml_ready/10_codeml_branch_models/nu_free_ratio"
mkdir -p "$OUTROOT/gene"
SUMMARY_BRANCHES="$OUTROOT/nu_free_ratio_tip_branches.all_genes.tsv"   # 全基因汇总
SUMMARY_LNL="$OUTROOT/nu_free_ratio_lnL.tsv"                           # 每基因 lnL 等
: > "$SUMMARY_BRANCHES"
: > "$SUMMARY_LNL"

# 依赖检查
for x in codeml nw_prune; do
  command -v "$x" >/dev/null 2>&1 || { echo "[FATAL] $x not found in PATH"; exit 1; }
done

# 写 ctl：free-ratio（branch model）
write_ctl(){  # $1 seqfile  $2 treefile  $3 outprefix
  cat > "$3.ctl" <<EOF
seqfile = $1
treefile = $2
outfile = $3.txt

noisy = 3
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 1         * free-ratio (每条分支独立 omega)
NSsites = 0
icode = 0         * 核基因标准遗传密码
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
EOF
}

# 从对齐里取样本顺序（PAML 将按此顺序编号 1..N）
emit_tip_order(){
  local fas="$1"
  grep '^>' "$fas" | sed 's/^>//; s/[ \t].*$//' | nl -ba
}

# 解析 free-ratio 的“dN & dS for each branch”表，抽出 tip 分支并加上样本名
# 输出列：gene  tip_idx  tip_name  t  N  S  omega  dN  dS  NdN  SdS
parse_fr_tips(){
  local gene="$1" res="$2" order_map="$3"
  # N: # of sequences
  local N; N=$(awk '/# of sequences/ {print $4; exit}' "$res" 2>/dev/null)
  [[ -z "${N:-}" ]] && return 0

  # 读顺序映射到数组 name[i]
  declare -A name
  while read -r idx nm; do name["$idx"]="$nm"; done < <(cat "$order_map")

  awk -v G="$gene" -v N="$N" -v MAP="$order_map" -v OFS="\t" '
    BEGIN{insec=0}
    $1=="dN" && $2=="&" && $3=="dS" && $4=="for" && $5=="each" && $6=="branch"{insec=1; getline; next}
    insec && NF>=9 {
      # 行格式：branch t N S dN/dS dN dS N*dN S*dS
      # 这里 $1 是形如 "28..1" 的分支名
      br=$1; gsub(/^ +| +$/,"",br)
      split(br,a,"..")
      left=a[1]+0; right=a[2]+0
      # tip 分支的判断：任一端点 <= N 且另一端 > N
      is_tip = ((left<=N && right>N) || (right<=N && left>N))
      if(!is_tip) next

      # 取 tip 的编号
      tip_idx = (left<=N ? left : right)

      # 填 tip 名：从 MAP 文件查
      tip_name="NA"
      # 为了速度在 awk 内读一次 MAP
      if(FNR==1){
        while( (getline L < MAP) > 0 ){
          split(L,x,"[ \t]+"); idx=x[1]; nm=x[2]
          if(nm!=""){ m[idx]=nm }
        }
        close(MAP)
      }
      if(tip_idx in m) tip_name=m[tip_idx]

      printf("%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
             G, tip_idx, tip_name, $(2), $(3), $(4), $(5), $(6), $(7), $(8), $(9))
    }' "$res"
}

# 抓 lnL/kappa/omega 行做简单记录
grab_lnL(){
  local gene="$1" res="$2"
  awk -v G="$gene" '
    /lnL\(ntime:/ {
      for(i=1;i<=NF;i++) if($i=="lnL"){lnl=$(i+2)}
    }
    /kappa \(ts\/tv\)/ {kappa=$4}
    /omega \(dN\/dS\)/ {omega=$4}
    END{
      printf("%s\t%s\t%s\t%s\n", G, (lnl?lnl:""), (kappa?kappa:""), (omega?omega:""))
    }' "$res"
}

echo -e "gene\ttip_idx\ttip_name\tt\tN\tS\tomega\tdN\tdS\tN*dN\tS*dS" >> "$SUMMARY_BRANCHES"
echo -e "gene\tlnL\tkappa\tomega_treewide" >> "$SUMMARY_LNL"

# ==== 主循环 ====
mapfile -t GENES < <(find "$NU_GENE_DIR" -type f -name "*.codon.fas" | sort)
echo "[plan] Found ${#GENES[@]} nuclear gene alignments."

for fas in "${GENES[@]}"; do
  gene=$(basename "$fas" .codon.fas)
  outdir="$OUTROOT/gene/$gene/FR"
  mkdir -p "$outdir"

  # 样本顺序映射（1..N → 名字）
  emit_tip_order "$fas" > "$outdir/tip_order.tsv"   # 两列：idx name

  # 用样本集合裁剪树
  cut -f2 "$outdir/tip_order.tsv" | tr '\n' ' ' | xargs -I{} bash -lc \
    'nw_prune -v "'"$TREEFILE"'" {} > "'"$outdir/tree.nwk"'"' || {
      echo "[WARN] prune failed for $gene"; continue;
    }

  # 写 ctl & 运行 codeml
  write_ctl "$fas" "$outdir/tree.nwk" "$outdir/result"
  ( cd "$outdir" && rm -f result.txt mlc lnf rst rst1 2NG.* 4fold.* run.log
    codeml result.ctl > run.log 2>&1 || true )

  if [[ ! -s "$outdir/result.txt" ]]; then
    echo "[WARN] no result for $gene"; continue
  fi

  # 解析 tip 分支、追加到总表
  parse_fr_tips "$gene" "$outdir/result.txt" "$outdir/tip_order.tsv" >> "$SUMMARY_BRANCHES"
  grab_lnL "$gene" "$outdir/result.txt" >> "$SUMMARY_LNL"
  echo "[ok] $gene : free-ratio done"
done

echo "[done] Free-ratio (nu) finished."
echo "  Tip-branches (all genes): $SUMMARY_BRANCHES"
echo "  lnL/kappa/omega summary : $SUMMARY_LNL"
