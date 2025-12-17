#erc for final oxphos genes
#bed for nuclear genes
#!/usr/bin/env bash
set -euo pipefail

REF=/work/cyu/stickleback_nuclear_only.fa
BED_DIR=/work/cyu/poolseq/PPalign_output/ann_nuclear_pergene_bed/named_symlinks_287
OUT_BED=/work/cyu/nuOXPHOS_genes_287.merged.bed

# 需要 .fai 来按参考顺序排序
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# 合并 287 个 BED → 一个去重、合并重叠、参考顺序的 BED
cat "$BED_DIR"/*.bed \
 | awk 'BEGIN{OFS="\t"} $1!~/^#/ && $2<$3' \
 | bedtools sort -faidx "${REF}.fai" -i - \
 | bedtools merge > "$OUT_BED"

echo "[BED] merged -> $OUT_BED  lines=$(wc -l < "$OUT_BED")"

#!/usr/bin/env bash
set -euo pipefail

BED=/work/cyu/nuOXPHOS_genes_287.merged.bed   # ← 若用你的 v2 总表，改成它
INPUT_DIR=/mnt/spareHD_2/nuclear_with_readgroups
OUTPUT_DIR=/mnt/spareHD_2/nuclear_with_readgroups_subset

mkdir -p "$OUTPUT_DIR"

# 并行度（按机器调）
THREADS=8

for BAM in "$INPUT_DIR"/*_rg.bam; do
  S=$(basename "$BAM" _rg.bam)
  OUT="$OUTPUT_DIR/${S}_subset.bam"

  echo "✂️  Subsetting $S ..."
  # 若无原索引，先建索引
  [ -f "${BAM}.bai" ] || samtools index -@ "$THREADS" "$BAM"

  # 截取 + 输出 BAM（保持坐标排序）
  samtools view -@ "$THREADS" -bh -L "$BED" "$BAM" -o "$OUT"

  # 建新索引
  samtools index -@ "$THREADS" "$OUT"
  echo "✅ Finished: $OUT"
done

echo "🎉 All BAMs subset to 287 nuclear OXPHOS regions!"


#!/usr/bin/env bash
set -euo pipefail

MT_BAMLIST=/work/cyu/poolseq/PPalign_output/mtDNA_bam/bamlist.txt
SUB_DIR=/mnt/spareHD_2/nuclear_with_readgroups_subset
NEW_BAMLIST=/mnt/spareHD_2/oxphos_gene_tree/bamlist_nuclear.txt

awk -v d="$SUB_DIR" '
  { if (match($0,/([0-9]+_[A-Z]+)/,a)) print d"/"a[1]"_subset.bam" }
' "$MT_BAMLIST" > "$NEW_BAMLIST"

echo "[lines] $(wc -l < "$NEW_BAMLIST")  ->  $NEW_BAMLIST"

# 快检：是否有缺失/空文件
xargs -a "$NEW_BAMLIST" -I{} bash -c '[ -s "{}" ] || echo "MISSING: {}"'



#mileup sync
#!/usr/bin/env bash
set -euo pipefail

REF=/work/cyu/stickleback_nuclear_only.fa
BED_DIR=/work/cyu/poolseq/PPalign_output/ann_nuclear_pergene_bed/named_symlinks_287
BAMLIST=/mnt/spareHD_2/oxphos_gene_tree/bamlist_nuclear.txt   # 子集后的 *_subset.bam 清单
MP_DIR=/mnt/spareHD_2/nu_287/mpileup
SYNC_DIR=/mnt/spareHD_2/nu_287/sync
CF_DIR=/mnt/spareHD_2/nu_287/cf_top2
mkdir -p "$MP_DIR" "$SYNC_DIR" "$CF_DIR"
[ -f "${REF}.fai" ] || samtools faidx "$REF"

for BED in "$BED_DIR"/*.bed; do
  GENE=$(basename "$BED" .bed)
  echo "🧬 mpileup $GENE"
  samtools mpileup \
    -B -f "$REF" -b "$BAMLIST" \
    -q 30 -Q 30 -d 5000 \
    -l "$BED" \
    > "$MP_DIR/${GENE}.mpileup"
done

# 快速 QC：看有多少空文件
echo "[empty mpileup files]: $(grep -L . "$MP_DIR"/*.mpileup | wc -l)"


2) mpileup → sync（PoPoolation2）
for MP in "$MP_DIR"/*.mpileup; do
  GENE=$(basename "$MP" .mpileup)
  echo "🔄 mpileup2sync $GENE"
  perl ~/popoolation2_1201/mpileup2sync.pl \
       --input "$MP" \
       --output "$SYNC_DIR/${GENE}.sync" \
       --fastq-type sanger \
       --min-qual 30
done








#mt
#!/usr/bin/env bash
set -euo pipefail

# ===== 必改路径 =====
GFF="/home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff"   # mt GFF
SYNC_ALL="/work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.sync"    # 用来取染色体ID
OUT_DIR="/work/cyu/poolseq/PPalign_output/ann_mt_pergene_bed"      # 输出目录
# 选择bed命名用哪个ID:  gene | locus | auto(优先gene, 其次locus_tag)
NAME_KEY="auto"
# ====================

mkdir -p "$OUT_DIR"/{all,pcg,trna,rrna}
CHR=$(awk 'NR==1{print $1; exit}' "$SYNC_ALL")
echo "[chr] use chromosome id: $CHR"

# 1) 从 GFF 抽 gene 级区间（0-based半开；4列：chr start end name）并分类
#    - name 按 NAME_KEY 选择：gene / locus_tag / auto
#    - type 分类: trna / rrna / pcg
awk -F'\t' -v OFS='\t' -v chr="$CHR" -v key="$NAME_KEY" '
  BEGIN{ IGNORECASE=1 }
  $0 ~ /^#/ { next }
  $3=="gene" {
    start0=$4-1; end=$5; attr=$9

    gene=""; locus=""; name=""
    if (match(attr, /gene=([^;]+)/, m))   gene=m[1]
    if (match(attr, /locus_tag=([^;]+)/, n)) locus=n[1]

    if (key=="gene" && gene!="")         name=gene
    else if (key=="locus" && locus!="")  name=locus
    else if (gene!="")                   name=gene
    else if (locus!="")                  name=locus
    else if (match(attr,/ID=([^;]+)/,p)) name=p[1]
    else name = sprintf("gene_%d_%d", $4, $5)

    # 分类：tRNA / rRNA（rns/rnl/12S/16S）/ PCG
    typ="pcg"
    if (gene ~ /^trn/i) typ="trna"
    else if (gene ~ /^(rns|rnl)$/i || gene ~ /(12s|16s)/i) typ="rrna"

    # 输出：chr start end name type
    print chr, start0, end, name, typ
  }
' "$GFF" > "$OUT_DIR/_all_genes.tsv"

echo "[count] genes total: $(wc -l < "$OUT_DIR/_all_genes.tsv")"

# 2) 写“总表BED”与按基因拆分的小BED
#    all_genes.bed（4列），并为每条记录写一个 all/<name>.bed
cut -f1-4 "$OUT_DIR/_all_genes.tsv" > "$OUT_DIR/all_genes.bed"

# 清理不安全文件名
sanitize(){ sed 's/[[:space:]]\+/_/g; s/[()]/_/g; s/[^A-Za-z0-9_.-]/_/g'; }

while IFS=$'\t' read -r c s e name typ; do
  safe=$(printf "%s\n" "$name" | sanitize)
  echo -e "$c\t$s\t$e\t$safe" > "$OUT_DIR/all/${safe}.bed"
  case "$typ" in
    trna) cp -f "$OUT_DIR/all/${safe}.bed" "$OUT_DIR/trna/${safe}.bed" ;;
    rrna) cp -f "$OUT_DIR/all/${safe}.bed" "$OUT_DIR/rrna/${safe}.bed" ;;
    *)    cp -f "$OUT_DIR/all/${safe}.bed" "$OUT_DIR/pcg/${safe}.bed"  ;;
  esac
done < "$OUT_DIR/_all_genes.tsv"

echo "[done] per-gene BEDs:"
echo "  - all  : $OUT_DIR/all/   (#$(ls -1 "$OUT_DIR/all"  | wc -l))"
echo "  - pcg  : $OUT_DIR/pcg/   (#$(ls -1 "$OUT_DIR/pcg"  | wc -l))"
echo "  - trna : $OUT_DIR/trna/  (#$(ls -1 "$OUT_DIR/trna" | wc -l))"
echo "  - rrna : $OUT_DIR/rrna/  (#$(ls -1 "$OUT_DIR/rrna" | wc -l))"


#mt
#!/usr/bin/env bash
set -euo pipefail

# ===== 路径（按需改） =====
SYNC_ALL="/work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.sync"               # 原始 fish.sync
BED_DIR="/work/cyu/poolseq/PPalign_output/ann_mt_pergene_bed/all"             # 你刚生成的单基因 BED 们
OUT_DIR="/work/cyu/poolseq/PPalign_output/ann_mt_pergene_sync"                # 新的 per-gene .sync 输出目录
# ========================

mkdir -p "$OUT_DIR"
echo "[0] 输出目录: $OUT_DIR"

# 0) 取染色体ID，并把 fish.sync 转成 bedtools 友好的“BED+附加列”
CHR=$(awk 'NR==1{print $1; exit}' "$SYNC_ALL")
SYNC_BED="$OUT_DIR/fish.sync.bed"
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $0}' "$SYNC_ALL" > "$SYNC_BED"
echo "[1] 写出: $SYNC_BED  (chr=$CHR)"

# 小工具：切一条 BED -> .sync（去掉前3个BED列，保留原始 sync 列）
cut_one() {  # $1 bedfile  $2 outpath
  bedtools intersect -wa -a "$SYNC_BED" -b "$1" \
  | sort -k1,1 -k2,2n -k3,3n -u \
  | awk 'BEGIN{OFS="\t"}{for(i=4;i<=NF;i++){printf i==NF? "%s\n":"%s\t",$i}}' \
  > "$2"
}

# 1) 批量切 per-gene .sync
shopt -s nullglob
n=0
for BED in "$BED_DIR"/*.bed; do
  gene=$(basename "$BED" .bed)
  out="$OUT_DIR/${gene}.sync"
  cut_one "$BED" "$out"
  ((++n))
  echo "   → ${gene}.sync"
done
shopt -u nullglob
echo "[2] 共生成 $n 个 .sync"

# 2) 快速 QC：检查每个 .sync 的首行是否为 CHR POS REF + 6计数模式
echo "[3] QC 首行格式（应为 [OK]）："
for S in "$OUT_DIR"/*.sync; do
  if awk 'NR==1{ ok=($3 ~ /^[ACGT]$/ && $4 ~ /^[0-9]+(:[0-9]+){5}$/); exit !ok }' "$S"
  then
    echo "[OK]   $(basename "$S")"
  else
    echo "[BAD]  $(basename "$S")"
    head -n1 "$S"
  fi
done








#cf mt top2

#!/usr/bin/env bash
set -euo pipefail

# ====== 路径与参数 ======
OUT_DIR="/work/cyu/poolseq/PPalign_output/ann_mt_pergene_sync"   # 你的 per-gene sync 新目录
SYNC_DIR="$OUT_DIR"                                              # 就用这个目录做输入
CF_DIR="/mnt/spareHD_2/mt_gene_tree/counts_top2"                 # 输出 cf
SCRIPT="/mnt/spareHD_2/mt_gene_tree/counts/sync2cf_pomo_fix.py"  # 你的 sync→cf 脚本
TOPN=2                                                           # PoMo 用 top2（推荐）

mkdir -p "$CF_DIR"

# POPS 顺序必须与 bamlist 一致（用 mtDNA 的 bamlist）
BAMLIST="/work/cyu/poolseq/PPalign_output/mtDNA_bam/bamlist.txt"
readarray -t POPS < <(awk -F'/' '{sub(/_mtDNA\.bam$/,"",$NF); print $NF}' "$BAMLIST")

echo "检测到 ${#POPS[@]} 个种群："
echo "${POPS[@]}"
echo

# ====== 遍历 per-gene sync ======
shopt -s nullglob
# 既兼容 *_fixed.sync 也兼容 *.sync
for S in "$SYNC_DIR"/*_fixed.sync "$SYNC_DIR"/*.sync; do
  # 如果同名的 *_fixed.sync 与 .sync 同时存在，优先 *_fixed.sync；下面这两行可避免把同一路径处理两次
  [[ "$S" == *"_fixed.sync" ]] || { [[ -f "${S%.sync}_fixed.sync" ]] && continue; }

  [[ -s "$S" ]] || { echo "⚠ 跳过空文件: $S"; continue; }

  BN=$(basename "$S")
  GENE=${BN%%_*}              # ATP6_fixed.sync → ATP6
  GENE=${GENE%.sync}          # 兜底

  O="$CF_DIR/${GENE}.cf"

  echo "🔄 生成 CF：$GENE"
  python "$SCRIPT" "$S" "$O" "${POPS[@]}" --topN "$TOPN"

  if [[ -s "$O" ]]; then
    NSITES=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="NSITES"){print $(i+1); exit}}' "$O")
    echo "✅ 写入 $O   (NSITES=${NSITES:-0})"
  else
    echo "❗ 失败：$O 为空"
  fi
done
shopt -u nullglob

echo -e "\n✅ counts(.cf) 全部输出到：$CF_DIR"


#!/usr/bin/env bash
set -euo pipefail

MASTER_TRE="/mnt/spareHD_2/oxphos_gene_tree/species_astral.tre"   # 核主树（AstraL）
CF_DIR="/mnt/spareHD_2/mt_gene_tree/counts_top2"                  # mt cf（top2）
OUT_DIR="/mnt/spareHD_2/mt_gene_tree/pruned_trees_top2"
THREADS=8

mkdir -p "$OUT_DIR/logs"

# 1) 提取 master 的物种名
tr '(),:;' '\n' < "$MASTER_TRE" | grep -vE '^$|^[0-9.]$|^[0-9]+\.[0-9]+$' \
| sort -u > /tmp/master.tips

for CF in "$CF_DIR"/*.cf; do
  GENE=$(basename "$CF" .cf)
  LOG="$OUT_DIR/logs/${GENE}.log"

  # 已有 treefile 则跳过（想重跑改为：if false && [[ -s ... ]] ）
  if [[ -s "$OUT_DIR/${GENE}_fixed.treefile" ]]; then
    echo "⏩ 跳过 $GENE（已存在输出）"
    continue
  fi

  echo "🌿 处理 $GENE ..."
  # 2) 从 CF 抓物种列（去 \r；匹配“CHROM  POS”）
  tr -d '\r' < "$CF" \
  | awk 'BEGIN{FS="[ \t]+"} $1 ~ /^#?CHROM$/ && $2=="POS"{for(i=3;i<=NF;i++) print $i; exit}' \
  | sort -u > /tmp/cf.taxa

  # 3) 与 master 取交集并检查数量
  comm -12 /tmp/master.tips /tmp/cf.taxa > /tmp/keep.taxa
  NUM_TAXA=$(wc -l < /tmp/keep.taxa)
  NSITES=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="NSITES"){print $(i+1); exit}}' "$CF")

  if (( NUM_TAXA < 4 )) || [[ -z "$NSITES" || "$NSITES" -lt 1 ]]; then
    echo "⚠️  跳过 $GENE ：taxa=$NUM_TAXA, NSITES=${NSITES:-0}"
    continue
  fi

  # 4) 固定 topology 重新估 branch length（PoMo）
  {
    echo "[info] $GENE taxa=$NUM_TAXA NSITES=$NSITES"
    iqtree2 -s "$CF" \
            -m GTR+P \
            -g "$MASTER_TRE" \
            -nt "$THREADS" \
            --safe \
            -pre "$OUT_DIR/${GENE}_fixed" \
            -quiet
  } >"$LOG" 2>&1 || { echo "❌ $GENE 失败，查看 $LOG"; continue; }

  if [[ -s "$OUT_DIR/${GENE}_fixed.treefile" ]]; then
    echo "✅ $GENE 完成（log: $LOG）"
  else
    echo "❗ $GENE 无 treefile（log: $LOG）"
  fi
done

echo -e "\n🎯 DONE → $OUT_DIR"



#nu cf top2
#!/usr/bin/env bash
set -euo pipefail

# —— 路径（按你的实际放置修改）——
BAMLIST_NU=/mnt/spareHD_2/oxphos_gene_tree/bamlist_nuclear.txt   # 你之前做好的核 bamlist（顺序=物种顺序）
SYNC_DIR=/mnt/spareHD_2/nu_287/sync                    # 核基因 *.sync 所在目录
CF_DIR=/mnt/spareHD_2/oxphos_gene_tree/counts_top2                # ← 新输出：核的 top2 cf
SCRIPT=/mnt/spareHD_2/mt_gene_tree/counts/sync2cf_pomo_fix.py     # 同一个脚本可通用
TOPN=2                                                            # ← 关键：top2

mkdir -p "$CF_DIR"

# 按 bamlist 提取物种列顺序（必须与 mpileup2sync 一致）
readarray -t POPS < <(awk -F'/' '{sub(/_subset\.bam$/,"",$NF); print $NF}' "$BAMLIST_NU")
echo "NPOP=${#POPS[@]}"; printf '%s ' "${POPS[@]}"; echo

shopt -s nullglob
for S in "$SYNC_DIR"/*.sync; do
  GENE=$(basename "$S" .sync)
  O="$CF_DIR/${GENE}.cf"

  [[ -s "$S" ]] || { echo "skip empty $S"; continue; }

  echo "🔄 PoMo cf  → $GENE (topN=$TOPN)"
  python "$SCRIPT" "$S" "$O" "${POPS[@]}" --topN "$TOPN"

  # 打印 NSITES
  awk 'NR==1{for(i=1;i<=NF;i++) if($i=="NSITES"){print FILENAME": NSITES="$(i+1)}}' "$O"
done
shopt -u nullglob

echo "✅ 核基因 PoMo cf(top2) 写入：$CF_DIR"






# mt prune tree &iqtree

MASTER_TRE="/mnt/spareHD_2/oxphos_gene_tree/species_astral.tre"
CF_DIR="/mnt/spareHD_2/mt_gene_tree/counts_top2"
OUT_DIR="/mnt/spareHD_2/mt_gene_tree/pruned_trees_te"
THREADS=8

mkdir -p "$OUT_DIR"

# 提取 master tree 物种名
tr '(),:;' '\n' < "$MASTER_TRE" | grep -vE '^$|^[0-9.]+$' | sort -u > /tmp/master.tips

for CF in "$CF_DIR"/*.cf; do
    GENE=$(basename "$CF" .cf)
    echo "🌿 处理 $GENE ..."

    # 从 .cf 文件头提取物种名
    awk 'NR==2 && $1=="CHROM"{for(i=3;i<=NF;i++) print $i}' "$CF" | sort -u > /tmp/cf.taxa

    # 取交集
    comm -12 /tmp/master.tips /tmp/cf.taxa > /tmp/keep.taxa
    NUM_TAXA=$(wc -l < /tmp/keep.taxa)

    if (( NUM_TAXA < 4 )); then
        echo "⚠️  $GENE 物种数太少 ($NUM_TAXA)，跳过"
        continue
    fi

    # 固定拓扑，PoMo 模型，重新估 branch length
    iqtree2 -s "$CF" \
            -m GTR+P \
            -te "$MASTER_TRE" \
            -nt "$THREADS" \
            -blmin 1e-12 -blmax 100 \
            -pre "$OUT_DIR/${GENE}_fixed" \
            --safe \
            -quiet \
    || { echo "❌ $GENE 运行失败，跳过"; continue; }

    echo "✅ $GENE 完成"
done



# nu prune tree & iqtree
MASTER_TRE="/mnt/spareHD_2/oxphos_gene_tree/species_astral.tre"
CF_DIR="/mnt/spareHD_2/oxphos_gene_tree/counts_top2"
OUT_DIR="/mnt/spareHD_2/oxphos_gene_tree/pruned_trees_te"
THREADS=8

mkdir -p "$OUT_DIR"

# 提取 master tree 物种名（去掉分支长度和支持值）
tr '(),:;' '\n' < "$MASTER_TRE" | grep -vE '^$|^[0-9.]+$' | sort -u > /tmp/master.tips

for CF in "$CF_DIR"/*.cf; do
    GENE=$(basename "$CF" .cf)
    echo "🌿 处理 $GENE ..."

    # 从 .cf 文件第2行提取物种名
    awk 'NR==2 && $1=="CHROM"{for(i=3;i<=NF;i++) print $i}' "$CF" | sort -u > /tmp/cf.taxa

    # 取交集
    comm -12 /tmp/master.tips /tmp/cf.taxa > /tmp/keep.taxa

    NUM_TAXA=$(wc -l < /tmp/keep.taxa)
    if (( NUM_TAXA < 4 )); then
        echo "⚠️  $GENE 物种数太少 ($NUM_TAXA)，跳过"
        continue
    fi

    # 用 IQ-TREE 固定 topology 重新估分支长度
    iqtree2 -s "$CF" \
            -m GTR+P \
            -te "$MASTER_TRE" \
            -nt "$THREADS" \
            --safe \
            -pre "$OUT_DIR/${GENE}_fixed" \
            -quiet \
    || { echo "❌ $GENE 运行失败，跳过"; continue; }

    echo "✅ $GENE 完成"
done

echo -e "\n🎯 所有基因树已剪枝并重新估分支长度，结果在 $OUT_DIR"





# merged tree
# 核：把 *_fixed.treefile 合并并在每棵树前加 [gene]
cd /mnt/spareHD_2/oxphos_gene_tree/pruned_trees
for f in $(ls *_fixed.treefile | sort); do
  gene=${f%%_*}
  printf '[%s]\n' "$gene"
  cat "$f"
  printf '\n'
done > /mnt/spareHD_2/oxphos_gene_tree/genes_astral_named.tre

# 线粒体：top2 目录
cd /mnt/spareHD_2/mt_gene_tree/pruned_trees
for f in $(ls *_fixed.treefile | sort); do
  gene=${f%%_*}          # ATP6/ND1/… 这些前缀
  printf '[%s]\n' "$gene"
  cat "$f"
  printf '\n'
done > /mnt/spareHD_2/mt_gene_tree/mt_genes_named.tre


cat \
  /mnt/spareHD_2/oxphos_gene_tree/genes_astral_named.tre \
  /mnt/spareHD_2/mt_gene_tree/mt_genes_named.tre \
  > /mnt/spareHD_2/all_genes_named.tre

# 再数一下总基因数（= 以 [gene] 开头的树数）
grep -c '^\[' /mnt/spareHD_2/all_genes_named.tre



# merged tree te
# 核：把 *_fixed.treefile 合并并在每棵树前加 [gene]
cd /mnt/spareHD_2/oxphos_gene_tree/pruned_trees_te
for f in $(ls *_fixed.treefile | sort); do
  gene=${f%%_*}
  printf '[%s]\n' "$gene"
  cat "$f"
  printf '\n'
done > /mnt/spareHD_2/oxphos_gene_tree/te_genes_astral_named.tre

# 线粒体：top2 目录
cd /mnt/spareHD_2/mt_gene_tree/pruned_trees_te
for f in $(ls *_fixed.treefile | sort); do
  gene=${f%%_*}          # ATP6/ND1/… 这些前缀
  printf '[%s]\n' "$gene"
  cat "$f"
  printf '\n'
done > /mnt/spareHD_2/mt_gene_tree/te_mt_genes_named.tre


cat \
  /mnt/spareHD_2/oxphos_gene_tree/te_genes_astral_named.tre \
  /mnt/spareHD_2/mt_gene_tree/te_mt_genes_named.tre \
  > /mnt/spareHD_2/te_all_genes_named.tre

# 再数一下总基因数（= 以 [gene] 开头的树数）
grep -c '^\[' /mnt/spareHD_2/all_genes_named.tre




MT_TE_DIR="/mnt/spareHD_2/mt_gene_tree/mt_gene_trees_te"
THREADS=8

mkdir -p "$MT_TE_DIR"

shopt -s nullglob
for ALN in "$ALIGN_DIR"/*.mt.aln.fasta; do
  GENE=$(basename "$ALN" .mt.aln.fasta)
  PRE="$MT_TE_DIR/${GENE}_te"

  if [[ -s "${PRE}.treefile" ]]; then
    echo "⏩ [mt] 跳过 $GENE（*_te.treefile 已存在）"
    continue
  fi

  echo "🌲 [mt] IQ-TREE 约束建树 $GENE ..."

  iqtree2 -s "$ALN" \
          -m GTR+G \        # 或者用 MFP 也行：-m MFP
          -te "$REF_TRE" \  # 拓扑固定为 species_astral
          -nt "$THREADS" \
          -keep-ident \
          -pre "$PRE" \
          --safe \
          -quiet

  if [[ -s "${PRE}.treefile" ]]; then
    echo "✅ [mt] $GENE → ${PRE}.treefile"
  else
    echo "❗ [mt] $GENE 失败，请查 ${PRE}.log"
  fi
done
shopt -u nullglob

echo -e "\n🎯 [mt] 约束 mt gene trees 在：$MT_TE_DIR"
