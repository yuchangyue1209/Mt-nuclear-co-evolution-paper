#!/usr/bin/env bash
set -euo pipefail

# ========= 配置 =========
REF="/work/cyu/stickleback_nuclear_only.fa"                                  # 参考基因组（可选：做翻译QC）
GFF="/work/cyu/stickleback_v5.gff3"                                          # GFF（未压缩；zcat -f 也兼容 .gz）
TABLE="/work/cyu/oxphos_from_ref_no_biomart/06_igv/final_genelist.csv"       # 你的基因表（含 stickleback_name 列）
OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv/dnds_annotations"        # 输出目录
mkdir -p "$OUTDIR"

echo "==> 1) 解析表格提取 stickleback_name 到 genes.list"
awk -F',' '
  NR==1 {for(i=1;i<=NF;i++){h=$i; gsub(/^[ \t]+|[ \t]+$/,"",h); if(h=="stickleback_name") col=i} next}
  col>0 && $col!="" {x=$col; gsub(/^[ \t]+|[ \t]+$/,"",x); print tolower(x)}
' "$TABLE" | sort -u > "$OUTDIR/genes.list"

NGEN=$(wc -l < "$OUTDIR/genes.list")
echo "🧬 目标基因数: $NGEN"
if [ "$NGEN" -eq 0 ]; then
  echo "❌ 未在 CSV 表头中找到 'stickleback_name' 或该列为空。"
  exit 1
fi

echo "==> 2) 从 GFF 抽取目标基因的转录本及其 CDS（保留 phase）"
# 生成候选转录本的全部 CDS 行（GFF3 风格）
awk -v GLIST="$OUTDIR/genes.list" 'BEGIN{
    FS=OFS="\t";
    while((getline g < GLIST)>0){ t[tolower(g)]=1 }
    close(GLIST)
}
$3=="gene"{
    name=gid="";
    if (match($9,/(^|;)Name=([^;]+)/,a)) name=a[2];
    if (match($9,/(^|;)ID=([^;]+)/,b))   gid=b[2];
    if (tolower(name) in t){ keepGene[gid]=1; gid2name[gid]=name }
    next
}
$3=="mRNA"{
    tx=pid="";
    if (match($9,/(^|;)ID=([^;]+)/,a))    tx=a[2];
    if (match($9,/(^|;)Parent=([^;]+)/,b)) pid=b[2];
    if (keepGene[pid]) keepTx[tx]=pid;
    next
}
$3=="CDS"{
    if (match($9,/(^|;)Parent=([^;]+)/,a)){
        tx=a[2];
        if (tx in keepTx) print $0
    }
}' "$GFF" > "$OUTDIR/target.CDS.gff3"

echo "候选 CDS 条数: $(wc -l < "$OUTDIR/target.CDS.gff3")"

echo "==> 3) 统计每个转录本 CDS 长度，并选择每个基因最长 CDS 的主转录本（纯 awk 合并）"
# 3A) 每转录本 CDS 总长度
awk 'BEGIN{FS=OFS="\t"}
     $3=="CDS"{L=$5-$4; match($9,/(^|;)Parent=([^;]+)/,a); tx=a[2]; len[tx]+=L}
     END{for(t in len) print t, len[t] }' \
  "$OUTDIR/target.CDS.gff3" > "$OUTDIR/tx_len.tsv"

# 3B) 为目标基因建立 tx -> geneID -> geneSymbol 映射
awk -v GLIST="$OUTDIR/genes.list" 'BEGIN{
  FS=OFS="\t";
  while((getline g < GLIST)>0){ tgt[tolower(g)]=1 } close(GLIST)
}
$3=="gene"{
  name=gid="";
  if (match($9,/(^|;)Name=([^;]+)/,a)) name=a[2];
  if (match($9,/(^|;)ID=([^;]+)/,b))   gid=b[2];
  gname[gid]=name; if (tolower(name) in tgt) keep[gid]=1; next
}
$3=="mRNA"{
  tx=pid="";
  if (match($9,/(^|;)ID=([^;]+)/,a))    tx=a[2];
  if (match($9,/(^|;)Parent=([^;]+)/,b)) pid=b[2];
  if (keep[pid]) print pid, tx, gname[pid]
}' "$GFF" > "$OUTDIR/gid_tx_name.targets.tsv"

# 3C) 选每个基因最长 CDS 的主转录本（不用 join）
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {len[$1]=$2; next}             # tx_len.tsv → len[tx]
     {gid=$1; tx=$2; gname=$3; L=len[tx]+0; if (L>bestL[gid]){bestL[gid]=L; bestTx[gid]=tx; bestName[gid]=gname}}
     END{for(g in bestTx) print g, bestTx[g], bestName[g], bestL[g] }' \
  "$OUTDIR/tx_len.tsv" "$OUTDIR/gid_tx_name.targets.tsv" \
| sort -k3,3 > "$OUTDIR/chosen_tx_per_gene.tsv"

NSEL=$(wc -l < "$OUTDIR/chosen_tx_per_gene.tsv")
echo "🏆 选出主转录本数：$NSEL（应接近 $NGEN）"

# 缺失基因报告：输入的 symbol 中哪些没有被选出（可能因 GFF 无 CDS 或别名不一致）
awk '{print tolower($1)}' "$OUTDIR/genes.list"                                   > "$OUTDIR/_all_symbols.tmp"
cut -f3 "$OUTDIR/chosen_tx_per_gene.tsv" | awk '{print tolower($1)}' | sort -u   > "$OUTDIR/_chosen_symbols.tmp"
comm -23 <(sort -u "$OUTDIR/_all_symbols.tmp") "$OUTDIR/_chosen_symbols.tmp"     > "$OUTDIR/missing_symbols.txt" || true
echo "🧭 缺失 symbol 列表：$OUTDIR/missing_symbols.txt（若非空，说明这些 symbol 在 GFF 里未匹配到带 CDS 的 mRNA）"
rm -f "$OUTDIR/_all_symbols.tmp" "$OUTDIR/_chosen_symbols.tmp"

echo "==> 4) 导出主转录本的 CDS 注释：GTF（含 phase）、BED6（0-based）、phase 明细 TSV"
cut -f2 "$OUTDIR/chosen_tx_per_gene.tsv" | sort -u > "$OUTDIR/chosen.tx.list"

# 4A) 仅保留主转录本的 CDS（保留原 phase），补 gene_name 属性
awk 'BEGIN{FS=OFS="\t"} FNR==NR{want[$1]=1; next}
     $3=="CDS"{ if (match($9,/(^|;)Parent=([^;]+)/,a)) {tx=a[2]; if (want[tx]) print $0} }' \
     "$OUTDIR/chosen.tx.list" "$OUTDIR/target.CDS.gff3" \
| awk -v MAP="$OUTDIR/chosen_tx_per_gene.tsv" 'BEGIN{
       FS=OFS="\t"; while((getline m < MAP)>0){ split(m,a,"\t"); tx2g[a[2]]=a[3] } close(MAP)
     }
     { match($9,/(^|;)Parent=([^;]+)/,a); tx=a[2]; g=(tx2g[tx]?tx2g[tx]:"NA");
       print $1,"Ensembl","CDS",$4,$5,$6,$7,$8,$9";gene_name="g }' \
> "$OUTDIR/oxphos_assembly.CDS.gtf"

# 4B) BED6 与 phase.tsv
: > "$OUTDIR/oxphos_assembly.CDS.bed"
: > "$OUTDIR/oxphos_assembly.CDS.phase.tsv"
awk -v MAP="$OUTDIR/chosen_tx_per_gene.tsv" '
BEGIN{
  FS=OFS="\t"; while((getline m < MAP)>0){ split(m,a,"\t"); tx2g[a[2]]=a[3] } close(MAP)
}
$3=="CDS"{
  if (match($9,/(^|;)Parent=([^;]+)/,p)){
    tx=p[2]; g=(tx2g[tx]?tx2g[tx]:"NA");
    bedStart=$4-1; bedEnd=$5; strand=$7; phase=$8; cnt[tx]++; name=g"|"tx"|CDS"cnt[tx];
    print $1,bedStart,bedEnd,name,0,strand >> "'"$OUTDIR/oxphos_assembly.CDS.bed"'";
    print $1,$4,$5,strand,phase,tx,g       >> "'"$OUTDIR/oxphos_assembly.CDS.phase.tsv"'";
  }
}' "$OUTDIR/oxphos_assembly.CDS.gtf"

echo "输出文件："
echo " - $OUTDIR/chosen_tx_per_gene.tsv           # 每基因选中的主转录本（geneID  txID  geneSymbol  cds_len）"
echo " - $OUTDIR/oxphos_assembly.CDS.gtf          # 给 dN/dS 用（GTF，含相位列）"
echo " - $OUTDIR/oxphos_assembly.CDS.bed          # 可视化 BED6（0-based）"
echo " - $OUTDIR/oxphos_assembly.CDS.phase.tsv    # 每段相位/转录本/基因明细"
echo " - $OUTDIR/missing_symbols.txt              # 未匹配到的 symbol（如非空需排查别名/注释）"

# ========= 5) 可选：翻译 QC（需要 gffread；若有 seqkit 更方便筛选） =========
if command -v gffread >/dev/null 2>&1; then
  echo "==> 5) 翻译 QC：导出主转录本的 CDS/蛋白序列（如需可进一步排查内部终止等）"
  # 直接全导出，再按 tx 过滤（兼容部分 gffread 不支持 -i 列表过滤的情况）
  gffread -g "$REF" -y "$OUTDIR/all.prot.fa" -x "$OUTDIR/all.cds.fa" -w "$OUTDIR/all.mrna.fa" "$GFF" || true
  if command -v seqkit >/dev/null 2>&1; then
    cut -f2 "$OUTDIR/chosen_tx_per_gene.tsv" > "$OUTDIR/chosen.tx.only"
    seqkit grep -f "$OUTDIR/chosen.tx.only" "$OUTDIR/all.prot.fa" > "$OUTDIR/chosen_tx.prot.fa" || true
    seqkit grep -f "$OUTDIR/chosen.tx.only" "$OUTDIR/all.cds.fa"  > "$OUTDIR/chosen_tx.cds.fa"  || true
    rm -f "$OUTDIR/all.prot.fa" "$OUTDIR/all.cds.fa" "$OUTDIR/all.mrna.fa" || true
  else
    echo "（提示）未检测到 seqkit，已导出 all.prot.fa / all.cds.fa，可手动筛选或安装 seqkit 后再筛。"
  fi
fi

echo "🎉 全流程完成。"


#delete