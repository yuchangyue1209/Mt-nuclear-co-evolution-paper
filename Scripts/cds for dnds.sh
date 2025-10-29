#!/usr/bin/env bash
set -euo pipefail

# ========= 配置 =========
REF="/work/cyu/stickleback_nuclear_only.fa"                                  # 可选：翻译QC用
GFF="/work/cyu/stickleback_v5.gff3"                                          # GFF（未压缩也可）
TABLE="/work/cyu/oxphos_from_ref_no_biomart/06_igv/final.csv"                # 你的最终表（含 stickleback_name）
OUTDIR="/work/cyu/oxphos_from_ref_no_biomart/06_igv/dnds_annotations_final"  # 输出目录（新的）
EXPECTED_N=291                                                                # 你预期的最终条数（可改）
mkdir -p "$OUTDIR"

# ====== 允许“保留两次”的基因（小写）======
# 把下面两个占位符换成你要保留两次的 gene symbol，小写，例如：("rpl5" "rps4x")
ALLOW_DUP=("symbolA" "symbolB")

# 是否排除 chrUn（建议 1）
EXCLUDE_CHRUN=1

echo "==> 0) 解析 final.csv: stickleback_name -> genes.list（小写去重）"
awk -F',' '
  NR==1 {
    for(i=1;i<=NF;i++){
      h=$i; gsub(/^[ \t]+|[ \t]+$/,"",h);
      if(h=="stickleback_name") col=i
    }
    next
  }
  col>0 && $col!="" {x=$col; gsub(/^[ \t]+|[ \t]+$/,"",x); print tolower(x)}
' "$TABLE" | sort -u > "$OUTDIR/genes.list"

NGEN=$(wc -l < "$OUTDIR/genes.list")
echo "🧬 目标基因（unique symbol）: $NGEN"

echo "==> 1) 从 GFF 抽取目标基因的转录本与其 CDS（保留 phase）"
awk -v GLIST="$OUTDIR/genes.list" 'BEGIN{
    FS=OFS="\t"; while((getline g<GLIST)>0){t[tolower(g)]=1} close(GLIST)
}
$3=="gene"{
  name=gid="";
  if(match($9,/(^|;)Name=([^;]+)/,a)) name=a[2];
  if(match($9,/(^|;)ID=([^;]+)/,b))   gid=b[2];
  if(t[tolower(name)]){ keepGene[gid]=1; gid2name[gid]=name }
  next
}
$3=="mRNA"{
  tx=pid="";
  if(match($9,/(^|;)ID=([^;]+)/,a))    tx=a[2];
  if(match($9,/(^|;)Parent=([^;]+)/,b)) pid=b[2];
  if(keepGene[pid]) keepTx[tx]=pid;
  next
}
$3=="CDS"{
  if(match($9,/(^|;)Parent=([^;]+)/,a)){
    tx=a[2];
    if(keepTx[tx]) print $0
  }
}' "$GFF" > "$OUTDIR/target.CDS.gff3"

echo "候选 CDS 条数: $(wc -l < "$OUTDIR/target.CDS.gff3")"

echo "==> 2) 统计每个转录本 CDS 长度；建立 tx→geneID→symbol 映射"
awk 'BEGIN{FS=OFS="\t"}
     $3=="CDS"{L=$5-$4; match($9,/(^|;)Parent=([^;]+)/,a); tx=a[2]; len[tx]+=L}
     END{for(t in len) print t,len[t] }' \
  "$OUTDIR/target.CDS.gff3" > "$OUTDIR/tx_len.tsv"

awk -v GLIST="$OUTDIR/genes.list" 'BEGIN{
  FS=OFS="\t"; while((getline g<GLIST)>0){t[tolower(g)]=1} close(GLIST)
}
$3=="gene"{
  name=gid="";
  if(match($9,/(^|;)Name=([^;]+)/,a)) name=a[2];
  if(match($9,/(^|;)ID=([^;]+)/,b))   gid=b[2];
  gname[gid]=name; if(t[tolower(name)]) keep[gid]=1; next
}
$3=="mRNA"{
  tx=pid="";
  if(match($9,/(^|;)ID=([^;]+)/,a))    tx=a[2];
  if(match($9,/(^|;)Parent=([^;]+)/,b)) pid=b[2];
  if(keep[pid]) print pid,tx,gname[pid]
}' "$GFF" > "$OUTDIR/gid_tx_name.targets.tsv"

echo "==> 3) 为每个 gene 选“最长 CDS”的主转录本（未去重）"
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR{len[$1]=$2; next}
     {gid=$1; tx=$2; sym=$3; L=len[tx]+0; if(L>bestL[gid]){bestL[gid]=L; bestTx[gid]=tx; bestName[gid]=sym}}
     END{for(g in bestTx) print g,bestTx[g],bestName[g],bestL[g] }' \
  "$OUTDIR/tx_len.tsv" "$OUTDIR/gid_tx_name.targets.tsv" \
| sort -k3,3 > "$OUTDIR/chosen_tx_per_gene.tsv"

NSEL=$(wc -l < "$OUTDIR/chosen_tx_per_gene.tsv")
echo "主转录本记录数（未按 symbol 去重）: $NSEL"

echo "==> 4) 构建 tx→(chr, chrUn标记)；做去重（优先非chrUn、再看CDS长）+ 白名单放行"
# tx → 是否 chrUn
awk 'BEGIN{FS=OFS="\t"}
     $3=="CDS" && match($9,/(^|;)Parent=([^;]+)/,a){
       tx=a[2]; chr=$1; un=(chr ~ /^chrUn/)?1:0;
       if(!(tx in seen)){print tx,un,chr; seen[tx]=1}
     }' "$OUTDIR/target.CDS.gff3" > "$OUTDIR/tx_chr_flag.tsv"   # tx  un_flag  chr

# 白名单
printf "%s\n" "${ALLOW_DUP[@]}" | awk '{print tolower($0)}' | sort -u > "$OUTDIR/dup_whitelist.symbols"

# 生成“最终映射”：白名单允许重复（可选地仍排除 chrUn），非白名单按策略唯一化
awk -v FLG="$OUTDIR/tx_chr_flag.tsv" -v WL="$OUTDIR/dup_whitelist.symbols" -v EXC="$EXCLUDE_CHRUN" '
BEGIN{
  FS=OFS="\t";
  while((getline f<FLG)>0){split(f,a,"\t"); tx_un[a[1]]=a[2]+0; tx_chr[a[1]]=a[3]} close(FLG);
  while((getline w<WL)>0){wl[w]=1} close(WL);
}
{
  gid=$1; tx=$2; sym=tolower($3); len=$4+0; un=(tx in tx_un?tx_un[tx]:1);
  if(wl[sym]){
    if(EXC && un==1) next;       # 白名单仍排除 chrUn（如需要）
    keep[++n]=$0
  } else {
    if(EXC && un==1) next;
    if(!(sym in best) || un<best_un[sym] || (un==best_un[sym] && len>best_len[sym])){
      best[sym]=$0; best_un[sym]=un; best_len[sym]=len
    }
  }
}
END{
  for(i=1;i<=n;i++) print keep[i];
  for(s in best)    print best[s];
}' "$OUTDIR/chosen_tx_per_gene.tsv" \
| sort -k3,3 > "$OUTDIR/chosen_tx_per_gene.final.tsv"

NFINAL=$(wc -l < "$OUTDIR/chosen_tx_per_gene.final.tsv")
echo "✅ 最终映射条数（含白名单放行、无 chrUn）: $NFINAL"

# 检查是否命中期望
if [ "$EXPECTED_N" -gt 0 ]; then
  if [ "$NFINAL" -ne "$EXPECTED_N" ]; then
    echo "⚠️ 最终计数($NFINAL) != 期望($EXPECTED_N)。请检查白名单/chrUn/重复策略。"
  fi
fi

echo "==> 5) 导出 final GTF（含 phase）、BED、phase.tsv"
cut -f2 "$OUTDIR/chosen_tx_per_gene.final.tsv" | sort -u > "$OUTDIR/chosen.tx.final.list"

# GTF（含 phase & gene_name）
awk 'BEGIN{FS=OFS="\t"} FNR==NR{want[$1]=1; next}
     $3=="CDS" && match($9,/(^|;)Parent=([^;]+)/,a){tx=a[2]; if(want[tx]) print}' \
     "$OUTDIR/chosen.tx.final.list" "$OUTDIR/target.CDS.gff3" \
| awk -v MAP="$OUTDIR/chosen_tx_per_gene.final.tsv" 'BEGIN{
       FS=OFS="\t"; while((getline m<MAP)>0){split(m,a,"\t"); tx2g[a[2]]=a[3]} close(MAP)
     }{
       match($9,/(^|;)Parent=([^;]+)/,a); tx=a[2]; g=(tx2g[tx]?tx2g[tx]:"NA");
       print $1,"Ensembl","CDS",$4,$5,$6,$7,$8,$9";gene_name="g
     }' > "$OUTDIR/oxphos_assembly.CDS.final.gtf"

# BED & phase.tsv
: > "$OUTDIR/oxphos_assembly.CDS.final.bed"
: > "$OUTDIR/oxphos_assembly.CDS.final.phase.tsv"
awk 'BEGIN{FS=OFS="\t"}
     $3=="CDS"{
       name="NA"; tx="";
       if(match($9,/gene_name=([^;]+)/,g)) name=g[1];
       if(match($9,/(^|;)Parent=([^;]+)/,p)) tx=p[2];
       cnt[tx]++; print $1,$4-1,$5,name"|"tx"|CDS"cnt[tx],0,$7 >> "'"$OUTDIR/oxphos_assembly.CDS.final.bed"'";
       print $1,$4,$5,$7,$8,tx,name >> "'"$OUTDIR/oxphos_assembly.CDS.final.phase.tsv"'";
     }' "$OUTDIR/oxphos_assembly.CDS.final.gtf"

echo "==> 6) 对账与质控报告"
# 是否含 chrUn
CHRUN_LINES=$(awk '$1 ~ /^chrUn/' "$OUTDIR/oxphos_assembly.CDS.final.gtf" | wc -l)
echo "chrUn 行数（应为 0）: $CHRUN_LINES"

# 列出最终 symbol 集合
cut -f3 "$OUTDIR/chosen_tx_per_gene.final.tsv" | awk '{print tolower($1)}' | sort -u > "$OUTDIR/final_symbols.list"
echo "最终 unique symbol 个数: $(wc -l < "$OUTDIR/final_symbols.list")"

# 与 CSV 对齐：CSV unique - 最终 = 被剔除者
awk -F',' '
  NR==1{for(i=1;i<=NF;i++){h=$i; gsub(/^[ \t]+|[ \t]+$/,"",h); if(h=="stickleback_name") col=i} next}
  col>0 && $col!="" {x=$col; gsub(/^[ \t]+|[ \t]+$/,"",x); print tolower(x)}
' "$TABLE" | sort -u > "$OUTDIR/csv_symbols.unique.txt"

comm -23 "$OUTDIR/csv_symbols.unique.txt" "$OUTDIR/final_symbols.list" > "$OUTDIR/symbols_dropped.txt" || true
echo "CSV unique 中被剔除的 symbol 数: $(wc -l < "$OUTDIR/symbols_dropped.txt")"

# 白名单基因各自保留的条数
printf "%s\n" "${ALLOW_DUP[@]}" | awk '{print tolower($0)}' \
| while read s; do
    c=$(cut -f3 "$OUTDIR/chosen_tx_per_gene.final.tsv" | awk -v s="$s" 'BEGIN{IGNORECASE=1} tolower($0)==s' | wc -l)
    echo "[whitelist] " "$s -> $c"
  done

echo "==> 7) （可选）翻译 QC（仅核参考会缺 chrM；不影响 dN/dS 注释使用）"
if command -v gffread >/dev/null 2>&1; then
  if [ -f "$REF" ]; then
    # 过滤到参考里存在的染色体
    if command -v samtools >/dev/null 2>&1; then
      samtools faidx "$REF" >/dev/null 2>&1 || true
      cut -f1 "${REF}.fai" | sort -u > "$OUTDIR/ref_contigs.list"
      awk 'NR==FNR{ok[$1]=1; next} ($1 in ok)' \
          "$OUTDIR/ref_contigs.list" \
          "$OUTDIR/oxphos_assembly.CDS.final.gtf" > "$OUTDIR/oxphos_assembly.CDS.final.inREF.gtf"
      gffread -g "$REF" -y "$OUTDIR/final.prot.fa" -x "$OUTDIR/final.cds.fa" "$OUTDIR/oxphos_assembly.CDS.final.inREF.gtf" || true
    else
      gffread -g "$REF" -y "$OUTDIR/final.prot.fa" -x "$OUTDIR/final.cds.fa" "$OUTDIR/oxphos_assembly.CDS.final.gtf" || true
    fi
  else
    echo "（提示）REF 不存在：$REF，跳过翻译 QC。"
  fi
fi

echo "🎉 完成。产物位于：$OUTDIR"
echo " - chosen_tx_per_gene.final.tsv"
echo " - oxphos_assembly.CDS.final.gtf / .bed / .phase.tsv"
echo " - final_symbols.list / symbols_dropped.txt"




OUT="/work/cyu/oxphos_from_ref_no_biomart/06_igv/dnds_annotations_final"
EXCLUDE_CHRUN=1   # 仍然排除 chrUn

# 1) 白名单（小写）
printf "%s\n" rpl5 rpl26 > "$OUT/dup_whitelist.symbols"

# 2) 重新做“去重+白名单放行”
awk 'BEGIN{FS=OFS="\t"}
  $3=="CDS" && match($9,/(^|;)Parent=([^;]+)/,a){
    tx=a[2]; chr=$1; un=(chr ~ /^chrUn/)?1:0;
    if(!(tx in seen)){print tx,un,chr; seen[tx]=1}
  }' "$OUT/target.CDS.gff3" > "$OUT/tx_chr_flag.tsv"

awk -v FLG="$OUT/tx_chr_flag.tsv" -v WL="$OUT/dup_whitelist.symbols" -v EXC="$EXCLUDE_CHRUN" '
BEGIN{
  FS=OFS="\t";
  while((getline f<FLG)>0){split(f,a,"\t"); tx_un[a[1]]=a[2]+0} close(FLG);
  while((getline w<WL)>0){wl[w]=1} close(WL);
}
{
  gid=$1; tx=$2; sym=tolower($3); len=$4+0; un=(tx in tx_un?tx_un[tx]:1);
  if(wl[sym]) { if(EXC && un==1) next; keep[++n]=$0 }
  else {
    if(EXC && un==1) next;
    if(!(sym in best) || un<best_un[sym] || (un==best_un[sym] && len>best_len[sym])){best[sym]=$0; best_un[sym]=un; best_len[sym]=len}
  }
}
END{for(i=1;i<=n;i++) print keep[i]; for(s in best) print best[s] }
' "$OUT/chosen_tx_per_gene.tsv" | sort -k3,3 > "$OUT/chosen_tx_per_gene.final.tsv"

# 3) 导出 final GTF / BED / phase
cut -f2 "$OUT/chosen_tx_per_gene.final.tsv" | sort -u > "$OUT/chosen.tx.final.list"
awk 'BEGIN{FS=OFS="\t"} FNR==NR{want[$1]=1; next}
     $3=="CDS" && match($9,/(^|;)Parent=([^;]+)/,a){tx=a[2]; if(want[tx]) print}' \
     "$OUT/chosen.tx.final.list" "$OUT/target.CDS.gff3" \
| awk -v MAP="$OUT/chosen_tx_per_gene.final.tsv" 'BEGIN{
       FS=OFS="\t"; while((getline m<MAP)>0){split(m,a,"\t"); tx2g[a[2]]=a[3]} close(MAP)
     }{match($9,/(^|;)Parent=([^;]+)/,a); tx=a[2]; g=tx2g[tx];
       print $1,"Ensembl","CDS",$4,$5,$6,$7,$8,$9";gene_name="g }' > "$OUT/oxphos_assembly.CDS.final.gtf"

: > "$OUT/oxphos_assembly.CDS.final.bed"
: > "$OUT/oxphos_assembly.CDS.final.phase.tsv"
awk 'BEGIN{FS=OFS="\t"}
     $3=="CDS"{
       if(match($9,/gene_name=([^;]+)/,g)) name=g[1];
       if(match($9,/(^|;)Parent=([^;]+)/,p)) tx=p[2];
       cnt[tx]++; print $1,$4-1,$5,name"|"tx"|CDS"cnt[tx],0,$7 >> "'"$OUT/oxphos_assembly.CDS.final.bed"'";
       print $1,$4,$5,$7,$8,tx,name >> "'"$OUT/oxphos_assembly.CDS.final.phase.tsv"'";
     }' "$OUT/oxphos_assembly.CDS.final.gtf"

# 4) 快速核对
echo "Final 映射条数（期望=291）：$(wc -l < "$OUT/chosen_tx_per_gene.final.tsv")"
echo "chrUn 行数（应=0）：$(awk '$1 ~ /^chrUn/' "$OUT/oxphos_assembly.CDS.final.gtf" | wc -l)"
for s in rpl5 rpl26; do
  c=$(cut -f3 "$OUT/chosen_tx_per_gene.final.tsv" | awk -v s="$s" 'BEGIN{IGNORECASE=1} tolower($0)==s' | wc -l)
  echo "$s -> $c"
done


#04
#!/usr/bin/env bash
set -euo pipefail

# ======= inputs =======
GFF3="/work/cyu/stickleback_v5.gff3"    # 若是 .gff3.gz 也行，脚本会自动 zcat
TXLIST="/work/cyu/oxphos_from_ref_no_biomart/06_igv/dnds_annotations_final/chosen.tx.final.list"

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

echo "[03] input GFF3: $GFF3"
echo "[03] tx list    : $TXLIST"
echo "[03] out dir    : $OUT"

# 0) 读入目标转录本集合
TX_N=$(grep -vc '^\s*$' "$TXLIST" || true)
echo "[03] transcripts in list: ${TX_N}"

# 1) 建 transcript->gene 映射（mRNA/transcript 的 ID 和 Parent）
T2G="${OUT}/t2g.tsv"
$READER "$GFF3" \
| awk -F'\t' '
  $0!~/^#/ && ($3=="mRNA"||$3=="transcript"){
    tid=""; gid="";
    if (match($9,/ID=([^;]+)/,I))   tid=I[1];
    if (match($9,/Parent=([^;]+)/,P)) gid=P[1];
    if (tid!="" && gid!="") print tid"\t"gid
  }' > "$T2G"
echo "[03] built t2g: $(wc -l < "$T2G") rows"

# 2) 由 txlist 标记需保留的 transcript 与其 gene
KEEP_G="${OUT}/keep_gene_ids.txt"
awk 'NR==FNR{want[$1]=1; next} ($1 in want){print $2}' "$TXLIST" "$T2G" \
  | sort -u > "$KEEP_G"
echo "[03] genes inferred from txlist: $(wc -l < "$KEEP_G")"

# 3) 从 GFF3 抽子集：保留 (gene in keep) + (mRNA/transcript in txlist) + (这些转录本的 exon/CDS)
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
echo "[03] subset gff3: $(wc -l < "$SUB_GFF3") lines -> $SUB_GFF3"

# 4) GFF3 -> GTF（gffread 会根据 gene/mRNA 层级写出 gene_id/transcript_id）
gffread -T -o "$CANON_GTF" "$SUB_GFF3"
echo "[03] wrote canonical GTF: $CANON_GTF"

# 5) 导出 gene_id ↔ gene symbol（从 gene 的 Name=）
$READER "$GFF3" \
| awk -F"\t" '
  $0!~/^#/ && $3=="gene" {
    gid=""; nm="";
    if (match($9,/ID=([^;]+)/,I))   gid=I[1];
    if (match($9,/Name=([^;]+)/,N)) nm=tolower(N[1]);
    if (gid!="") {
      if(nm=="") nm=gid;
      print gid"\t"nm
    }
  }' > "$ID2NAME"
echo "[03] wrote: $ID2NAME"

# 6) 生成 transcript_id ↔ gene_id ↔ symbol 三列表
awk -F'\t' '$3=="transcript"{
  if (match($9,/transcript_id "([^"]+)"/,t) && match($9,/gene_id "([^"]+)"/,g))
    print t[1]"\t"g[1]
}' "$CANON_GTF" \
| sort -u \
| awk 'NR==FNR{s[$1]=$2; next}{gid=$2; sym=(gid in s?s[gid]:gid); print $1"\t"$2"\t"sym}' \
     "$ID2NAME" - > "$TXMAP"
echo "[03] wrote: $TXMAP"

# 7) 快速验收
echo -n "[QC] transcripts kept: "
awk -F'\t' '$3=="transcript"{if (match($9,/transcript_id "([^"]+)"/,m)) print m[1] }' "$CANON_GTF" | sort -u | wc -l
echo -n "[QC] CDS lines: "
awk -F'\t' '$3=="CDS"{c++} END{print c+0}' "$CANON_GTF"




