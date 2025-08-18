#nu dnds
#00_call_variants_hap1.sh

#!/usr/bin/env bash
set -euo pipefail

REF="/work/cyu/stickleback_nuclear_only.fa"
BED="/work/cyu/nuOXPHOS_genes_with_complex_core.bed"
BAM_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
OUT_VCF="/mnt/spareHD_2/oxphos_codeml_ready/01_vcf_hap1"
mkdir -p "$OUT_VCF"

# 1) 先确保有 .fai
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# 2) 用参考顺序排序 BED（只需一次）
BED_SORTED="${OUT_VCF}/nuOXPHOS_genes_with_complex_core.sorted.bed"
bedtools sort -faidx "${REF}.fai" -i "$BED" | bedtools merge > "$BED_SORTED"

for BAM in "$BAM_DIR"/*_rg.bam; do
  SAMPLE=$(basename "$BAM" _rg.bam)
  echo "[00] $SAMPLE → VCF (hap1)"

  bcftools mpileup -Ou -f "$REF" -T "$BED_SORTED" -q 30 -Q 25 \
    -a FORMAT/AD,FORMAT/DP "$BAM" \
  | bcftools call -mv -Ou --ploidy 1 \
  | bcftools norm -f "$REF" -m -any -Ou \
  | bcftools sort -Oz -o "$OUT_VCF/${SAMPLE}.oxphos.hap1.vcf.gz"

  bcftools index -f "$OUT_VCF/${SAMPLE}.oxphos.hap1.vcf.gz"
done

#01
#!/usr/bin/env bash
set -euo pipefail

# ---- paths (keep new prefixes) ----
REF="/work/cyu/stickleback_nuclear_only.fa"
BED="/work/cyu/nuOXPHOS_genes_with_complex_core.bed"
BAM_DIR="/mnt/spareHD_2/nuclear_with_readgroups"
VCF_DIR="/mnt/spareHD_2/oxphos_codeml_ready/01_vcf_hap1"
MASK_DIR="/mnt/spareHD_2/oxphos_codeml_ready/02_masks"

# ---- thresholds (light mask for OXPHOS) ----
MIN_DP=6          # depth <6 -> mask
AF_LO=0.2         # 0.2 < AF < 0.8 -> mask (not fixed)
AF_HI=0.8
INDEL_PAD=3       # ±3 bp around indels -> mask

mkdir -p "$MASK_DIR"

# ensure REF.fai
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# sort/merge BED to reference order (safe even if already sorted)
BED_SORTED="${MASK_DIR}/nuOXPHOS_genes_with_complex_core.sorted.bed"
bedtools sort -faidx "${REF}.fai" -i "$BED" | bedtools merge > "$BED_SORTED"

for VCF in "$VCF_DIR"/*.oxphos.hap1.vcf.gz; do
  SAMPLE=$(basename "$VCF" .oxphos.hap1.vcf.gz)
  # prefer sorted BAM if present
  if [ -f "${BAM_DIR}/${SAMPLE}_rg.sorted.bam" ]; then
    BAM="${BAM_DIR}/${SAMPLE}_rg.sorted.bam"
  else
    BAM="${BAM_DIR}/${SAMPLE}_rg.bam"
  fi

  echo "[01] ${SAMPLE} → building hap1 mask"

  TMP_AF="$(mktemp)"
  TMP_DP="$(mktemp)"
  TMP_INDEL="$(mktemp)"

  # 1) not-fixed sites (0.2<AF<0.8) AND variant-level low depth (DP<MIN_DP)
  #    AD format expected from 00 step (REF,ALT counts)
  bcftools query -f '%CHROM\t%POS0\t%POS\t[%AD]\n' "$VCF" \
  | awk -v OFS='\t' -v mdp="$MIN_DP" -v lo="$AF_LO" -v hi="$AF_HI" '
      {
        split($4,a,","); dp=a[1]+a[2];
        if (dp==0) next;
        if (dp<mdp) {print $1,$2,$3; next}
        af=a[2]/dp;
        if (af>lo && af<hi) print $1,$2,$3
      }' > "$TMP_AF"

  # 2) genome-level low depth inside target intervals (even if not variant)
  samtools depth -b "$BED_SORTED" "$BAM" \
  | awk -v mdp="$MIN_DP" 'BEGIN{OFS="\t"} $3<mdp {print $1,$2-1,$2}' > "$TMP_DP"

  # 3) ±INDEL_PAD around indels
  bcftools view -v indels -H "$VCF" \
  | awk -v p="$INDEL_PAD" 'BEGIN{OFS="\t"}
      { s=$2-1-p; if(s<0)s=0; e=$2+p; print $1,s,e }' > "$TMP_INDEL"

  # 4) merge all -> final mask
  cat "$TMP_AF" "$TMP_DP" "$TMP_INDEL" \
  | bedtools sort -faidx "${REF}.fai" \
  | bedtools merge > "$MASK_DIR/${SAMPLE}.hap1.mask.bed"

  rm -f "$TMP_AF" "$TMP_DP" "$TMP_INDEL"
  echo "[01] ${SAMPLE} ✓ mask -> ${MASK_DIR}/${SAMPLE}.hap1.mask.bed"
done

echo "[01] All masks done."



#02
# 02_make_consensus_hap1.sh
#!/usr/bin/env bash
set -euo pipefail

# ---- paths ----
REF="/work/cyu/stickleback_nuclear_only.fa"
VCF_DIR="/mnt/spareHD_2/oxphos_codeml_ready/01_vcf_hap1"
MASK_DIR="/mnt/spareHD_2/oxphos_codeml_ready/02_masks"
CONS_DIR="/mnt/spareHD_2/oxphos_codeml_ready/03_consensus_hap1"

mkdir -p "$CONS_DIR"

# 参考必须有 .fai
[ -f "${REF}.fai" ] || samtools faidx "$REF"

for VCF in "$VCF_DIR"/*.oxphos.hap1.vcf.gz; do
  SAMPLE=$(basename "$VCF" .oxphos.hap1.vcf.gz)
  MASK="$MASK_DIR/${SAMPLE}.hap1.mask.bed"
  OUT="$CONS_DIR/${SAMPLE}.hap1.consensus.fa"

  echo "[02] ${SAMPLE} → consensus (apply mask)"

  # 小检查：mask 没有也继续，只是会少 N 掩蔽
  if [ ! -s "$MASK" ]; then
    echo "[02] WARN: mask not found for $SAMPLE -> $MASK (will proceed without -m)"
    bcftools consensus -f "$REF" -s "$SAMPLE" \
      -M N "$VCF" > "$OUT"
  else
    bcftools consensus -f "$REF" -s "$SAMPLE" \
      -m "$MASK" -M N "$VCF" > "$OUT"
  fi

  # 给共识 FASTA 建索引，后面 gffread 会更快
  samtools faidx "$OUT"
done

echo "[02] All consensus done -> $CONS_DIR"




#03
#!/usr/bin/env bash
set -euo pipefail

# ---- inputs ----
GFF3_GZ="/work/cyu/stickleback_v5.gff3.gz"
BED="/work/cyu/nuOXPHOS_genes_with_complex_core.bed"

# ---- outputs ----
OUT="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/gff3_fix"
CANON_GTF="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/nuOXPHOS_canonical.gtf"
mkdir -p "$OUT"

# 0) 从 BED 的第4列得到目标基因名（小写、去重）
BED_GENES="$OUT/bed.genes"
cut -f4 "$BED" | awk 'NF>0{print tolower($1)}' | sort -u > "$BED_GENES"

# 1) 基于 GFF3 gene 行的 Name= 匹配出 gene_id（72 个左右）
KEEP_IDS="$OUT/keep_gene_ids.txt"
ID2NAME="$OUT/id2name.tsv"
zcat "$GFF3_GZ" \
| awk -F'\t' -v L="$BED_GENES" '
  BEGIN{while((getline g<L)>0){ok[g]=1}}
  $0!~/^#/ && $3=="gene" {
    gid=""; nm="";
    if (match($9,/ID=([^;]+)/,I))   gid=I[1];
    if (match($9,/Name=([^;]+)/,N)) nm=tolower(N[1]);
    if (gid!="" && nm!="" && ok[nm]) {
      print gid > "'"$KEEP_IDS"'";
      print gid"\t"nm > "'"$ID2NAME"'";
    }
  }'
echo "[03a] genes matched: $(wc -l < "$KEEP_IDS")"

# 2) 建 transcript -> gene 映射（mRNA/transcript）
T2G="$OUT/t2g.tsv"
zcat "$GFF3_GZ" \
| awk -F'\t' '
  $0!~/^#/ && ($3=="mRNA"||$3=="transcript"){
    if (match($9,/ID=([^;]+)/,I) && match($9,/Parent=([^;]+)/,P))
      print I[1]"\t"P[1]
  }' > "$T2G"

# 3) 取子集：gene + 其mRNA + 子孙(exon/CDS)
SUB_GFF3="$OUT/nuOXPHOS_subset.gff3"
zcat "$GFF3_GZ" \
| awk -F'\t' -v KG="$KEEP_IDS" -v T2G="$T2G" '
  BEGIN{
    while((getline k<KG)>0){keepG[k]=1}
    while((getline t<T2G)>0){split(t,a,"\t"); t2g[a[1]]=a[2]}
  }
  $0~/^#/ {next}
  {
    feat=$3; attr=$9; id=""; par="";
    if (match(attr,/ID=([^;]+)/,I)) id=I[1];
    if (match(attr,/Parent=([^;]+)/,P)) par=P[1];

    keep=0;
    if (feat=="gene" && keepG[id]) keep=1;
    else if ((feat=="mRNA"||feat=="transcript") && keepG[par]) keep=1;
    else if (par!="") {
      gid=(par in t2g)? t2g[par] : "";
      if ((gid!="" && keepG[gid]) || keepG[par]) keep=1;
    }
    if (keep) print
  }' > "$SUB_GFF3"

# 4) GFF3 -> GTF
SUB_GTF="$OUT/nuOXPHOS_subset.gtf"
gffread -T -o "$SUB_GTF" "$SUB_GFF3"

# 5) 统计每个 transcript 的 CDS 总长
TX_LEN="$OUT/tx_len.tsv"
awk -F'\t' '$3=="CDS"{
  len=$5-$4+1
  if (match($9,/transcript_id "([^"]+)"/,m_tx) && match($9,/gene_id "([^"]+)"/,m_g)){
    L[m_tx[1]]+=len; G[m_tx[1]]=m_g[1]
  }
} END{
  for (tx in L) print G[tx]"\t"tx"\t"L[tx]
}' "$SUB_GTF" > "$TX_LEN"

# 6) 选“每基因最长CDS”的代表转录本
KEEP_TX="$OUT/keep_tx.txt"
awk '{g=$1;tx=$2;len=$3; if(!(g in best) || len>bestlen[g]){best[g]=tx; bestlen[g]=len}}
     END{for(g in best) print best[g]}' "$TX_LEN" > "$KEEP_TX"

# 7) 生成 canonical GTF（仅保留目标基因 + 代表转录本 的 transcript/exon/CDS）
awk -F'\t' -v KT="$KEEP_TX" -v KG="$KEEP_IDS" '
BEGIN{
  while((getline x<KT)>0){keepTx[x]=1}
  while((getline y<KG)>0){keepG[y]=1}
}
$0!~/^#/ && ($3=="transcript"||$3=="exon"||$3=="CDS"){
  attr=$9; gid=""; tid=""
  if (match(attr,/gene_id "([^"]+)"/,m_g)) gid=m_g[1]; else next
  if (!keepG[gid]) next
  if (match(attr,/transcript_id "([^"]+)"/,m_t)){tid=m_t[1]; if (!keepTx[tid]) next}
  print
}' "$SUB_GTF" > "$CANON_GTF"

# 8) 快速验收
echo -n "[03a] transcripts: "
awk -F'\t' '$3=="transcript"{if(match($9,/transcript_id "([^"]+)"/,m)) print m[1]}' "$CANON_GTF" | sort -u | wc -l
echo -n "[03a] CDS lines: "
awk -F'\t' '$3=="CDS"{c++} END{print c+0}' "$CANON_GTF"

#!/usr/bin/env bash
set -euo pipefail

GTF="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/nuOXPHOS_canonical.gtf"
CONS="/mnt/spareHD_2/oxphos_codeml_ready/03_consensus_hap1"
CDS="/mnt/spareHD_2/oxphos_codeml_ready/04_cds_by_sample"
PEP="/mnt/spareHD_2/oxphos_codeml_ready/05_pep_by_sample"
mkdir -p "$CDS" "$PEP"

for F in "$CONS"/*.hap1.consensus.fa; do
  S=$(basename "$F" .hap1.consensus.fa)
  echo "[03b] $S → CDS/PEP"
  gffread -g "$F" -S "$GTF" -x "$CDS/${S}.hap1.cds.fa" -y "$PEP/${S}.hap1.pep.fa"
done

echo "[03b] done -> $CDS ; $PEP"


#04
#!/usr/bin/env bash
set -euo pipefail

GTF="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/nuOXPHOS_canonical.gtf"
CDS_DIR="/mnt/spareHD_2/oxphos_codeml_ready/04_cds_by_sample"
PEP_DIR="/mnt/spareHD_2/oxphos_codeml_ready/05_pep_by_sample"
ALIGN_DIR="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
mkdir -p "$ALIGN_DIR"

# transcript -> gene_symbol
TXMAP="/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/tx_gene_symbol.tsv"
awk -F'\t' '$1!~/^#/ && $3=="transcript" {
  if (match($9,/transcript_id "([^"]+)"/,t) && match($9,/gene_id "([^"]+)"/,g))
    print t[1]"\t"g[1]
}' "$GTF" \
| sort -u \
| awk 'NR==FNR{sym[$1]=$2;next}{g=$2;s=(g in sym?sym[g]:g);print $1"\t"$2"\t"s}' \
     "/mnt/spareHD_2/oxphos_codeml_ready/00_gtf_prep/id2name.tsv" - > "$TXMAP"

GENE_LIST="$ALIGN_DIR/genes72.list"
cut -f3 "$TXMAP" | sort -u > "$GENE_LIST"

extract_one() {  # fasta里取一个序列
  awk -v id="$2" 'BEGIN{RS=">";FS="\n"} NR>1{hdr=$1;if(hdr==id){printf(">%s\n",$1);for(i=2;i<=NF;i++)print $i;exit}}' "$1"
}

while read -r G; do
  [[ -z "$G" ]] && continue
  outd="$ALIGN_DIR/$G"; mkdir -p "$outd"
  pep_all="$outd/$G.pep.faa"; cds_all="$outd/$G.cds.fna"
  : > "$pep_all"; : > "$cds_all"

  TX=$(awk -v g="$G" '$3==g{print $1;exit}' "$TXMAP") || true
  [[ -z "$TX" ]] && { echo "[skip] $G no TX"; rm -rf "$outd"; continue; }

  for PEP in "$PEP_DIR"/*.hap1.pep.fa; do
    S=$(basename "$PEP" .hap1.pep.fa)
    CDS="$CDS_DIR/${S}.hap1.cds.fa"
    pseq=$(extract_one "$PEP" "$TX"); [[ -z "$pseq" ]] && continue
    cseq=$(extract_one "$CDS" "$TX"); [[ -z "$cseq" ]] && continue
    echo "$pseq" | sed "1 s/^>.*/>${S}/" >> "$pep_all"
    echo "$cseq" | sed "1 s/^>.*/>${S}/" >> "$cds_all"
  done

  n=$(grep -c '^>' "$pep_all" || true)
  (( n<2 )) && { echo "[skip] $G seqs=$n"; rm -rf "$outd"; continue; }

  mafft --maxiterate 1000 --localpair "$pep_all" > "$outd/$G.pep.aln.faa"
  pal2nal.pl "$outd/$G.pep.aln.faa" "$cds_all" -output fasta -nogap > "$outd/$G.codon.fas"
  echo "[04] $G ✓ ($n seqs)"
done < "$GENE_LIST"

echo "[04] all done -> $ALIGN_DIR"

#05
#!/usr/bin/env bash
set -euo pipefail
BED="/work/cyu/nuOXPHOS_genes_with_complex_core.bed"
ALIGN_DIR="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
LIST_DIR="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists"
mkdir -p "$LIST_DIR"

GENE72="$ALIGN_DIR/genes72.list"
awk 'BEGIN{OFS="\t"}{print tolower($4),$7,$8,$9}' "$BED" | sort -u > "$LIST_DIR/bed_meta.tsv"
awk 'NR==FNR{ok[$1]=1;next} ($1 in ok)' "$GENE72" "$LIST_DIR/bed_meta.tsv" > "$LIST_DIR/bed_meta_72.tsv"

for C in I II III IV V; do
  awk -v c="$C" '$2==c{print $1}' "$LIST_DIR/bed_meta_72.tsv" | sort -u > "$LIST_DIR/complex_${C}.genes"
done
awk '$4=="contact"{print $1}'    "$LIST_DIR/bed_meta_72.tsv" | sort -u > "$LIST_DIR/contact.genes"
awk '$4=="noncontact"{print $1}' "$LIST_DIR/bed_meta_72.tsv" | sort -u > "$LIST_DIR/noncontact.genes"
awk '$3=="core"{print $1}'       "$LIST_DIR/bed_meta_72.tsv" | sort -u > "$LIST_DIR/core.genes"
awk '$3=="noncore"{print $1}'    "$LIST_DIR/bed_meta_72.tsv" | sort -u > "$LIST_DIR/noncore.genes"

echo "[05a] lists in $LIST_DIR"
wc -l "$LIST_DIR"/*.genes
 33 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/complex_I.genes
   5 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/complex_II.genes
   9 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/complex_III.genes
  12 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/complex_IV.genes
  12 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/complex_V.genes
  58 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/contact.genes
  27 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/core.genes
  13 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/noncontact.genes
  44 /mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/noncore.genes

#05b python
#!/usr/bin/env python3
import os, glob
from collections import OrderedDict

ALIGN_ROOT = "/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
LIST_DIR   = "/mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists"
OUT_DIR    = "/mnt/spareHD_2/oxphos_codeml_ready/07_concat/fasta"
os.makedirs(OUT_DIR, exist_ok=True)

def read_fasta(path):
    seqs = OrderedDict()
    if not os.path.exists(path): return seqs
    with open(path) as f:
        name=None; buf=[]
        for line in f:
            line=line.rstrip()
            if not line: continue
            if line.startswith(">"):
                if name is not None: seqs[name]="".join(buf)
                name=line[1:].strip(); buf=[]
            else:
                buf.append(line)
        if name is not None: seqs[name]="".join(buf)
    return seqs

# 收集样本名 & 各基因对齐
all_samples=set()
gene_aln={}
for gdir in sorted(glob.glob(os.path.join(ALIGN_ROOT,"*"))):
    g=os.path.basename(gdir)
    aln=os.path.join(gdir, f"{g}.codon.fas")
    if os.path.exists(aln):
        seqs=read_fasta(aln)
        if seqs:
            gene_aln[g]=seqs
            all_samples.update(seqs.keys())
all_samples=sorted(all_samples)

def concat_one(list_path):
    cat=os.path.splitext(os.path.basename(list_path))[0]
    genes=[x.strip() for x in open(list_path) if x.strip()]
    glen={}
    for g in genes:
        if g in gene_aln:
            glen[g]=len(next(iter(gene_aln[g].values())))
        else:
            print(f"[warn] missing gene alignment: {g}")
    out_fa=os.path.join(OUT_DIR, f"{cat}.codon.concat.fas")
    with open(out_fa,"w") as out:
        for s in all_samples:
            parts=[]
            for g in genes:
                L=glen.get(g,0)
                seq=gene_aln.get(g,{}).get(s)
                parts.append(seq if seq is not None else '-'*L)
            if parts:
                out.write(f">{s}\n{''.join(parts)}\n")
    print(f"[05b] {cat} -> {out_fa}")

for lst in sorted(glob.glob(os.path.join(LIST_DIR,"*.genes"))):
    concat_one(lst)


#
#!/usr/bin/env bash
set -euo pipefail
CAT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/fasta"
PAIR_DIR="/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/_pairs"
mkdir -p "$PAIR_DIR"

# 选一个存在的拼接FA作为样本来源（取目录里的第一个）
ANY_FA=$(ls "$CAT_DIR"/*.codon.concat.fas | head -n1)
[ -z "${ANY_FA:-}" ] && { echo "ERROR: 没找到 *.codon.concat.fas"; exit 1; }

# 样本名列表
grep '^>' "$ANY_FA" | sed 's/^>//' | sort > "$PAIR_DIR/samples.list"

# 生成所有无序两两配对
awk 'NR==FNR{a[++n]=$0; next} END{for(i=1;i<=n;i++)for(j=i+1;j<=n;j++)print a[i]"\t"a[j]}' \
    "$PAIR_DIR/samples.list" "$PAIR_DIR/samples.list" > "$PAIR_DIR/pairs.txt"

echo "[pairs] $(wc -l < "$PAIR_DIR/pairs.txt") pairs -> $PAIR_DIR/pairs.txt"

#gene level
#!/usr/bin/env bash
set -euo pipefail

ALIGN_DIR="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
PAIR_DIR="/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/_pairs"
OUT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/gene_level"
mkdir -p "$OUT_DIR"

CTL="$OUT_DIR/pairwise.nu.gene.ctl"
cat > "$CTL" <<'EOF'
      seqfile = pair.fas
      treefile =
      outfile = pair.out
        noisy = 9
      verbose = 0
      runmode = -2
       seqtype = 1
     CodonFreq = 2
         clock = 0
         model = 0
        NSsites = 0
      icode = 0
   fix_kappa = 0
       kappa = 2
   fix_omega = 0
       omega = 0.4
EOF

subset_two () {
  local aln="$1" a="$2" b="$3"
  python - "$aln" "$a" "$b" <<'PY'
import sys
aln,a,b=sys.argv[1:4]
need={a,b}; seq={a:"", b:""}; cur=None
for line in open(aln):
    if line.startswith('>'): cur=line[1:].strip()
    else:
        if cur in need: seq[cur]+=line.strip()
if any(len(seq[k])==0 for k in need): sys.exit(1)
print(f">{a}\n{seq[a]}\n>{b}\n{seq[b]}")
PY
}

OUT_TSV="$OUT_DIR/nu_gene_pairwise_dnds.tsv"
echo -e "gene\tA\tB\tomega\tdN\tdS" > "$OUT_TSV"

PAIRS="$PAIR_DIR/pairs.txt"
for ALN in "$ALIGN_DIR"/*/*.codon.fas; do
  G=$(basename "$(dirname "$ALN")")
  tmpd="$OUT_DIR/_tmp_${G}"; mkdir -p "$tmpd"
  while IFS=$'\t' read -r A B; do
    subset_two "$ALN" "$A" "$B" > "$tmpd/pair.fas" || continue
    ( cd "$tmpd" && codeml "$CTL" >/dev/null ) || continue
    vals=$(awk '
      /dN =/ {for(i=1;i<=NF;i++){if($i=="dN") dn=$(i+2); if($i=="dS") ds=$(i+2)}}
      /omega|w/ {
        for(i=1;i<=NF;i++){
          if($i=="omega"||$i=="w"){
            if($(i+1)=="(dN/dS)" && $(i+2)=="=") w=$(i+3);
            else if($(i+1)=="=") w=$(i+2);
          }
        }
      }
      END{
        if(w==""||w=="nan"){ if(ds>0) w=dn/ds; else w="NA" }
        if(dn=="") dn="NA"; if(ds=="") ds="NA";
        print w"\t"dn"\t"ds
      }' "$tmpd/pair.out")
    echo -e "${G}\t${A}\t${B}\t${vals}" >> "$OUT_TSV"
  done < "$PAIRS"
  rm -rf "$tmpd"
  echo "[gene] $G done"
done

echo "[gene] -> $OUT_TSV"
#set level
#!/usr/bin/env bash
set -euo pipefail

CAT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/fasta"
PAIR_DIR="/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/_pairs"
OUT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/category_level"
mkdir -p "$OUT_DIR"

CTL="$OUT_DIR/pairwise.nu.cat.ctl"
cat > "$CTL" <<'EOF'
      seqfile = pair.fas
      treefile =
      outfile = pair.out
        noisy = 9
      verbose = 0
      runmode = -2
       seqtype = 1
     CodonFreq = 2
         clock = 0
         model = 0
        NSsites = 0
      icode = 0
   fix_kappa = 0
       kappa = 2
   fix_omega = 0
       omega = 0.4
EOF

subset_two () {
  local aln="$1" a="$2" b="$3"
  python - "$aln" "$a" "$b" <<'PY'
import sys
aln,a,b=sys.argv[1:4]
need={a,b}; seq={a:"", b:""}; cur=None
for line in open(aln):
    if line.startswith('>'): cur=line[1:].strip()
    else:
        if cur in need: seq[cur]+=line.strip()
if any(len(seq[k])==0 for k in need): sys.exit(1)
print(f">{a}\n{seq[a]}\n>{b}\n{seq[b]}")
PY
}

OUT_TSV="$OUT_DIR/nu_category_pairwise_dnds.tsv"
echo -e "category\tA\tB\tomega\tdN\tdS" > "$OUT_TSV"

PAIRS="$PAIR_DIR/pairs.txt"
for ALN in "$CAT_DIR"/*.codon.concat.fas; do
  CAT=$(basename "$ALN" .codon.concat.fas)
  tmpd="$OUT_DIR/_tmp_${CAT}"; mkdir -p "$tmpd"
  while IFS=$'\t' read -r A B; do
    subset_two "$ALN" "$A" "$B" > "$tmpd/pair.fas" || continue
    ( cd "$tmpd" && codeml "$CTL" >/dev/null ) || continue
    vals=$(awk '
      /dN =/ {for(i=1;i<=NF;i++){if($i=="dN") dn=$(i+2); if($i=="dS") ds=$(i+2)}}
      /omega|w/ {
        for(i=1;i<=NF;i++){
          if($i=="omega"||$i=="w"){
            if($(i+1)=="(dN/dS)" && $(i+2)=="=") w=$(i+3);
            else if($(i+1)=="=") w=$(i+2);
          }
        }
      }
      END{
        if(w==""||w=="nan"){ if(ds>0) w=dn/ds; else w="NA" }
        if(dn=="") dn="NA"; if(ds=="") ds="NA";
        print w"\t"dn"\t"ds
      }' "$tmpd/pair.out")
    echo -e "${CAT}\t${A}\t${B}\t${vals}" >> "$OUT_TSV"
  done < "$PAIRS"
  rm -rf "$tmpd"
  echo "[cat] $CAT done"
done

echo "[cat] -> $OUT_TSV"
#filter 0.001-2
#!/usr/bin/env bash
set -euo pipefail
filter_ds () {
  local IN="$1"; local OUT="$2"; local LO="${3:-0.001}"; local HI="${4:-2}"
  awk -F'\t' -v lo="$LO" -v hi="$HI" 'BEGIN{OFS="\t"}
    NR==1{print; next}
    {
      dn=$5; ds=$6; w=$4;
      if (w=="NA" || w=="" || w=="nan") { if (ds>0) w=dn/ds; else w="NA" }
      if (ds!="" && ds!="NA" && ds>=lo && ds<=hi) { $4=w; print }
    }' "$IN" > "$OUT"
  echo "[filter] $(wc -l < "$OUT") rows -> $OUT"
}

BASE="/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72"
filter_ds "$BASE/gene_level/nu_gene_pairwise_dnds.tsv" \
          "$BASE/gene_level/nu_gene_pairwise_dnds.dsFiltered.tsv"

filter_ds "$BASE/category_level/nu_category_pairwise_dnds.tsv" \
          "$BASE/category_level/nu_category_pairwise_dnds.dsFiltered.tsv"




#whole oxphos level
#!/usr/bin/env bash
set -euo pipefail

ALIGN_ROOT="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
LIST="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/nu_all.genes"
OUT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/fasta"
OUT="$OUT_DIR/nu_all.codon.concat.fas"
mkdir -p "$OUT_DIR"

python - <<'PY'
import os, glob
from collections import OrderedDict

ALIGN_ROOT="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
LIST="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/lists/nu_all.genes"
OUT="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/fasta/nu_all.codon.concat.fas"

def read_fa(p):
    d=OrderedDict(); name=None; buf=[]
    if not os.path.exists(p): return d
    for ln in open(p):
        ln=ln.rstrip()
        if not ln: continue
        if ln.startswith(">"):
            if name is not None: d[name]="".join(buf)
            name=ln[1:].strip(); buf=[]
        else:
            buf.append(ln)
    if name is not None: d[name]="".join(buf)
    return d

# gene list
genes=[x.strip() for x in open(LIST) if x.strip()]
# collect sample IDs from all existing gene alignments
gene_aln={}
samples=set()
missing=[]
for g in genes:
    aln=os.path.join(ALIGN_ROOT,g,f"{g}.codon.fas")
    if os.path.exists(aln):
        seqs=read_fa(aln)
        if seqs:
            gene_aln[g]=seqs
            samples.update(seqs.keys())
        else:
            missing.append(g)
    else:
        missing.append(g)

samples=sorted(samples)

# determine fragment lengths per gene (use first seq length); keep only genes that exist
fraglen={}
for g,seqs in gene_aln.items():
    fraglen[g]=len(next(iter(seqs.values())))

included=[g for g in genes if g in fraglen]

# write concatenation
with open(OUT,"w") as w:
    for s in samples:
        parts=[]
        for g in included:
            L=fraglen[g]
            seq=gene_aln[g].get(s)
            parts.append(seq if seq is not None else "-"*L)
        w.write(f">{s}\n{''.join(parts)}\n")

print(f"[nu_all] wrote: {OUT}")
print(f"[nu_all] samples: {len(samples)}  genes_included: {len(included)}  genes_missing: {len(missing)}")
if missing:
    print("[nu_all] missing (no alignment found):", ", ".join(missing))
PY

# 快检
ls -lh "$OUT"
grep -c '^>' "$OUT"



#!/usr/bin/env bash
set -euo pipefail
CAT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/fasta"
ALN="$CAT_DIR/nu_all.codon.concat.fas"
PAIR_DIR="/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/_pairs"
PAIRS="$PAIR_DIR/pairs.txt"
OUT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/whole_system"
mkdir -p "$OUT_DIR"

CTL="$OUT_DIR/pairwise.nu_all.ctl"
cat > "$CTL" <<'EOF'
      seqfile = pair.fas
      treefile =
      outfile = pair.out
        noisy = 9
      verbose = 0
      runmode = -2
       seqtype = 1
     CodonFreq = 2
         clock = 0
         model = 0
        NSsites = 0
      icode = 0
   fix_kappa = 0
       kappa = 2
   fix_omega = 0
       omega = 0.4
EOF

subset_two () {
  local aln="$1" a="$2" b="$3"
  python - "$aln" "$a" "$b" <<'PY'
import sys
aln,a,b=sys.argv[1:4]
need={a,b}; seq={a:"", b:""}; cur=None
for line in open(aln):
    if line.startswith('>'): cur=line[1:].strip()
    else:
        if cur in need: seq[cur]+=line.strip()
if any(len(seq[k])==0 for k in need): sys.exit(1)
print(f">{a}\n{seq[a]}\n>{b}\n{seq[b]}")
PY
}

OUT_TSV="$OUT_DIR/nu_all_pairwise_dnds.tsv"
echo -e "category\tA\tB\tomega\tdN\tdS" > "$OUT_TSV"

tmpd="$OUT_DIR/_tmp"; mkdir -p "$tmpd"
while IFS=$'\t' read -r A B; do
  subset_two "$ALN" "$A" "$B" > "$tmpd/pair.fas" || continue
  ( cd "$tmpd" && codeml "$CTL" >/dev/null ) || continue
  vals=$(awk '
    /dN =/ {for(i=1;i<=NF;i++){if($i=="dN") dn=$(i+2); if($i=="dS") ds=$(i+2)}}
    /omega|w/ {
      for(i=1;i<=NF;i++){
        if($i=="omega"||$i=="w"){
          if($(i+1)=="(dN/dS)" && $(i+2)=="=") w=$(i+3);
          else if($(i+1)=="=") w=$(i+2);
        }
      }
    }
    END{
      if(w==""||w=="nan"){ if(ds>0) w=dn/ds; else w="NA" }
      if(dn=="") dn="NA"; if(ds=="") ds="NA";
      print w"\t"dn"\t"ds
    }' "$tmpd/pair.out")
  echo -e "nu_all\t${A}\t${B}\t${vals}" >> "$OUT_TSV"
done < "$PAIRS"
rm -rf "$tmpd"
echo "[nu_all] -> $OUT_TSV"





fix_filter () {
  local IN="$1"; local OUT="$2"; local LO="${3:-0.001}"; local HI="${4:-2}"
  awk -F'\t' -v lo="$LO" -v hi="$HI" 'BEGIN{OFS="\t"}
    NR==1{print; next}
    {
      dn=$5; ds=$6; w=$4;
      if (w=="NA" || w=="" || w=="nan") { if (ds>0) w=dn/ds; else w="NA" }
      if (ds!="" && ds!="NA" && ds>=lo && ds<=hi) { $4=w; print }
    }' "$IN" > "$OUT"
  echo "[filter] $(wc -l < "$OUT") rows -> $OUT"
}

# mt
fix_filter "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_pairwise_13/whole_system/mt_all_pairwise_dnds.tsv" \
           "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_pairwise_13/whole_system/mt_all_pairwise_dnds.dsFiltered.tsv"

# nu
fix_filter "/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/whole_system/nu_all_pairwise_dnds.tsv" \
           "/mnt/spareHD_2/oxphos_codeml_ready/08_codeml_pairwise_72/whole_system/nu_all_pairwise_dnds.dsFiltered.tsv"
