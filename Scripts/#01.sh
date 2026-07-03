#01 
MT_GFF="/home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff"
MT_GTF_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_gtf"
mkdir -p "$MT_GTF_DIR"

# 取出 13 个目标基因（名称按需要小写/大写都行）
cat > "$MT_GTF_DIR/mt13.list" <<EOF
ND1
ND2
ND3
ND4
ND4L
ND5
ND6
CYTB
COX1
COX2
COX3
ATP6
ATP8
EOF

# 先转 GTF，再按 Name/ID 过滤到 13 基因（保留 transcript/exon/CDS 即可）
gffread -T -o "$MT_GTF_DIR/mt_all.gtf" "$MT_GFF"

awk -F'\t' -v L="$MT_GTF_DIR/mt13.list" '
BEGIN{while((getline g<L)>0){ok[tolower(g)]=1}}
$1!~/^#/ && ($3=="transcript"||$3=="exon"||$3=="CDS"){
  nm=""; if (match($9,/gene_name "([^"]+)"/,m)) nm=tolower(m[1]);
  else if (match($9,/gene_id "([^"]+)"/,m2)) nm=tolower(m2[1]);
  if (nm in ok) print
}' "$MT_GTF_DIR/mt_all.gtf" > "$MT_GTF_DIR/mt13_canonical.gtf"

echo "[GTF] CDS lines: $(awk -F'\t' '$3=="CDS"{c++} END{print c+0}' "$MT_GTF_DIR/mt13_canonical.gtf")"

#02revise header
# 路径
CONS_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"     # 你已有的 mt 共识（每个样本一个 .fasta）
MT_GTF="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_gtf/mt13_canonical.gtf"  # 你的 mt 13 基因 GTF
CDS_DIR="$CONS_DIR/_mt_cds_by_sample"
PEP_DIR="$CONS_DIR/_mt_pep_by_sample"
mkdir -p "$CDS_DIR" "$PEP_DIR"

# GTF 里的染色体ID（先确认一下）
MT_CONTIG=$(awk -F'\t' '$1!~/^#/ {print $1; exit}' "$MT_GTF")
echo "[check] GTF contig = $MT_CONTIG"

# 只处理“每样本未对齐的共识”，即 *_consensus.fasta
for F in "$CONS_DIR"/*_consensus.fasta; do
  S=$(basename "$F" _consensus.fasta)
  echo "[CDS/PEP] $S"

  # 1) 临时把FASTA第一行的头改成与GTF一致（只改第一行）
  TMP="$CONS_DIR/__tmp_${S}.fa"
  sed -e "1s/^>.*/>${MT_CONTIG}/" "$F" > "$TMP"
  samtools faidx "$TMP" >/dev/null 2>&1 || true

  # 2) 抽CDS
  gffread -g "$TMP" -S "$MT_GTF" -x "$CDS_DIR/${S}.mt.cds.fa"

  # 3) 用线粒体码表(2)翻译蛋白
  seqkit translate -T 2 -w 0 "$CDS_DIR/${S}.mt.cds.fa" > "$PEP_DIR/${S}.mt.pep.fa"

  # 4) 清理
  rm -f "$TMP" "$TMP.fai"
done

# 快检：任取一个样本，应该≈13条序列、全为3的倍数、蛋白无*
S=$(basename "$(ls "$CONS_DIR"/*_consensus.fasta | head -n1)" _consensus.fasta)
echo -n "[check] CDS count: "; grep -c '^>' "$CDS_DIR/${S}.mt.cds.fa"
echo -n "[check] PEP count: "; grep -c '^>' "$PEP_DIR/${S}.mt.pep.fa"
awk '/^>/{if(NR>1){print (len%3? "BAD3":"OK")}; len=0; next}{len+=length($0)} END{print (len%3? "BAD3":"OK")}' "$CDS_DIR/${S}.mt.cds.fa" | sort | uniq -c
grep -n '\*' "$PEP_DIR/${S}.mt.pep.fa" || true

#03
SAN_CDS="$CONS_DIR/_mt_cds_by_sample_sanitized"
SAN_PEP="$CONS_DIR/_mt_pep_by_sample_sanitized"
mkdir -p "$SAN_CDS" "$SAN_PEP"

for IN in "$CDS_DIR"/*.mt.cds.fa; do
  S=$(basename "$IN" .mt.cds.fa)
  OUT="$SAN_CDS/${S}.mt.cds.fa"
  # 对每条记录：去掉尾端 1–2 个碱基，使长度%3==0
  awk '
    /^>/{
      if (seq!="") {
        r=length(seq)%3; if(r>0) seq=substr(seq,1,length(seq)-r);
        print seq
      }
      print; seq=""; next
    }
    {seq=seq $0}
    END{
      if (seq!="") {
        r=length(seq)%3; if(r>0) seq=substr(seq,1,length(seq)-r);
        print seq
      }
    }
  ' "$IN" > "$OUT"

  # 重新按 mt 表 2 翻译
  seqkit translate -T 2 -w 0 "$OUT" > "$SAN_PEP/${S}.mt.pep.fa"
done

# 复查一个样本应全 OK、且无“内部 *”
S=10_THE
awk '/^>/{if(NR>1){print (len%3? "BAD3":"OK")}; len=0; next}{len+=length($0)} END{print (len%3? "BAD3":"OK")}' \
  "$SAN_CDS/${S}.mt.cds.fa" | sort | uniq -c
awk 'BEGIN{RS=">"; ORS=""} NR>1{
  split($0,a,"\n"); hdr=a[1]; seq=""
  for(i=2;i<=length(a);i++) seq=seq a[i]
  if (seq ~ /\*/ && substr(seq,length(seq),1)!="*") print "INTERNAL_STOP\t"hdr"\n"
}' "$SAN_PEP/${S}.mt.pep.fa"


#
ALIGN_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13"
mkdir -p "$ALIGN_DIR"
GENES=(ND1 ND2 COX1 COX2 ATP8 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 CYTB)

extract_one () {
  awk -v id="$2" 'BEGIN{RS=">";FS="\n"} NR>1{hdr=$1; if(hdr==id){printf(">%s\n",id); for(i=2;i<=NF;i++)print $i}}' "$1"
}

for G in "${GENES[@]}"; do
  outd="$ALIGN_DIR/$G"; mkdir -p "$outd"
  pep_all="$outd/$G.pep.faa"; cds_all="$outd/$G.cds.fna"
  : > "$pep_all"; : > "$cds_all"

  for P in "$SAN_PEP"/*.mt.pep.fa; do
    S=$(basename "$P" .mt.pep.fa)
    C="$SAN_CDS/${S}.mt.cds.fa"
    pseq=$(extract_one "$P" "$G"); cseq=$(extract_one "$C" "$G")
    [[ -z "$pseq" || -z "$cseq" ]] && continue
    echo "$pseq" | sed "1 s/^>.*/>${S}/" >> "$pep_all"
    echo "$cseq" | sed "1 s/^>.*/>${S}/" >> "$cds_all"
  done

  n=$(grep -c '^>' "$pep_all" || true)
  (( n<2 )) && { echo "[skip] $G sequences=$n"; rm -rf "$outd"; continue; }

  mafft --maxiterate 1000 --localpair "$pep_all" > "$outd/$G.pep.aln.faa"
  pal2nal.pl "$outd/$G.pep.aln.faa" "$cds_all" -output fasta -nogap -codontable 2 > "$outd/$G.codon.fas"
  echo "[align] $G ✓ ($n seqs)"
done


#!/usr/bin/env bash
set -euo pipefail
BASE="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"
LIST_DIR="$BASE/_mt_concat/lists"
mkdir -p "$LIST_DIR"

cat > "$LIST_DIR/mt_complex_I.genes" <<EOF
ND1
ND2
ND3
ND4
ND4L
ND5
ND6
EOF
echo "CYTB" > "$LIST_DIR/mt_complex_III.genes"
printf "COX1\nCOX2\nCOX3\n" > "$LIST_DIR/mt_complex_IV.genes"
printf "ATP6\nATP8\n" > "$LIST_DIR/mt_complex_V.genes"

echo "[05a] lists -> $LIST_DIR"

#py
#!/usr/bin/env python3
# 05b_mt_concat.py
import os, glob
from collections import OrderedDict

BASE = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"
ALIGN_ROOT = f"{BASE}/_mt_gene_align_13"
LIST_DIR   = f"{BASE}/_mt_concat/lists"
OUT_DIR    = f"{BASE}/_mt_concat/fasta"
os.makedirs(OUT_DIR, exist_ok=True)

def read_fa(p):
    d=OrderedDict(); name=None; buf=[]
    if not os.path.exists(p): return d
    for line in open(p):
        line=line.rstrip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: d[name]="".join(buf)
            name=line[1:].strip(); buf=[]
        else: buf.append(line)
    if name is not None: d[name]="".join(buf)
    return d

# 收集每个基因的 codon 对齐
all_samples=set(); gene_aln={}
for gdir in sorted(glob.glob(os.path.join(ALIGN_ROOT,"*"))):
    g=os.path.basename(gdir)
    aln=os.path.join(gdir,f"{g}.codon.fas")
    if os.path.exists(aln):
        seqs=read_fa(aln)
        if seqs:
            gene_aln[g]=seqs
            all_samples.update(seqs.keys())
all_samples=sorted(all_samples)

def concat_one(lst_path):
    cat=os.path.splitext(os.path.basename(lst_path))[0]
    genes=[x.strip() for x in open(lst_path) if x.strip()]
    # 每个基因片段长度（按现有对齐第一条的长度）
    glen={}
    for g in genes:
        if g in gene_aln:
            glen[g]=len(next(iter(gene_aln[g].values())))
        else:
            print("[warn] missing gene in alignments:", g)
    out=os.path.join(OUT_DIR, f"{cat}.codon.concat.fas")
    with open(out,"w") as w:
        for s in all_samples:
            parts=[]
            for g in genes:
                L=glen.get(g,0)
                seq=gene_aln.get(g,{}).get(s)
                parts.append(seq if seq is not None else '-'*L)
            w.write(f">{s}\n{''.join(parts)}\n")
    print("[05b]", cat, "->", out)

for lst in sorted(glob.glob(os.path.join(LIST_DIR,"*.genes"))):
    concat_one(lst)


#
#!/usr/bin/env bash
set -euo pipefail
BASE="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"
SAN_PEP="$BASE/_mt_pep_by_sample_sanitized"
PAIRS="$BASE/_mt_pairs/pairs.txt"
mkdir -p "$(dirname "$PAIRS")"

ls "$SAN_PEP"/*.mt.pep.fa \
| sed 's#.*/##; s/.mt\.pep\.fa$//' \
| sort > "$BASE/_mt_pairs/samples.list"

awk 'NR==FNR{a[++n]=$0; next} END{for(i=1;i<=n;i++)for(j=i+1;j<=n;j++)print a[i]"\t"a[j]}' \
    "$BASE/_mt_pairs/samples.list" "$BASE/_mt_pairs/samples.list" > "$PAIRS"

echo "[pairs] -> $PAIRS  ($(wc -l < "$PAIRS") lines)"




#
#!/usr/bin/env bash
set -euo pipefail
BASE="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"
CAT_DIR="$BASE/_mt_concat/fasta"
OUT_DIR="$BASE/_mt_pairwise_13/category_level"
PAIRS="$BASE/_mt_pairs/pairs.txt"
mkdir -p "$OUT_DIR"

CTL="$OUT_DIR/pairwise.mt.ctl"
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
      icode = 2
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
import sys
if any(len(seq[k])==0 for k in need): sys.exit(1)
print(f">{a}\n{seq[a]}\n>{b}\n{seq[b]}")
PY
}

OUT_TSV="$OUT_DIR/mt_category_pairwise_dnds.tsv"
echo -e "category\tA\tB\tomega\tdN\tdS" > "$OUT_TSV"

for ALN in "$CAT_DIR"/mt_complex_*.codon.concat.fas; do
  CAT=$(basename "$ALN" .codon.concat.fas)
  tmpd="$OUT_DIR/_tmp_${CAT}"; mkdir -p "$tmpd"
  while IFS=$'\t' read -r A B; do
    subset_two "$ALN" "$A" "$B" > "$tmpd/pair.fas" || continue
    ( cd "$tmpd" && codeml "$CTL" >/dev/null ) || continue
    vals=$(awk '/dN =/ {for(i=1;i<=NF;i++){if($i=="dN") dn=$(i+2); if($i=="dS") ds=$(i+2)}}
                /omega \(dN\/dS\) =/ {w=$4}
                /w =/ {w=$3}
                END{if(w=="")w="NA"; if(dn=="")dn="NA"; if(ds=="")ds="NA"; print w"\t"dn"\t"ds}' "$tmpd/pair.out")
    echo -e "${CAT}\t${A}\t${B}\t${vals}" >> "$OUT_TSV"
  done < "$PAIRS"
  rm -rf "$tmpd"
done

echo "[06b] -> $OUT_TSV"
#ratio
IN=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_pairwise_13/category_level/mt_category_pairwise_dnds.tsv
OUT=${IN%.tsv}.withOmega.tsv

awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{print; next}
{
  dn=$5; ds=$6; w=$4;
  if (w=="NA" || w=="" || w=="nan") {
    if (ds>0) w = dn/ds; else w="NA";
  }
  $4=w; print
}' "$IN" > "$OUT"

echo "[fix] wrote -> $OUT"


#whole oxphos level
# 路径
BASE="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"
ALIGN_ROOT="$BASE/_mt_gene_align_13"      # 逐基因对齐在这（每个基因一个目录+*.codon.fas）
LIST_DIR="$BASE/_mt_concat/lists"
OUT_DIR="$BASE/_mt_concat/fasta"
mkdir -p "$LIST_DIR" "$OUT_DIR"

# 写 mt 全基因清单（13 个）
cat > "$LIST_DIR/mt_all.genes" <<EOF
ND1
ND2
COX1
COX2
ATP8
ATP6
COX3
ND3
ND4L
ND4
ND5
ND6
CYTB
EOF

# 立即拼接（独立的小脚本，专门做 mt_all）
python - <<'PY'
import os, glob
from collections import OrderedDict
BASE="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"
ALIGN_ROOT=f"{BASE}/_mt_gene_align_13"
LIST=f"{BASE}/_mt_concat/lists/mt_all.genes"
OUT=f"{BASE}/_mt_concat/fasta/mt_all.codon.concat.fas"

def read_fa(p):
    d=OrderedDict(); name=None; buf=[]
    if not os.path.exists(p): return d
    for line in open(p):
        line=line.rstrip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: d[name]="".join(buf)
            name=line[1:].strip(); buf=[]
        else: buf.append(line)
    if name is not None: d[name]="".join(buf)
    return d

# 收集每基因对齐
gene_aln={}
samples=set()
for gdir in sorted(glob.glob(os.path.join(ALIGN_ROOT,"*"))):
    g=os.path.basename(gdir)
    aln=os.path.join(gdir,f"{g}.codon.fas")
    if os.path.exists(aln):
        seqs=read_fa(aln)
        if seqs:
            gene_aln[g]=seqs
            samples.update(seqs.keys())
samples=sorted(samples)

# 读取 mt_all 清单
genes=[x.strip() for x in open(LIST) if x.strip()]
# 记录每基因片段长度
glen={}
for g in genes:
    if g in gene_aln:
        glen[g]=len(next(iter(gene_aln[g].values())))
    else:
        print(f"[warn] missing gene alignment: {g}")

with open(OUT,"w") as w:
    for s in samples:
        parts=[]
        for g in genes:
            L=glen.get(g,0)
            seq=gene_aln.get(g,{}).get(s)
            parts.append(seq if seq is not None else '-'*L)
        w.write(f">{s}\n{''.join(parts)}\n")
print("[mt_all] ->", OUT)
PY

# 简短检查
ls -lh "$OUT_DIR/mt_all.codon.concat.fas"
grep -c '^>' "$OUT_DIR/mt_all.codon.concat.fas"




#!/usr/bin/env bash
set -euo pipefail
BASE="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus"
ALN="$BASE/_mt_concat/fasta/mt_all.codon.concat.fas"
PAIRS="$BASE/_mt_pairs/pairs.txt"   # 你刚才做 mt 的配对文件
OUT_DIR="$BASE/_mt_pairwise_13/whole_system"
mkdir -p "$OUT_DIR"

CTL="$OUT_DIR/pairwise.mt_all.ctl"
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
      icode = 2
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

OUT_TSV="$OUT_DIR/mt_all_pairwise_dnds.tsv"
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
  echo -e "mt_all\t${A}\t${B}\t${vals}" >> "$OUT_TSV"
done < "$PAIRS"
rm -rf "$tmpd"
echo "[mt_all] -> $OUT_TSV"



#into tsv
#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models"
OUT="${ROOT}/_summary"
mkdir -p "$OUT"

FIT="${OUT}/sites_models_fits.tsv"
BEST="${OUT}/best_model_per_tag.tsv"
LRT="${OUT}/lrt_results.tsv"

echo -e "domain\tlevel\ttag\tmodel\tlnL\tnp\tAIC\tkappa\tomega\tresult_path" > "$FIT"

parse_one() {
  local f="$1"   # /.../nu/gene/atp5mf/M2a/result.txt
  local path="$f"

  # 解析 domain/level/tag/model
  local domain level tag model
  domain=$(echo "$path" | awk -F'/' '{for(i=1;i<=NF;i++) if($(i-1)=="09_codeml_sites_models"){print $(i)} }')
  level=$( echo "$path" | awk -F'/' '{for(i=1;i<=NF;i++) if($(i-1)=="09_codeml_sites_models"){print $(i+1)} }')
  tag=$(   echo "$path" | awk -F'/' '{for(i=1;i<=NF;i++) if($(i-1)=="gene" || $(i-1)=="set"){print $(i)} }')
  model=$( echo "$path" | awk -F'/' '{print $(NF-1)}')

  # 提取 lnL、np、kappa、omega
  # lnL/np 行形如： lnL(ntime: 51  np: 56):   -363.353466      +0.000000
  # 用 awk 的正则抓 np 与 lnL（最后一个数字字段可能是 +0.000000，取倒数第二个浮点数更稳）
  read -r lnL np < <(awk '
    match($0,/lnL\(ntime:[^)]*np:[[:space:]]*([0-9]+)\):[[:space:]]*(-?[0-9]+\.[0-9]+)/,m){
      printf "%s %s\n", m[2], m[1]
    }' "$f" | tail -n1)

  # kappa/omega
  kappa=$(awk '/^kappa/ {print $3}' "$f" | tail -n1)
  omega=$(awk '/^omega[[:space:]=]/ {print $3}' "$f" | tail -n1)

  # 容错：可能为空
  [[ -z "${lnL:-}" ]] && lnL="NA"
  [[ -z "${np:-}"  ]] && np="NA"
  [[ -z "${kappa:-}" ]] && kappa="NA"
  [[ -z "${omega:-}" ]] && omega="NA"

  # AIC
  if [[ "$lnL" != "NA" && "$np" != "NA" ]]; then
    AIC=$(python - <<PY
lnL=$lnL
np=int($np)
print(2*np - 2*lnL)
PY
)
  else
    AIC="NA"
  fi

  echo -e "${domain}\t${level}\t${tag}\t${model}\t${lnL}\t${np}\t${AIC}\t${kappa}\t${omega}\t${path}"
}

export -f parse_one

# 收集 result.txt
mapfile -t FILES < <(find "$ROOT" -type f -name "result.txt" | sort)
if ((${#FILES[@]}==0)); then
  echo "No result.txt found under $ROOT"
  exit 1
fi

# 逐个解析
for f in "${FILES[@]}"; do
  parse_one "$f" >> "$FIT"
done

echo "[✓] wrote $FIT"

# 选 AIC 最优模型
python - "$FIT" "$BEST" <<'PY'
import sys, csv, math
fi, fo = sys.argv[1:3]
rows = list(csv.DictReader(open(fi), delimiter='\t'))
out = open(fo,'w')
print("domain\tlevel\ttag\tbest_model\tbest_AIC\tbest_lnL\tbest_np\tbest_path", file=out)

key = lambda r:(r['domain'], r['level'], r['tag'])
from collections import defaultdict
grp=defaultdict(list)
for r in rows:
    grp[key(r)].append(r)

for (dom,lev,tag), lst in grp.items():
    best=None; bestA=math.inf
    for r in lst:
        try:
            a=float(r['AIC'])
        except:
            continue
        if a<bestA:
            bestA=a; best=r
    if best is None:
        print(f"{dom}\t{lev}\t{tag}\tNA\tNA\tNA\tNA\tNA", file=out)
    else:
        print(f"{dom}\t{lev}\t{tag}\t{best['model']}\t{best['AIC']}\t{best['lnL']}\t{best['np']}\t{best['result_path']}", file=out)
out.close()
PY

echo "[✓] wrote $BEST"

# LRT：M2a vs M1a（df=2），M8 vs M7（df=2），分别对 nu 与 mt 做 FDR
python - "$FIT" "$LRT" <<'PY'
import sys, csv, math
from collections import defaultdict

fi, fo = sys.argv[1:3]
rows = list(csv.DictReader(open(fi), delimiter='\t'))

# 按 domain/level/tag 汇总模型
grp=defaultdict(dict)
for r in rows:
    grp[(r['domain'],r['level'],r['tag'])][r['model']]=r

def get_lnL_np(d, m):
    if m not in d: return None
    r=d[m]
    try:
        lnL=float(r['lnL']); np=int(r['np'])
    except:
        return None
    return (lnL,np,r['result_path'])

tests=[]
for key, dd in grp.items():
    dom,lev,tag = key
    for pair,df in [(('M2a','M1a'),2),(('M8','M7'),2)]:
        m1=get_lnL_np(dd, pair[0])
        m0=get_lnL_np(dd, pair[1])
        if not m1 or not m0: 
            tests.append((dom,lev,tag,f"{pair[0]}>{pair[1]}","NA","NA","NA"))
            continue
        lnL1,np1,_=m1; lnL0,np0,_=m0
        if any(math.isnan(x) for x in [lnL1,lnL0]): 
            tests.append((dom,lev,tag,f"{pair[0]}>{pair[1]}","NA","NA","NA"))
            continue
        LR=2*(lnL1 - lnL0)
        if LR<0: LR=0.0
        # p = 1 - CDF_chi2(LR,df)
        try:
            import mpmath as mp
            p = 1 - mp.gammainc(df/2, 0, LR/2)/mp.gamma(df/2)
            p=float(p)
        except Exception:
            p="NA"
        tests.append((dom,lev,tag,f"{pair[0]}>{pair[1]}",f"{LR:.6g}",str(df),p))

# FDR（分别对 nu 和 mt）
def bh_fdr(ps):
    # 输入 list of (idx, p or 'NA')
    vals=[(i,p) for i,p in ps if isinstance(p,(float,int))]
    m=len(vals)
    if m==0: return {i:'NA' for i,_ in ps}
    vals_sorted=sorted(vals, key=lambda x:x[1])
    q=[None]*m
    prev=1.0
    for rank,(i,p) in enumerate(vals_sorted, start=1):
        qval=min(prev, p*m/rank)
        q[rank-1]=qval
        prev=qval
    # 回填
    out={}
    for (rank,(i,p)),qv in zip(enumerate(vals_sorted, start=1), q):
        out[i]=qv
    for i,p in ps:
        if i not in out: out[i]='NA'
    return out

# 拆分域
nu_idx=[]; mt_idx=[]
for idx,(dom,lev,tag,cmp,LR,df,p) in enumerate(tests):
    if p=="NA": 
        pass
    else:
        p=float(p)
    if dom=="nu": nu_idx.append((idx,p))
    elif dom=="mt": mt_idx.append((idx,p))

nu_fdr=bh_fdr(nu_idx)
mt_fdr=bh_fdr(mt_idx)

# 输出
with open(fo,'w') as out:
    print("domain\tlevel\ttag\tcontrast\tLR\tdf\tp\tFDR", file=out)
    for idx,(dom,lev,tag,cmp,LR,df,p) in enumerate(tests):
        FDR = nu_fdr[idx] if dom=="nu" else (mt_fdr[idx] if dom=="mt" else 'NA')
        print(f"{dom}\t{lev}\t{tag}\t{cmp}\t{LR}\t{df}\t{p}\t{FDR}", file=out)
PY

echo "[✓] wrote $LRT"
echo "Done."

#py for dn ds omega
#!/usr/bin/env python3
import os, glob, re

OUTROOT = "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models"
SUMMARY_DIR = os.path.join(OUTROOT, "_summary")
os.makedirs(SUMMARY_DIR, exist_ok=True)
OUT_TSV = os.path.join(SUMMARY_DIR, "dnds_from_mlc.tsv")

DOMAINS = ["nu", "mt"]
LEVELS  = ["gene", "set"]
MODELS  = ["M0", "M1a", "M2a", "M7", "M8"]

# --- 正则准备（适配两种 lnL 打印风格）
RX_LNL_EQ     = re.compile(r"lnL\s*=\s*([\-0-9.eE]+)")
RX_LNL_NP     = re.compile(r"lnL\([^)]*np:\s*([0-9]+)\)\s*:\s*([\-0-9.eE]+)")
RX_KAPPA      = re.compile(r"kappa.*=\s*([0-9.]+)")
RX_TL_DN      = re.compile(r"tree length for dN:\s*([\-0-9.eE]+)")
RX_TL_DS      = re.compile(r"tree length for dS:\s*([\-0-9.eE]+)")

def parse_result(path):
    """返回 dict: {'lnL','np','kappa','dN_tree','dS_tree','omega_tree'}，若缺失则 'NA'。"""
    try:
        txt = open(path, "r", errors="ignore").read()
    except Exception:
        return dict(lnL="NA", np="NA", kappa="NA", dN_tree="NA", dS_tree="NA", omega_tree="NA")

    lnL, np = "NA", "NA"

    # 1) lnL 与 np
    m2 = RX_LNL_NP.search(txt)
    if m2:
        np = m2.group(1)
        lnL = m2.group(2)
    else:
        m1 = RX_LNL_EQ.search(txt)
        if m1:
            lnL = m1.group(1)

    # 2) kappa
    kappa = "NA"
    m = RX_KAPPA.search(txt)
    if m:
        kappa = m.group(1)

    # 3) 树的 dN/dS 长度
    def _get(rx):
        m = rx.search(txt)
        return m.group(1) if m else "NA"

    dN_tree = _get(RX_TL_DN)
    dS_tree = _get(RX_TL_DS)

    # 4) omega_tree = dN/dS
    omega_tree = "NA"
    try:
        if dN_tree != "NA" and dS_tree != "NA":
            dn = float(dN_tree); ds = float(dS_tree)
            if ds > 0:
                omega_tree = f"{dn/ds:.6g}"
            else:
                omega_tree = "NA"
    except Exception:
        omega_tree = "NA"

    return dict(lnL=lnL, np=np, kappa=kappa,
                dN_tree=dN_tree, dS_tree=dS_tree, omega_tree=omega_tree)

rows = []
for dom in DOMAINS:
    for lvl in LEVELS:
        for model in MODELS:
            pattern = os.path.join(OUTROOT, dom, lvl, "*", model, "result.txt")
            for res in sorted(glob.glob(pattern)):
                tag = os.path.basename(os.path.dirname(os.path.dirname(res)))  # gene or set name
                stat = parse_result(res)
                rows.append([
                    dom, lvl, tag, model,
                    stat["lnL"], stat["np"], stat["kappa"],
                    stat["dN_tree"], stat["dS_tree"], stat["omega_tree"],
                    res
                ])

# 写出
with open(OUT_TSV, "w") as out:
    out.write("\t".join([
        "domain","level","tag","model",
        "lnL","np","kappa",
        "dN_tree_len","dS_tree_len","omega_tree",
        "result_path"
    ]) + "\n")
    for r in rows:
        out.write("\t".join(r) + "\n")

print(f"[done] Wrote {OUT_TSV} with {len(rows)} lines.")

