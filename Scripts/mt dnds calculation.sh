#mt dnds calculation
TREEFILE="/mnt/spareHD_2/oxphos_gene_tree/unrooted_NJ_tree_pruned.standard.tree"
MT_GENE_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13"
MT_CAT_DIR="/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_concat/fasta"

# ================= 输出路径 =================
OUTROOT="/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/mt"
mkdir -p "$OUTROOT"

# ================= 模型列表 =================
declare -A NSMAP=(
    [M0]=0
    [M1a]=1
    [M2a]=2
    [M7]=7
    [M8]=8
)

# ================= 运行函数 =================
run_codeml() {
    local align=$1
    local tag=$2
    local outdir=$3
    mkdir -p "$outdir"
    cp "$TREEFILE" "$outdir/tree.nwk"

    for model in "${!NSMAP[@]}"; do
        local model_dir="$outdir/$model"
        mkdir -p "$model_dir"

        cat > "$model_dir/codeml.ctl" <<EOF
seqfile = $align
treefile = $outdir/tree.nwk
outfile = $model_dir/result.txt

noisy = 9
verbose = 1
runmode = 0
seqtype = 1       * codon sequences
CodonFreq = 2
clock = 0
aaDist = 0
model = 0         * same omega distribution for all branches
NSsites = ${NSMAP[$model]}
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
EOF
        (cd "$model_dir" && codeml codeml.ctl)
    done
}

# ================= gene 级别 =================
echo "[plan] Searching mtDNA gene alignments..."
mapfile -t MT_GENES < <(find "$MT_GENE_DIR" -type f -name "*.codon.fas" | sort)
echo "[plan] Found ${#MT_GENES[@]} gene alignments."

for fas in "${MT_GENES[@]}"; do
    gene=$(basename "$fas" .codon.fas)
    run_codeml "$fas" "$gene" "$OUTROOT/gene/$gene"
done

# ================= concat 级别 =================
echo "[plan] Searching mtDNA set/whole alignments..."
mapfile -t MT_CATS < <(find "$MT_CAT_DIR" -type f -name "*.codon.concat.fas" | sort)
echo "[plan] Found ${#MT_CATS[@]} set/whole alignments."

for fas in "${MT_CATS[@]}"; do
    setname=$(basename "$fas" .codon.concat.fas)
    run_codeml "$fas" "$setname" "$OUTROOT/set/$setname"
done

echo "[done] All mtDNA codeml runs finished."





#nu
#!/usr/bin/env bash
# 无 set -euo，个别失败不会中断整个批处理

# ========= 输入 =========
TREEFILE="/mnt/spareHD_2/oxphos_gene_tree/unrooted_NJ_tree_pruned.standard.tree"

# Gene-level 对齐：每个子目录一个基因，文件名为 <gene>.codon.fas
NU_GENE_DIR="/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"

# Set/whole 对齐：*codon.concat.fas（排除 mt_*）
NU_CAT_DIR="/mnt/spareHD_2/oxphos_codeml_ready/07_concat/fasta"

# ========= 输出 =========
OUTROOT="/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/nu"
mkdir -p "$OUTROOT"
LOGERR="$OUTROOT/error_runs.log"
: > "$LOGERR"   # 清空旧的错误日志

# ========= 模型（站点模型，皆为 model=0）=========
declare -A NSMAP=(
  [M0]=0
  [M1a]=1
  [M2a]=2
  [M7]=7
  [M8]=8
)
MODELS=(M0 M1a M2a M7 M8)   # 固定迭代顺序

run_codeml() {
  local align="$1"   # 多序列密码子比对
  local tag="$2"     # 标签（基因名或集合名）
  local outdir="$3"  # 输出根目录
  mkdir -p "$outdir"

  # 复制树（统一拓扑；不 prune）
  cp -f "$TREEFILE" "$outdir/tree.nwk"

  # 简单样本数检查（至少3个序列再跑）
  local ntaxa
  ntaxa=$(grep -c '^>' "$align" 2>/dev/null || echo 0)
  if [ "$ntaxa" -lt 3 ]; then
    echo "[skip] $tag : ntaxa=$ntaxa (<3)" | tee -a "$LOGERR"
    return 0
  fi

  for model in "${MODELS[@]}"; do
    local model_dir="$outdir/$model"
    mkdir -p "$model_dir"

    # 写 ctl（核密码表 icode=0；site models 由 NSsites 指定）
    cat > "$model_dir/codeml.ctl" <<EOF
seqfile = $align
treefile = $outdir/tree.nwk
outfile = $model_dir/result.txt

noisy = 3
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 0
NSsites = ${NSMAP[$model]}
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
EOF

    # 跑 codeml；不因退出码而终止批处理
    ( cd "$model_dir" && codeml codeml.ctl > run.log 2>&1 )
    # 结果存在性检查（有些“error: end of tree file”也能正常产出 result.txt）
    if [ ! -s "$model_dir/result.txt" ]; then
      echo "[WARN] no output: $tag $model" | tee -a "$LOGERR"
    else
      # 简要提示 lnL
      awk '/lnL/ && /time/ {for(i=1;i<=NF;i++) if($i=="lnL") print "[ok] '"$tag $model"': lnL="$(i+2)}' "$model_dir/result.txt" | tail -n1
    fi
  done
}

echo "==== NU gene-level ===="
# 递归找所有 *.codon.fas（每基因一个目录）
mapfile -t NU_GENE_ALN < <(find "$NU_GENE_DIR" -type f -name "*.codon.fas" | sort)
echo "[plan] genes: ${#NU_GENE_ALN[@]}"

for fas in "${NU_GENE_ALN[@]}"; do
  gene="$(basename "$fas" .codon.fas)"
  run_codeml "$fas" "$gene" "$OUTROOT/gene/$gene"
done

echo "==== NU set/whole ===="
# 只拿 nu 的 concat，排除 mt_* 前缀
mapfile -t NU_CAT_ALN < <(find "$NU_CAT_DIR" -maxdepth 1 -type f -name "*codon.concat.fas" ! -name "mt_*" | sort)
echo "[plan] sets: ${#NU_CAT_ALN[@]}"

for fas in "${NU_CAT_ALN[@]}"; do
  setname="$(basename "$fas" .codon.concat.fas)"
  run_codeml "$fas" "$setname" "$OUTROOT/set/$setname"
done

echo "[done] NU codeml runs. Error log: $LOGERR"
