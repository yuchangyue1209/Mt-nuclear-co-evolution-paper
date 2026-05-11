#ld test for candidates

#!/usr/bin/env bash
set -euo pipefail

REF="/work/cyu/stickleback_nuclear_only.fa"
LDX="/work/cyu/LDx.pl"
BAM="/mnt/spareHD_2/nuclear_marked_duplicates/25_RS_dedup.bam"

OUTDIR="/work/cyu/ldx_test_RS"
mkdir -p "$OUTDIR"

CHR="chrXXI"
START=3390264
END=3394799
REGION="${CHR}:${START}-${END}"

SAM_OUT="${OUTDIR}/hccsb.sam"
VCF_OUT="${OUTDIR}/hccsb.vcf"
LD_OUT="${OUTDIR}/hccsb.ld"

echo "=== RS test ==="
echo "[info] REGION = $REGION"

# 检查工具
command -v samtools >/dev/null || { echo "samtools not found"; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools not found"; exit 1; }

# 建 index
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
  echo "[info] indexing BAM"
  samtools index "$BAM"
fi

# 提取 reads
echo "[step1] extracting reads"
samtools view -h "$BAM" "$REGION" > "$SAM_OUT"

# call SNP
echo "[step2] calling SNPs"
bcftools mpileup -f "$REF" -r "$REGION" -q 30 -Q 30 -d 5000 -Ou "$BAM" \
| bcftools call -mv -Ov -o "$VCF_OUT"

# 跑 LDx
echo "[step3] running LDx"
perl "$LDX" -l 10 -q 20 -s 500 "$SAM_OUT" "$VCF_OUT" > "$LD_OUT"

echo "[done]"
echo "Output: $LD_OUT"





# 
cd /work/cyu

python3 - <<'PY'
import csv

infile = "/work/cyu/TableS2.csv"
outfile = "/work/cyu/subunit_genes_dedup.csv"
bedfile = "/work/cyu/subunit_genes_dedup.bed"

rows = []
seen = set()

with open(infile, newline='', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    reader.fieldnames = [x.strip() for x in reader.fieldnames]

    for row in reader:
        row = {k.strip(): (v.strip() if v is not None else "") for k, v in row.items()}
        if row.get("role", "") != "subunit":
            continue

        gene = row.get("stickleback_name", "")
        chrom = row.get("chr_filled", "")
        start = row.get("region_start", "")
        end = row.get("region_end", "")
        complex_ = row.get("complex", "")

        if not gene or not chrom or not start or not end:
            continue

        start = int(start)
        end = int(end)
        if end < start:
            continue

        key = (gene, chrom, start, end)
        if key in seen:
            continue
        seen.add(key)

        rows.append((gene, chrom, start, end, complex_))

rows.sort(key=lambda x: (x[1], x[2], x[3], x[0]))

with open(outfile, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene", "chr", "start", "end", "complex"])
    w.writerows(rows)

with open(bedfile, "w", newline="") as f:
    for gene, chrom, start, end, complex_ in rows:
        f.write(f"{chrom}\t{start-1}\t{end}\t{gene}\n")

print("wrote", outfile, "n=", len(rows))
print("wrote", bedfile)
PY



cat > /work/cyu/make_subunit_mpileup.sh <<'EOF'
#!/bin/bash

MANIFEST="/work/cyu/subunit_manifest.tsv"
REF="/work/cyu/stickleback_nuclear_only.fa"
OUT_DIR="/work/cyu/subunit_mpileup"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r POP BAM GENE CHR START END COMPLEX MPILEUP; do
    REGION="${CHR}:${START}-${END}"
    LOG="${LOG_DIR}/${POP}__${GENE}.log"

    if [[ -s "$MPILEUP" ]]; then
        continue
    fi

    echo "[INFO] ${POP} ${GENE} ${REGION}"

    samtools mpileup \
        -r "$REGION" \
        -f "$REF" \
        -q 20 \
        -Q 20 \
        "$BAM" > "$MPILEUP" 2> "$LOG"

    status=$?
    if [[ $status -ne 0 ]]; then
        echo "[ERROR] ${POP} ${GENE} failed; see $LOG"
        rm -f "$MPILEUP"
        continue
    fi
done
EOF

chmod +x /work/cyu/make_subunit_mpileup.sh



cat > /work/cyu/run_ldx_all_subunits.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REF="/work/cyu/stickleback_nuclear_only.fa"
LDX="/work/cyu/LDx.pl"
GENE_CSV="/work/cyu/subunit_genes_dedup.csv"
BAM_DIR="/mnt/spareHD_2/nuclear_marked_duplicates"

OUTDIR="/work/cyu/ldx_all_subunits"
SAM_DIR="${OUTDIR}/sam"
VCF_DIR="${OUTDIR}/vcf"
LD_DIR="${OUTDIR}/ld"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "$OUTDIR" "$SAM_DIR" "$VCF_DIR" "$LD_DIR" "$LOG_DIR"

command -v samtools >/dev/null || { echo "samtools not found"; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools not found"; exit 1; }
[[ -f "$REF" ]] || { echo "Reference not found: $REF"; exit 1; }
[[ -f "$LDX" ]] || { echo "LDx.pl not found: $LDX"; exit 1; }

[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

for BAM in "${BAM_DIR}"/*.bam; do
  [[ -f "${BAM}.bai" || -f "${BAM%.bam}.bai" ]] || samtools index "$BAM"
done

tail -n +2 "$GENE_CSV" | while IFS=',' read -r GENE CHR START END COMPLEX; do
  REGION="${CHR}:${START}-${END}"

  for BAM in "${BAM_DIR}"/*.bam; do
    POP=$(basename "$BAM" .bam)

    PREFIX="${POP}__${GENE}"
    SAM_OUT="${SAM_DIR}/${PREFIX}.sam"
    VCF_OUT="${VCF_DIR}/${PREFIX}.vcf"
    LD_OUT="${LD_DIR}/${PREFIX}.ld"
    LOG_OUT="${LOG_DIR}/${PREFIX}.log"

    [[ -s "$LD_OUT" ]] && continue

    echo "=== ${PREFIX} ===" > "$LOG_OUT"
    echo "[info] REGION = ${REGION}" >> "$LOG_OUT"

    echo "[step1] extracting reads" >> "$LOG_OUT"
    if ! samtools view -h "$BAM" "$REGION" > "$SAM_OUT" 2>> "$LOG_OUT"; then
      echo "[ERROR] samtools view failed" >> "$LOG_OUT"
      rm -f "$SAM_OUT" "$VCF_OUT" "$LD_OUT"
      continue
    fi

    echo "[step2] calling SNPs" >> "$LOG_OUT"
    if ! bcftools mpileup -f "$REF" -r "$REGION" -q 30 -Q 30 -d 5000 -Ou "$BAM" 2>> "$LOG_OUT" \
      | bcftools call -mv -Ov -o "$VCF_OUT" 2>> "$LOG_OUT"; then
      echo "[ERROR] bcftools call failed" >> "$LOG_OUT"
      rm -f "$VCF_OUT" "$LD_OUT"
      continue
    fi

    NSNP=$(grep -vc '^#' "$VCF_OUT" || true)
    echo "[info] SNPs in VCF = ${NSNP}" >> "$LOG_OUT"

    if [[ "${NSNP}" -eq 0 ]]; then
      echo "[warn] no SNPs found; skip LDx" >> "$LOG_OUT"
      : > "$LD_OUT"
      continue
    fi

    echo "[step3] running LDx" >> "$LOG_OUT"
    if ! perl "$LDX" -l 10 -q 20 -s 500 "$SAM_OUT" "$VCF_OUT" > "$LD_OUT" 2>> "$LOG_OUT"; then
      echo "[ERROR] LDx failed" >> "$LOG_OUT"
      rm -f "$LD_OUT"
      continue
    fi

    echo "[done] Output: $LD_OUT" >> "$LOG_OUT"
  done
done
EOF

chmod +x /work/cyu/run_ldx_all_subunits.sh


#chmod +x /work/cyu/prune_ld_single.py
#!/usr/bin/env python3
import sys
from collections import defaultdict

if len(sys.argv) != 5:
    sys.stderr.write(
        "Usage: python3 prune_ld_single.py <ld_file> <pop> <gene> <deltaAF_lookup.tsv>\n"
    )
    sys.exit(1)

ld_file = sys.argv[1]
pop = sys.argv[2]
gene = sys.argv[3]
lookup_file = sys.argv[4]

R2_THRESHOLD = 0.8

# load |deltaAF| for this pop + gene
delta = {}
with open(lookup_file) as f:
    header = next(f)
    for line in f:
        p, g, chrom, pos, daf = line.rstrip("\n").split("\t")
        if p == pop and g == gene:
            try:
                delta[pos] = abs(float(daf))
            except ValueError:
                continue

# build LD graph from strong-LD pairs
graph = defaultdict(set)
all_nodes = set()

with open(ld_file) as f:
    for line in f:
        cols = line.strip().split()
        if len(cols) < 15:
            continue

        p1 = cols[0]
        p2 = cols[1]

        try:
            r2 = float(cols[14])
        except ValueError:
            continue

        if r2 >= R2_THRESHOLD:
            graph[p1].add(p2)
            graph[p2].add(p1)
            all_nodes.add(p1)
            all_nodes.add(p2)

# find connected components = LD clusters
visited = set()
clusters = []

def dfs(start):
    comp = []
    stack = [start]
    while stack:
        node = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        comp.append(node)
        for nei in graph[node]:
            if nei not in visited:
                stack.append(nei)
    return comp

for node in all_nodes:
    if node not in visited:
        clusters.append(dfs(node))

# keep SNP with max |deltaAF| in each cluster
# tie-break: smallest genomic position
for comp in clusters:
    ranked = sorted(
        comp,
        key=lambda x: (-delta.get(x, -1.0), int(x))
    )
    best = ranked[0]
    print(best)




#淡水ld的分布
import os
import pandas as pd

ld_dir = "/work/cyu/ldx_all_subunits/ld"
freshwater_pops = {
    "BEA","BOOT","ECHO","FG","GOS","JOE","LAW","LB","LG",
    "MUC","PYE","ROB","SL","SR","SWA","THE","TL","WB","WK","WT"
}

summary = []

if not os.path.exists(ld_dir):
    print(f"Error: {ld_dir} not found.")
else:
    for f in os.listdir(ld_dir):
        if not f.endswith(".ld"):
            continue

        # 灵活解析文件名
        name_part = f.replace(".ld", "")
        if "__" in name_part:
            parts = name_part.split("__")
        elif "_" in name_part:
            parts = name_part.split("_", 1)
        else:
            parts = name_part.split()

        if len(parts) < 2:
            continue
            
        pop, gene = parts[0], parts[1]

        if pop in freshwater_pops:
            with open(os.path.join(ld_dir, f)) as ld:
                for line in ld:
                    cols = line.strip().split()
                    if len(cols) >= 15:
                        try:
                            summary.append([pop, gene, float(cols[14])])
                        except:
                            continue

if not summary:
    print("No data found. Check your file names and population IDs.")
else:
    df = pd.DataFrame(summary, columns=['Pop', 'Gene', 'R2'])
    print(f"\nTotal records: {len(df)}")
    
    # 统计不同阈值
    for t in [0.2, 0.5, 0.6, 0.7, 0.8]:
        count = len(df[df['R2'] >= t])
        print(f"R2 >= {t}: {count} pairs")

    # 保存结果
    df[df['R2'] >= 0.5].to_csv("ld_freshwater_filtered.csv", index=False)
    print("\nFile saved: ld_freshwater_filtered.csv")


#heatmap
import os
import pandas as pd

ld_dir = "/work/cyu/ldx_all_subunits/ld"

# ✅ freshwater populations（20个）
freshwater_pops = {
    "BEA","BOOT","ECHO","FG","GOS","JOE","LAW","LB","LG",
    "MUC","PYE","ROB","SL","SR","SWA","THE","TL","WB","WK","WT"
}

summary = []

# Step 1️⃣ 提取 high LD pairs
for f in os.listdir(ld_dir):
    if not f.endswith(".ld"):
        continue

    pop_gene = f.replace(".ld", "")
    pop, gene = pop_gene.split("__")

    if pop not in freshwater_pops:
        continue

    with open(os.path.join(ld_dir, f)) as ld:
        for line in ld:
            cols = line.strip().split()
            if len(cols) < 15:
                continue

            r2 = float(cols[14])
            if r2 >= 0.8:
                summary.append([pop, gene, r2])

# 转成 DataFrame
ld_df = pd.DataFrame(summary, columns=['Pop', 'Gene', 'R2'])

print(f"共找到 {len(ld_df)} 条 high LD pairs（freshwater）")

# Step 2️⃣ gene-level 汇总
ld_summary = (
    ld_df
    .groupby(['Pop', 'Gene'])
    .agg(
        n_high_LD_pairs=('R2', 'count'),
        mean_r2=('R2', 'mean'),
        max_r2=('R2', 'max')
    )
    .reset_index()
)

# Step 3️⃣ 保存
ld_summary.to_csv("ld_summary_freshwater_gene_level.csv", index=False)

print("✅ 已生成：ld_summary_freshwater_gene_level.csv")



#heatmap
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =========================
# 1. 参数
# =========================
ld_dir = "/work/cyu/ldx_all_subunits/ld"
r2_cutoff = 0.8

freshwater = {
    "BEA", "BOOT", "ECHO", "FG", "GOS", "JOE", "LAW", "LB", "LG",
    "MUC", "PYE", "ROB", "SL", "SR", "SWA", "THE", "TL", "WB", "WK", "WT"
}

# =========================
# 2. cluster 定义
# =========================
cluster_map = {
    # C1
    "SR": "C1", "TL": "C1", "WK": "C1", "LB": "C1", "MUC": "C1", "SWA": "C1",
    # C2
    "BEA": "C2", "THE": "C2",
    # C3
    "ROB": "C3", "WB": "C3", "LG": "C3", "SL": "C3",
    "LAW": "C3", "BOOT": "C3", "JOE": "C3", "FG": "C3",
    # C4
    "GOS": "C4", "PYE": "C4", "ECHO": "C4", "WT": "C4"
}

cluster_order = ["C1", "C2", "C3", "C4"]

# 固定顺序：按 cluster 排列
pop_order = [
    "LB", "MUC", "SR", "SWA", "TL", "WK",      # C1
    "BEA", "THE",                              # C2
    "BOOT", "FG", "JOE", "LAW", "LG", "ROB", "SL", "WB",  # C3
    "ECHO", "GOS", "PYE", "WT"                 # C4
]

# =========================
# 3. 读取 .ld 文件
# =========================
records = []

for f in os.listdir(ld_dir):
    if not f.endswith(".ld"):
        continue

    prefix = f.replace(".ld", "")
    if "__" not in prefix:
        continue

    pop, gene = prefix.split("__", 1)

    if pop not in freshwater:
        continue

    fpath = os.path.join(ld_dir, f)

    with open(fpath) as infile:
        for line in infile:
            cols = line.strip().split()
            if len(cols) < 15:
                continue
            try:
                r2 = float(cols[14])
            except ValueError:
                continue
            records.append([gene, pop, r2])

df = pd.DataFrame(records, columns=["Gene", "Pop", "R2"])

if df.empty:
    raise ValueError("没有读到任何 LD 数据，请检查 ld_dir。")

print("Total records:", len(df))
print("Unique genes:", df["Gene"].nunique())
print("Unique pops:", df["Pop"].nunique())

# =========================
# 4. 汇总
# =========================
summary = (
    df.groupby(["Gene", "Pop"])
      .agg(
          total_pairs=("R2", "size"),
          high_pairs=("R2", lambda x: (x >= r2_cutoff).sum()),
          LD_density=("R2", lambda x: (x >= r2_cutoff).sum() / len(x)),
          mean_r2=("R2", "mean")
      )
      .reset_index()
)

# 补全所有 gene × pop 组合
all_genes = sorted(summary["Gene"].unique())
full_index = pd.MultiIndex.from_product([all_genes, pop_order], names=["Gene", "Pop"])

summary = (
    summary.set_index(["Gene", "Pop"])
           .reindex(full_index, fill_value=0)
           .reset_index()
)

# =========================
# 5. 构建三个矩阵
# =========================
count_mat = (
    summary.pivot(index="Gene", columns="Pop", values="high_pairs")
           .fillna(0)
           .astype(int)
)

density_mat = (
    summary.pivot(index="Gene", columns="Pop", values="LD_density")
           .fillna(0)
)

meanr2_mat = (
    summary.pivot(index="Gene", columns="Pop", values="mean_r2")
           .fillna(0)
)

# gene 排序：按 high LD pairs 总数从高到低
gene_order = count_mat.sum(axis=1).sort_values(ascending=False).index.tolist()

count_mat = count_mat.loc[gene_order, pop_order]
density_mat = density_mat.loc[gene_order, pop_order]
meanr2_mat = meanr2_mat.loc[gene_order, pop_order]

# =========================
# 6. cluster 分界线和标签位置
# =========================
cluster_sizes = [6, 2, 8, 4]  # C1, C2, C3, C4
boundaries = np.cumsum(cluster_sizes)[:-1]  # [6, 8, 16]

starts = np.cumsum([0] + cluster_sizes[:-1])
midpoints = [s + n / 2 for s, n in zip(starts, cluster_sizes)]

# =========================
# 7. 统一绘图函数（无 title）
# =========================
def draw_heatmap(data, outfile, fmt, cbar_label, annot_size=7):
    plt.figure(figsize=(14, 10))

    ax = sns.heatmap(
        data,
        annot=True,
        fmt=fmt,
        cmap="YlGnBu",
        linewidths=0.5,
        cbar_kws={"label": cbar_label},
        annot_kws={"size": annot_size}
    )

    # cluster 分界线
    for b in boundaries:
        ax.axvline(b, color="black", linewidth=2)

    # 顶部 C1/C2/C3/C4 标签
    for label, mid in zip(cluster_order, midpoints):
        ax.text(
            mid,
            -0.6,
            label,
            ha="center",
            va="center",
            fontsize=14,
            fontweight="bold",
            color="black"
        )

    plt.xlabel("Population")
    plt.ylabel("Gene")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(fontsize=7)

    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ saved: {outfile}")

# =========================
# 8. 三张图
# =========================

# 图1：High LD pairs count
draw_heatmap(
    count_mat,
    "ld_count_clustered.png",
    "d",
    "Count of Pairs (R2 >= 0.8)",
    annot_size=7
)

# 图2：LD density
draw_heatmap(
    density_mat,
    "ld_density_clustered.png",
    ".2f",
    "LD density",
    annot_size=6
)

# 图3：Mean R²
draw_heatmap(
    meanr2_mat,
    "mean_r2_clustered.png",
    ".2f",
    "Mean R²",
    annot_size=6
)


#大表
import os
import pandas as pd
from collections import defaultdict

ld_dir = "/work/cyu/ldx_all_subunits/ld"

freshwater = {
    "BEA","BOOT","ECHO","FG","GOS","JOE","LAW","LB","LG",
    "MUC","PYE","ROB","SL","SR","SWA","THE","TL","WB","WK","WT"
}

R2_THRESHOLD = 0.8

# 👉 如果你有 deltaAF lookup（推荐用）
delta_file = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"

delta_df = pd.read_csv(delta_file, sep="\t")
delta_df["pos"] = delta_df["pos"].astype(str)

delta_map = {
    (row["pop"], row["pos"]): abs(row["deltaAF"])
    for _, row in delta_df.iterrows()
}

# =========================
# 输出容器
# =========================
block_records = []
leader_records = []
summary_records = []

# =========================
# 主循环
# =========================
for f in os.listdir(ld_dir):

    if not f.endswith(".ld"):
        continue

    pop, gene = f.replace(".ld", "").split("__")

    if pop not in freshwater:
        continue

    path = os.path.join(ld_dir, f)

    # ===== 构建 LD graph =====
    graph = defaultdict(set)
    nodes = set()

    with open(path) as infile:
        for line in infile:
            cols = line.strip().split()
            if len(cols) < 15:
                continue

            p1, p2 = cols[0], cols[1]

            try:
                r2 = float(cols[14])
            except:
                continue

            if r2 >= R2_THRESHOLD:
                graph[p1].add(p2)
                graph[p2].add(p1)
                nodes.update([p1, p2])

    # ===== 找 LD blocks（connected components）=====
    visited = set()
    blocks = []

    for node in nodes:
        if node in visited:
            continue

        stack = [node]
        comp = []

        while stack:
            curr = stack.pop()
            if curr in visited:
                continue

            visited.add(curr)
            comp.append(curr)
            stack.extend(graph[curr] - visited)

        blocks.append(comp)

    # ===== 每个 block 处理 =====
    for i, block in enumerate(blocks):

        # leader SNP（最大 deltaAF）
        scored = sorted(
            block,
            key=lambda x: (-delta_map.get((pop, x), -1), int(x))
        )

        leader = scored[0]

        # 记录 block
        block_records.append({
            "Pop": pop,
            "Gene": gene,
            "Block_ID": f"{pop}_{gene}_{i}",
            "Block_size": len(block),
            "Leader_SNP": leader
        })

        # 记录 leader
        leader_records.append({
            "Pop": pop,
            "Gene": gene,
            "Leader_SNP": leader,
            "Block_size": len(block)
        })

    # ===== gene-level summary =====
    summary_records.append({
        "Pop": pop,
        "Gene": gene,
        "n_blocks": len(blocks),
        "mean_block_size": sum(len(b) for b in blocks) / len(blocks) if blocks else 0,
        "max_block_size": max(len(b) for b in blocks) if blocks else 0
    })

# =========================
# 保存
# =========================
block_df = pd.DataFrame(block_records)
leader_df = pd.DataFrame(leader_records)
summary_df = pd.DataFrame(summary_records)

block_df.to_csv("ld_blocks_all.csv", index=False)
leader_df.to_csv("ld_leader_snps.csv", index=False)
summary_df.to_csv("ld_summary_blocks.csv", index=False)

print("✅ 完成")

#3，py
3.py
import os
import pandas as pd
from collections import defaultdict

ld_dir = "/work/cyu/ldx_all_subunits/ld"
R2_THRESHOLD = 0.8

freshwater = {
    "BEA","BOOT","ECHO","FG","GOS","JOE","LAW","LB","LG",
    "MUC","PYE","ROB","SL","SR","SWA","THE","TL","WB","WK","WT"
}

cluster_map = {
    "SR":"C1_AK","TL":"C1_AK","WK":"C1_AK","LB":"C1_AK","MUC":"C1_AK","SWA":"C1_AK",
    "BEA":"C2_Recent","THE":"C2_Recent",
    "ROB":"C3_MarineLike","WB":"C3_MarineLike","LG":"C3_MarineLike","SL":"C3_MarineLike",
    "LAW":"C3_MarineLike","BOOT":"C3_MarineLike","JOE":"C3_MarineLike","FG":"C3_MarineLike",
    "GOS":"C4_GOS","PYE":"C4_GOS","ECHO":"C4_GOS","WT":"C4_GOS"
}

delta_file = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
delta_df = pd.read_csv(delta_file, sep="\t", compression="gzip")
delta_df["pop"] = delta_df["pop"].astype(str).str.upper()
delta_df["pos"] = delta_df["pos"].astype(str)

delta_map = {
    (row["pop"], row["pos"]): abs(row["deltaAF"])
    for _, row in delta_df.iterrows()
}

delta_pos_set = set(zip(delta_df["pop"], delta_df["pos"]))

block_rows = []

for f in os.listdir(ld_dir):
    if not f.endswith(".ld"):
        continue

    pop, gene = f.replace(".ld", "").split("__", 1)
    pop = pop.upper()
    if pop not in freshwater:
        continue

    graph = defaultdict(set)
    nodes = set()

    with open(os.path.join(ld_dir, f)) as infile:
        for line in infile:
            cols = line.strip().split()
            if len(cols) < 15:
                continue
            p1, p2 = cols[0], cols[1]
            try:
                r2 = float(cols[14])
            except ValueError:
                continue
            if r2 >= R2_THRESHOLD:
                graph[p1].add(p2)
                graph[p2].add(p1)
                nodes.update([p1, p2])

    visited = set()
    blocks = []

    for node in nodes:
        if node in visited:
            continue
        stack = [node]
        comp = []
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            comp.append(curr)
            stack.extend(graph[curr] - visited)
        blocks.append(comp)

    for i, block in enumerate(blocks):
        block = sorted(block, key=lambda x: int(x))
        block_hits = [x for x in block if (pop, x) in delta_pos_set]

        if block_hits:
            leader = sorted(
                block_hits,
                key=lambda x: (-delta_map[(pop, x)], int(x))
            )[0]
            leader_source = "deltaAF"
            leader_deltaAF = delta_map[(pop, leader)]
            leader_has_delta = True
        else:
            leader = block[0]
            leader_source = "min_pos"
            leader_deltaAF = None
            leader_has_delta = False

        block_rows.append({
            "Pop": pop,
            "Gene": gene,
            "Cluster": cluster_map[pop],
            "Block_ID": f"{pop}_{gene}_{i}",
            "Block_size": len(block),
            "n_deltaAF_snps_in_block": len(block_hits),
            "block_has_deltaAF": len(block_hits) > 0,
            "Leader_SNP": leader,
            "Leader_source": leader_source,
            "Leader_deltaAF": leader_deltaAF,
            "Leader_has_deltaAF": leader_has_delta
        })

block_df = pd.DataFrame(block_rows)
block_df.to_csv("ld_blocks_overlap_deltaAF.csv", index=False)

print(block_df.head())
print("\nTotal blocks:", len(block_df))
print("Blocks with ≥1 deltaAF SNP:", block_df["block_has_deltaAF"].sum())








#ld汇总 freshwater
import os
import pandas as pd
from collections import defaultdict

# =========================
# 1. 参数
# =========================
ld_dir = "/work/cyu/ldx_all_subunits/ld"
delta_file = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
R2_THRESHOLD = 0.8

freshwater = {
    "BEA","BOOT","ECHO","FG","GOS","JOE","LAW","LB","LG",
    "MUC","PYE","ROB","SL","SR","SWA","THE","TL","WB","WK","WT"
}

cluster_map = {
    "SR":"C1_AK","TL":"C1_AK","WK":"C1_AK","LB":"C1_AK","MUC":"C1_AK","SWA":"C1_AK",
    "BEA":"C2_Recent","THE":"C2_Recent",
    "ROB":"C3_MarineLike","WB":"C3_MarineLike","LG":"C3_MarineLike","SL":"C3_MarineLike",
    "LAW":"C3_MarineLike","BOOT":"C3_MarineLike","JOE":"C3_MarineLike","FG":"C3_MarineLike",
    "GOS":"C4_GOS","PYE":"C4_GOS","ECHO":"C4_GOS","WT":"C4_GOS"
}

# =========================
# 2. 读入 deltaAF 大表
# =========================
df = pd.read_csv(delta_file, sep="\t", compression="gzip")
df["pop"] = df["pop"].astype(str).str.upper()
df["pos"] = df["pos"].astype(str)

# 只保留 freshwater
df = df[df["pop"].isin(freshwater)].copy()

# 为了选 leader，建 abs(deltaAF) lookup
val_map = df.set_index(["pop", "pos"])["deltaAF"].abs().to_dict()

# deltaAF universe
delta_pos_set = set(zip(df["pop"], df["pos"]))

# =========================
# 3. 容器
# =========================
block_rows = []
summary_rows = []

# pruning
redundant_pairs = set()
leader_pairs = set()

# =========================
# 4. 遍历所有 .ld 文件，构建 population-specific LD blocks
# =========================
for ld_file in os.listdir(ld_dir):
    if not ld_file.endswith(".ld"):
        continue

    prefix = ld_file.replace(".ld", "")
    if "__" not in prefix:
        continue

    pop, gene = prefix.split("__", 1)
    pop = pop.upper()

    if pop not in freshwater:
        continue

    graph = defaultdict(set)
    nodes = set()

    with open(os.path.join(ld_dir, ld_file)) as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) < 15:
                continue

            p1, p2 = cols[0], cols[1]

            try:
                r2 = float(cols[14])
            except ValueError:
                continue

            if r2 >= R2_THRESHOLD:
                graph[p1].add(p2)
                graph[p2].add(p1)
                nodes.update([p1, p2])

    visited = set()
    blocks = []

    for node in nodes:
        if node in visited:
            continue

        cluster = []
        stack = [node]

        while stack:
            curr = stack.pop()
            if curr in visited:
                continue

            visited.add(curr)
            cluster.append(curr)
            stack.extend(graph[curr] - visited)

        blocks.append(sorted(cluster, key=lambda x: int(x)))

    # =========================
    # 5. 每个 block 选 leader（只在 deltaAF universe 内）
    # =========================
    for i, block in enumerate(blocks):
        block_hits = [x for x in block if (pop, x) in delta_pos_set]

        if len(block_hits) > 0:
            leader = sorted(
                block_hits,
                key=lambda x: (-val_map[(pop, x)], int(x))
            )[0]
            leader_source = "deltaAF"
            leader_deltaAF = val_map[(pop, leader)]
            leader_has_deltaAF = True
        else:
            # block 完全不在 deltaAF universe 中，只能 fallback
            leader = block[0]
            leader_source = "min_pos"
            leader_deltaAF = None
            leader_has_deltaAF = False

        # 记录 leader
        if leader_has_deltaAF:
            leader_pairs.add((pop, leader))

        # 只有 block_hits 里多于1个时才标记冗余
        if len(block_hits) > 1:
            for snp_pos in block_hits:
                if snp_pos != leader:
                    redundant_pairs.add((pop, snp_pos))

        block_rows.append({
            "Pop": pop,
            "Gene": gene,
            "Cluster": cluster_map[pop],
            "Block_ID": f"{pop}_{gene}_{i}",
            "Block_size": len(block),
            "n_deltaAF_snps_in_block": len(block_hits),
            "block_has_deltaAF": len(block_hits) > 0,
            "Leader_SNP": leader,
            "Leader_source": leader_source,
            "Leader_deltaAF": leader_deltaAF,
            "Leader_has_deltaAF": leader_has_deltaAF
        })

    summary_rows.append({
        "Pop": pop,
        "Gene": gene,
        "Cluster": cluster_map[pop],
        "n_blocks": len(blocks),
        "total_snps_in_blocks": sum(len(b) for b in blocks) if blocks else 0,
        "mean_block_size": (sum(len(b) for b in blocks) / len(blocks)) if blocks else 0,
        "max_block_size": max(len(b) for b in blocks) if blocks else 0
    })

# =========================
# 6. block-level 大表
# =========================
block_df = pd.DataFrame(block_rows)

# 补一个 pop×gene summary，用来跑 block-level model
pop_gene_block = (
    block_df.groupby(["Pop", "Gene", "Cluster"])
    .agg(
        n_blocks=("Block_ID", "count"),
        total_snps_in_blocks=("Block_size", "sum"),
        mean_block_size=("Block_size", "mean"),
        max_block_size=("Block_size", "max"),
        prop_blocks_with_deltaAF=("block_has_deltaAF", "mean"),
        n_leader_with_deltaAF=("Leader_has_deltaAF", "sum")
    )
    .reset_index()
)

# =========================
# 7. 对 deltaAF 大表做 pruning
# =========================
df["is_ld_leader"] = df.apply(lambda r: (r["pop"], r["pos"]) in leader_pairs, axis=1)
df["is_ld_redundant"] = df.apply(lambda r: (r["pop"], r["pos"]) in redundant_pairs, axis=1)

df["deltaAF_pruned"] = df["deltaAF"]
df.loc[df["is_ld_redundant"], "deltaAF_pruned"] = pd.NA

# 严格保留版：去掉冗余 SNP
df_keep = df[~df["is_ld_redundant"]].copy()

# =========================
# 8. 保存
# =========================
block_df.to_csv("ld_blocks_big_table.csv", index=False)
print("✅ saved: ld_blocks_big_table.csv")

pop_gene_block.to_csv("ld_block_pop_gene_summary.csv", index=False)
print("✅ saved: ld_block_pop_gene_summary.csv")

df.to_csv("deltaAF_long.noAMO.ldPruned_masked.tsv.gz", sep="\t", index=False, compression="gzip")
print("✅ saved: deltaAF_long.noAMO.ldPruned_masked.tsv.gz")

df_keep.to_csv("deltaAF_long.noAMO.ldPruned_kept.tsv.gz", sep="\t", index=False, compression="gzip")
print("✅ saved: deltaAF_long.noAMO.ldPruned_kept.tsv.gz")

# =========================
# 9. 打印一些 summary
# =========================
print("\n=== Summary ===")
print("Total LD blocks:", len(block_df))
print("Blocks with >=1 deltaAF SNP:", int(block_df["block_has_deltaAF"].sum()))
print("Blocks with deltaAF-supported leader:", int(block_df["Leader_has_deltaAF"].sum()))
print("Redundant deltaAF SNPs removed:", int(df["is_ld_redundant"].sum()))
print("Leader deltaAF SNPs kept:", int(df["is_ld_leader"].sum()))
print("Rows in original deltaAF table:", len(df))
print("Rows in pruned-kept table:", len(df_keep))



#ld prune model
#!/usr/bin/env Rscript
# ============================================================
# LD-pruned version
# deltaAF ~ mtCluster(manual) + treePC1 + treePC2
# per SNP, within region
# NO permutation
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ----------------------------
# Config
# ----------------------------
DELTA_FILE   <- "/work/cyu/ldx_all_subunits/ld/deltaAF_long.noAMO.ldPruned_kept.tsv.gz"
CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
OUTDIR       <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0_ldPruned"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

safe_zread <- function(f){
  fread(cmd = paste("zcat", shQuote(f)))
}

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

lm_cluster_ftest <- function(dt){
  dt <- dt[
    is.finite(deltaAF) &
      is.finite(treePC1) &
      is.finite(treePC2) &
      !is.na(mtCluster)
  ]

  if (nrow(dt) < 6) {
    return(list(p = NA_real_, F = NA_real_, n = nrow(dt)))
  }

  if (length(unique(dt$mtCluster)) < 2) {
    return(list(p = NA_real_, F = NA_real_, n = nrow(dt)))
  }

  fit_full <- try(lm(deltaAF ~ mtCluster + treePC1 + treePC2, data = dt), silent = TRUE)
  fit_red  <- try(lm(deltaAF ~             treePC1 + treePC2, data = dt), silent = TRUE)

  if (inherits(fit_full, "try-error") || inherits(fit_red, "try-error")) {
    return(list(p = NA_real_, F = NA_real_, n = nrow(dt)))
  }

  a <- try(anova(fit_red, fit_full), silent = TRUE)
  if (inherits(a, "try-error") || nrow(a) < 2) {
    return(list(p = NA_real_, F = NA_real_, n = nrow(dt)))
  }

  list(
    p = as.numeric(a$`Pr(>F)`[2]),
    F = as.numeric(a$F[2]),
    n = nrow(dt)
  )
}









