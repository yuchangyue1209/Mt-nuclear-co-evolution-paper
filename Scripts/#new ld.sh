#new ld
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


#exclude amo recent
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



cat > /work/cyu/ldx_all_subunits/ld_prune_new_delta_noAMO_noLB_keepMarine_R2_0.2.py <<'PY'
#!/usr/bin/env python3

import os
import pandas as pd
from collections import defaultdict

# =========================
# Input / output
# =========================
ld_dir = "/work/cyu/ldx_all_subunits/ld"

delta_file = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_noAMO_noLB/deltaAF_long.noAMO_noLB.tsv.gz"

out_dir = "/work/cyu/ldx_all_subunits/ld"
os.makedirs(out_dir, exist_ok=True)

# LD threshold
R2_THRESHOLD = 0.2

# Keep marine + AK/BC freshwater
# Exclude AMO, LB, PACH, FRED, SC, CH
keep_pops = {
    "RS", "SAY",

    "FG", "LG", "SL", "SR",
    "TL", "WB", "WK", "WT",

    "BEA", "BOOT", "ECHO",
    "GOS", "JOE", "LAW",
    "MUC", "PYE", "ROB",
    "SWA", "THE"
}

cluster_map = {
    "RS": "C3",
    "SAY": "C3",

    "SR": "C1",
    "TL": "C1",
    "WK": "C1",
    "MUC": "C1",
    "SWA": "C1",

    "BEA": "C2",
    "THE": "C2",

    "ROB": "C3",
    "WB": "C3",
    "LG": "C3",
    "SL": "C3",
    "LAW": "C3",
    "BOOT": "C3",
    "JOE": "C3",
    "FG": "C3",

    "GOS": "C4",
    "PYE": "C4",
    "ECHO": "C4",
    "WT": "C4"
}

# =========================
# Read new deltaAF
# =========================
df = pd.read_csv(delta_file, sep="\t", compression="gzip")

df["pop"] = df["pop"].astype(str).str.upper()
df["gene"] = df["gene"].astype(str).str.lower()
df["pos"] = df["pos"].astype(str)

df = df[df["pop"].isin(keep_pops)].copy()

val_map = df.set_index(["pop", "gene", "pos"])["deltaAF"].abs().to_dict()
delta_key_set = set(zip(df["pop"], df["gene"], df["pos"]))

# =========================
# Containers
# =========================
block_rows = []
summary_rows = []

leader_keys = set()
redundant_keys = set()

# =========================
# Build LD blocks from existing .ld files
# =========================
for ld_file in os.listdir(ld_dir):

    if not ld_file.endswith(".ld"):
        continue

    prefix = ld_file.replace(".ld", "")
    if "__" not in prefix:
        continue

    pop, gene = prefix.split("__", 1)
    pop = pop.upper()
    gene = gene.lower()

    if pop not in keep_pops:
        continue

    path = os.path.join(ld_dir, ld_file)

    graph = defaultdict(set)
    nodes = set()

    with open(path) as f:
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

        comp = []
        stack = [node]

        while stack:
            curr = stack.pop()
            if curr in visited:
                continue

            visited.add(curr)
            comp.append(curr)
            stack.extend(graph[curr] - visited)

        blocks.append(sorted(comp, key=lambda x: int(x)))

    for i, block in enumerate(blocks):

        block_hits = [
            pos for pos in block
            if (pop, gene, pos) in delta_key_set
        ]

        if len(block_hits) > 0:
            leader = sorted(
                block_hits,
                key=lambda x: (-val_map[(pop, gene, x)], int(x))
            )[0]

            leader_source = "deltaAF"
            leader_deltaAF = val_map[(pop, gene, leader)]
            leader_has_deltaAF = True

            leader_keys.add((pop, gene, leader))

            if len(block_hits) > 1:
                for pos in block_hits:
                    if pos != leader:
                        redundant_keys.add((pop, gene, pos))

        else:
            leader = block[0] if len(block) > 0 else None
            leader_source = "min_pos"
            leader_deltaAF = None
            leader_has_deltaAF = False

        block_rows.append({
            "Pop": pop,
            "Gene": gene,
            "Cluster": cluster_map.get(pop, "Other"),
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
        "Cluster": cluster_map.get(pop, "Other"),
        "n_blocks": len(blocks),
        "total_snps_in_blocks": sum(len(b) for b in blocks) if blocks else 0,
        "mean_block_size": (sum(len(b) for b in blocks) / len(blocks)) if blocks else 0,
        "max_block_size": max(len(b) for b in blocks) if blocks else 0
    })

# =========================
# Save LD block tables
# =========================
block_df = pd.DataFrame(block_rows)
summary_df = pd.DataFrame(summary_rows)

if not block_df.empty:
    pop_gene_block = (
        block_df.groupby(["Pop", "Gene", "Cluster"], dropna=False)
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
else:
    pop_gene_block = pd.DataFrame()

block_out = os.path.join(out_dir, "ld_blocks_big_table.newDelta.noAMO_noLB.keepMarine.R2_0.2.csv")
summary_out = os.path.join(out_dir, "ld_block_pop_gene_summary.newDelta.noAMO_noLB.keepMarine.R2_0.2.csv")

block_df.to_csv(block_out, index=False)
pop_gene_block.to_csv(summary_out, index=False)

# =========================
# Prune deltaAF table
# =========================
df["is_ld_leader"] = df.apply(
    lambda r: (r["pop"], r["gene"], str(r["pos"])) in leader_keys,
    axis=1
)

df["is_ld_redundant"] = df.apply(
    lambda r: (r["pop"], r["gene"], str(r["pos"])) in redundant_keys,
    axis=1
)

df["deltaAF_pruned"] = df["deltaAF"]
df.loc[df["is_ld_redundant"], "deltaAF_pruned"] = pd.NA

df_keep = df[~df["is_ld_redundant"]].copy()

masked_out = os.path.join(out_dir, "deltaAF_long.noAMO_noLB.ldPruned_masked.keepMarine.R2_0.2.tsv.gz")
kept_out = os.path.join(out_dir, "deltaAF_long.noAMO_noLB.ldPruned_kept.keepMarine.R2_0.2.tsv.gz")

df.to_csv(masked_out, sep="\t", index=False, compression="gzip")
df_keep.to_csv(kept_out, sep="\t", index=False, compression="gzip")

# =========================
# Summary
# =========================
print("\n=== LD pruning summary ===")
print("R2 threshold:", R2_THRESHOLD)
print("Input deltaAF rows:", len(df))
print("Rows after pruning:", len(df_keep))
print("Total LD blocks:", len(block_df))
print("Blocks with >=1 deltaAF SNP:", int(block_df["block_has_deltaAF"].sum()) if len(block_df) else 0)
print("Blocks with deltaAF-supported leader:", int(block_df["Leader_has_deltaAF"].sum()) if len(block_df) else 0)
print("Redundant deltaAF SNPs removed:", int(df["is_ld_redundant"].sum()))
print("Leader deltaAF SNPs kept:", int(df["is_ld_leader"].sum()))

print("\n[OK] wrote:")
print(block_out)
print(summary_out)
print(masked_out)
print(kept_out)
PY

python /work/cyu/ldx_all_subunits/ld_prune_new_delta_noAMO_noLB_keepMarine_R2_0.2.py




#heatmap
cat > /work/cyu/ldx_all_subunits/plot_LD_heatmap_R2_0.2_keepMarine_bigFont.py <<'PY'
#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =========================
# Parameters
# =========================
ld_dir = "/work/cyu/ldx_all_subunits/ld"
out_dir = "/work/cyu/ldx_all_subunits/ld"
os.makedirs(out_dir, exist_ok=True)

r2_cutoff = 0.2

keep_pops = {
    # Marine
    "RS", "SAY",

    # Alaska freshwater
    "FG", "LG", "SL", "SR",
    "TL", "WB", "WK", "WT",

    # BC freshwater
    "BEA", "BOOT", "ECHO",
    "GOS", "JOE", "LAW",
    "MUC", "PYE", "ROB",
    "SWA", "THE"
}

cluster_map = {
    # C1
    "SR": "C1",
    "TL": "C1",
    "WK": "C1",
    "MUC": "C1",
    "SWA": "C1",

    # C2
    "BEA": "C2",
    "THE": "C2",

    # C3, including marine RS/SAY
    "RS": "C3",
    "SAY": "C3",
    "ROB": "C3",
    "WB": "C3",
    "LG": "C3",
    "SL": "C3",
    "LAW": "C3",
    "BOOT": "C3",
    "JOE": "C3",
    "FG": "C3",

    # C4
    "GOS": "C4",
    "PYE": "C4",
    "ECHO": "C4",
    "WT": "C4"
}

cluster_order = ["C1", "C2", "C3", "C4"]

pop_order = [
    # C1
    "MUC", "SR", "SWA", "TL", "WK",

    # C2
    "BEA", "THE",

    # C3
    "RS", "SAY", "BOOT", "FG", "JOE", "LAW", "LG", "ROB", "SL", "WB",

    # C4
    "ECHO", "GOS", "PYE", "WT"
]

# =========================
# Font sizes
# =========================
GENE_FONT_SIZE = 9
POP_FONT_SIZE = 12
NUMBER_FONT_SIZE = 8
AXIS_LABEL_SIZE = 15
CLUSTER_FONT_SIZE = 18
CBAR_LABEL_SIZE = 13
CBAR_TICK_SIZE = 11

# =========================
# Read LD files
# =========================
records = []
files_read = []

for f in os.listdir(ld_dir):
    if not f.endswith(".ld"):
        continue

    prefix = f.replace(".ld", "")

    if "__" not in prefix:
        continue

    pop, gene = prefix.split("__", 1)
    pop = pop.upper()
    gene = gene.lower()

    if pop not in keep_pops:
        continue

    files_read.append(f)

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
    raise ValueError("No LD records found. Check ld_dir and file names.")

# =========================
# Debug checks
# =========================
print("\n================ LD input check ================")
print("Total .ld files read:", len(files_read))
print("Total LD pairs:", len(df))
print("Unique genes:", df["Gene"].nunique())
print("Unique populations:", df["Pop"].nunique())

print("\nPopulations read:")
print(sorted(df["Pop"].unique()))

print("\nLD pair counts by population:")
print(df["Pop"].value_counts().sort_index())

if "RS" in df["Pop"].unique():
    rs_df = df[df["Pop"] == "RS"]
    print("\nRS check:")
    print("RS total LD pairs:", len(rs_df))
    print("RS pairs with R2 >= 0.2:", int((rs_df["R2"] >= r2_cutoff).sum()))
    print("RS R2 summary:")
    print(rs_df["R2"].describe())
else:
    print("\nWARNING: RS was not read into df.")

print("================================================\n")

# =========================
# Summarize per gene x pop
# =========================
summary = (
    df.groupby(["Gene", "Pop"])
      .agg(
          total_pairs=("R2", "size"),
          high_pairs=("R2", lambda x: (x >= r2_cutoff).sum()),
          LD_density=("R2", lambda x: (x >= r2_cutoff).sum() / len(x)),
          mean_r2=("R2", "mean"),
          max_r2=("R2", "max")
      )
      .reset_index()
)

summary_out = os.path.join(
    out_dir,
    "ld_summary_gene_pop_keepMarine_R2_0.2.bigFont.csv"
)
summary.to_csv(summary_out, index=False)

print("[OK] wrote:", summary_out)

# =========================
# Complete gene x pop matrix
# =========================
all_genes = sorted(summary["Gene"].unique())

full_index = pd.MultiIndex.from_product(
    [all_genes, pop_order],
    names=["Gene", "Pop"]
)

summary_full = (
    summary.set_index(["Gene", "Pop"])
           .reindex(full_index, fill_value=0)
           .reset_index()
)

# =========================
# Matrices
# =========================
count_mat = (
    summary_full
    .pivot(index="Gene", columns="Pop", values="high_pairs")
    .fillna(0)
    .astype(int)
)

density_mat = (
    summary_full
    .pivot(index="Gene", columns="Pop", values="LD_density")
    .fillna(0)
)

meanr2_mat = (
    summary_full
    .pivot(index="Gene", columns="Pop", values="mean_r2")
    .fillna(0)
)

maxr2_mat = (
    summary_full
    .pivot(index="Gene", columns="Pop", values="max_r2")
    .fillna(0)
)

# Sort genes by total high-LD pair count
gene_order = count_mat.sum(axis=1).sort_values(ascending=False).index.tolist()

count_mat = count_mat.loc[gene_order, pop_order]
density_mat = density_mat.loc[gene_order, pop_order]
meanr2_mat = meanr2_mat.loc[gene_order, pop_order]
maxr2_mat = maxr2_mat.loc[gene_order, pop_order]

# Save matrices
count_mat.to_csv(os.path.join(out_dir, "FigS3_LD_count_matrix_keepMarine_R2_0.2.bigFont.csv"))
density_mat.to_csv(os.path.join(out_dir, "FigS3_LD_density_matrix_keepMarine_R2_0.2.bigFont.csv"))
meanr2_mat.to_csv(os.path.join(out_dir, "FigS3_mean_R2_matrix_keepMarine.bigFont.csv"))
maxr2_mat.to_csv(os.path.join(out_dir, "FigS3_max_R2_matrix_keepMarine.bigFont.csv"))

# =========================
# Cluster boundaries
# =========================
cluster_sizes = [
    5,   # C1
    2,   # C2
    10,  # C3, including RS/SAY
    4    # C4
]

boundaries = np.cumsum(cluster_sizes)[:-1]
starts = np.cumsum([0] + cluster_sizes[:-1])
midpoints = [s + n / 2 for s, n in zip(starts, cluster_sizes)]

# =========================
# Plot function
# =========================
def draw_heatmap(data, outfile, fmt, cbar_label, annot=False, annot_size=NUMBER_FONT_SIZE):
    n_genes = data.shape[0]
    n_pops = data.shape[1]

    # Bigger figure to show all gene names and numbers clearly
    fig_height = max(16, n_genes * 0.34)
    fig_width = max(18, n_pops * 0.72)

    plt.figure(figsize=(fig_width, fig_height))

    ax = sns.heatmap(
        data,
        annot=annot,
        fmt=fmt,
        cmap="YlGnBu",
        linewidths=0.35,
        linecolor="white",
        cbar_kws={"label": cbar_label},
        annot_kws={
            "size": annot_size,
            "weight": "bold"
        }
    )

    # Colorbar font
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
    cbar.set_label(cbar_label, fontsize=CBAR_LABEL_SIZE)

    # Cluster boundaries
    for b in boundaries:
        ax.axvline(b, color="black", linewidth=2.2)

    # Cluster labels
    for label, mid in zip(cluster_order, midpoints):
        ax.text(
            mid,
            -0.75,
            label,
            ha="center",
            va="center",
            fontsize=CLUSTER_FONT_SIZE,
            fontweight="bold",
            color="black"
        )

    plt.xlabel("Population", fontsize=AXIS_LABEL_SIZE)
    plt.ylabel("Gene", fontsize=AXIS_LABEL_SIZE)

    # Force all tick labels
    ax.set_xticks(np.arange(n_pops) + 0.5)
    ax.set_xticklabels(
        data.columns.tolist(),
        rotation=45,
        ha="right",
        fontsize=POP_FONT_SIZE,
        fontweight="bold"
    )

    ax.set_yticks(np.arange(n_genes) + 0.5)
    ax.set_yticklabels(
        data.index.tolist(),
        rotation=0,
        fontsize=GENE_FONT_SIZE,
        fontstyle="italic"
    )

    # Make tick marks cleaner
    ax.tick_params(axis="both", length=0)

    # Avoid clipping gene names and top cluster labels
    plt.subplots_adjust(left=0.25, right=0.98, top=0.96, bottom=0.08)

    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()

    print("[OK] saved:", outfile)

# =========================
# Save heatmaps
# =========================
draw_heatmap(
    count_mat,
    os.path.join(out_dir, "FigS3_LD_count_keepMarine_R2_0.2_bigFont_geneItalic.png"),
    "d",
    "Number of SNP pairs with R² ≥ 0.2",
    annot=True,
    annot_size=NUMBER_FONT_SIZE
)

draw_heatmap(
    density_mat,
    os.path.join(out_dir, "FigS3_LD_density_keepMarine_R2_0.2_bigFont_geneItalic.png"),
    ".2f",
    "LD density (proportion of pairs with R² ≥ 0.2)",
    annot=False
)

draw_heatmap(
    meanr2_mat,
    os.path.join(out_dir, "FigS3_mean_R2_keepMarine_bigFont_geneItalic.png"),
    ".2f",
    "Mean R²",
    annot=False
)

draw_heatmap(
    maxr2_mat,
    os.path.join(out_dir, "FigS3_max_R2_keepMarine_bigFont_geneItalic.png"),
    ".2f",
    "Maximum R²",
    annot=False
)

print("\nDONE")
PY

python /work/cyu/ldx_all_subunits/plot_LD_heatmap_R2_0.2_keepMarine_bigFont.py


#SELECTION LD
library(data.table)
library(ggplot2)

# =========================
# Input
# =========================
infile <- "/work/cyu/ldx_all_subunits/ld/ld_summary_gene_pop_keepMarine_R2_0.2.bigFont.csv"

outdir <- "/work/cyu/ldx_all_subunits/ld/marine_vs_freshwater_LD"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dat <- fread(infile)

# Check columns
print(colnames(dat))
head(dat)

# =========================
# Standardize names
# =========================
setnames(dat, old = c("Gene", "Pop"), new = c("gene", "pop"), skip_absent = TRUE)

dat[, gene := tolower(gene)]
dat[, pop := toupper(pop)]

marine_pops <- c("RS", "SAY")

freshwater_pops <- c(
  "FG", "LG", "SL", "SR", "TL", "WB", "WK", "WT",
  "BEA", "BOOT", "ECHO", "GOS", "JOE", "LAW",
  "MUC", "PYE", "ROB", "SWA", "THE"
)

keep_pops <- c(marine_pops, freshwater_pops)

dat <- dat[pop %in% keep_pops]

# Habitat group
dat[, habitat := ifelse(pop %in% marine_pops, "Marine", "Freshwater")]
dat[, habitat := factor(habitat, levels = c("Marine", "Freshwater"))]

# Region / cluster labels, optional
dat[, region := fifelse(pop %in% marine_pops, "Marine",
                 fifelse(pop %in% c("FG","LG","SL","SR","TL","WB","WK","WT"), "AK_FW",
                 fifelse(pop %in% c("BEA","BOOT","ECHO","GOS","JOE","LAW","MUC","PYE","ROB","SWA","THE"), "BC_FW", NA_character_)))]

dat[, mtCluster := fifelse(pop %in% c("SR","TL","WK","MUC","SWA"), "C1",
                    fifelse(pop %in% c("BEA","THE"), "C2",
                    fifelse(pop %in% c("RS","SAY","ROB","WB","LG","SL","LAW","BOOT","JOE","FG"), "C3",
                    fifelse(pop %in% c("GOS","PYE","ECHO","WT"), "C4", NA_character_))))]

# =========================
# Complete missing gene-pop combinations as 0
# Important because missing combinations mean no LD pairs recorded
# =========================
all_genes <- sort(unique(dat$gene))
full_grid <- CJ(gene = all_genes, pop = keep_pops)

dat_full <- merge(full_grid, dat, by = c("gene", "pop"), all.x = TRUE)

dat_full[is.na(total_pairs), total_pairs := 0]
dat_full[is.na(high_pairs), high_pairs := 0]
dat_full[is.na(LD_density), LD_density := 0]
dat_full[is.na(mean_r2), mean_r2 := 0]
dat_full[is.na(max_r2), max_r2 := 0]

dat_full[, habitat := ifelse(pop %in% marine_pops, "Marine", "Freshwater")]
dat_full[, habitat := factor(habitat, levels = c("Marine", "Freshwater"))]

dat_full[, region := fifelse(pop %in% marine_pops, "Marine",
                      fifelse(pop %in% c("FG","LG","SL","SR","TL","WB","WK","WT"), "AK_FW",
                      fifelse(pop %in% c("BEA","BOOT","ECHO","GOS","JOE","LAW","MUC","PYE","ROB","SWA","THE"), "BC_FW", NA_character_)))]

dat_full[, mtCluster := fifelse(pop %in% c("SR","TL","WK","MUC","SWA"), "C1",
                         fifelse(pop %in% c("BEA","THE"), "C2",
                         fifelse(pop %in% c("RS","SAY","ROB","WB","LG","SL","LAW","BOOT","JOE","FG"), "C3",
                         fifelse(pop %in% c("GOS","PYE","ECHO","WT"), "C4", NA_character_))))]

# =========================
# Population-level summary
# =========================
pop_sum <- dat_full[, .(
  n_genes = .N,
  total_LD_pairs = sum(total_pairs),
  total_high_pairs = sum(high_pairs),
  mean_LD_density = mean(LD_density),
  median_LD_density = median(LD_density),
  mean_mean_r2 = mean(mean_r2),
  mean_max_r2 = mean(max_r2)
), by = .(pop, habitat, region, mtCluster)]

fwrite(pop_sum, file.path(outdir, "population_level_LD_summary.tsv"), sep = "\t")

print(pop_sum[order(habitat, -total_high_pairs)])

# =========================
# Habitat-level summary
# =========================
habitat_sum <- dat_full[, .(
  n_gene_pop = .N,
  n_genes = uniqueN(gene),
  n_pops = uniqueN(pop),
  total_LD_pairs = sum(total_pairs),
  total_high_pairs = sum(high_pairs),
  mean_LD_density = mean(LD_density),
  median_LD_density = median(LD_density),
  mean_mean_r2 = mean(mean_r2),
  mean_max_r2 = mean(max_r2)
), by = habitat]

fwrite(habitat_sum, file.path(outdir, "habitat_level_LD_summary.tsv"), sep = "\t")

print(habitat_sum)

# =========================
# Gene-level Marine vs Freshwater comparison
# For each gene:
# mean Marine LD_density across RS/SAY
# mean Freshwater LD_density across 19 FW populations
# This is the cleanest comparison.
# =========================
gene_habitat <- dat_full[, .(
  mean_LD_density = mean(LD_density),
  median_LD_density = median(LD_density),
  mean_high_pairs = mean(high_pairs),
  sum_high_pairs = sum(high_pairs),
  mean_mean_r2 = mean(mean_r2),
  mean_max_r2 = mean(max_r2)
), by = .(gene, habitat)]

gene_wide <- dcast(
  gene_habitat,
  gene ~ habitat,
  value.var = c(
    "mean_LD_density",
    "median_LD_density",
    "mean_high_pairs",
    "sum_high_pairs",
    "mean_mean_r2",
    "mean_max_r2"
  )
)

# Difference: Marine - Freshwater
gene_wide[, diff_mean_LD_density := mean_LD_density_Marine - mean_LD_density_Freshwater]
gene_wide[, diff_mean_high_pairs := mean_high_pairs_Marine - mean_high_pairs_Freshwater]
gene_wide[, diff_mean_r2 := mean_mean_r2_Marine - mean_mean_r2_Freshwater]
gene_wide[, diff_max_r2 := mean_max_r2_Marine - mean_max_r2_Freshwater]

fwrite(gene_wide, file.path(outdir, "gene_level_marine_vs_freshwater_LD.tsv"), sep = "\t")

# =========================
# Statistical tests
# Paired across genes
# =========================

cat("\n================ Paired Wilcoxon tests across genes ================\n")

test_density <- wilcox.test(
  gene_wide$mean_LD_density_Marine,
  gene_wide$mean_LD_density_Freshwater,
  paired = TRUE,
  exact = FALSE
)

test_highpairs <- wilcox.test(
  gene_wide$mean_high_pairs_Marine,
  gene_wide$mean_high_pairs_Freshwater,
  paired = TRUE,
  exact = FALSE
)

test_meanr2 <- wilcox.test(
  gene_wide$mean_mean_r2_Marine,
  gene_wide$mean_mean_r2_Freshwater,
  paired = TRUE,
  exact = FALSE
)

test_maxr2 <- wilcox.test(
  gene_wide$mean_max_r2_Marine,
  gene_wide$mean_max_r2_Freshwater,
  paired = TRUE,
  exact = FALSE
)

print(test_density)
print(test_highpairs)
print(test_meanr2)
print(test_maxr2)

test_out <- data.table(
  metric = c("mean_LD_density", "mean_high_pairs", "mean_mean_r2", "mean_max_r2"),
  p_value = c(
    test_density$p.value,
    test_highpairs$p.value,
    test_meanr2$p.value,
    test_maxr2$p.value
  ),
  marine_mean = c(
    mean(gene_wide$mean_LD_density_Marine),
    mean(gene_wide$mean_high_pairs_Marine),
    mean(gene_wide$mean_mean_r2_Marine),
    mean(gene_wide$mean_max_r2_Marine)
  ),
  freshwater_mean = c(
    mean(gene_wide$mean_LD_density_Freshwater),
    mean(gene_wide$mean_high_pairs_Freshwater),
    mean(gene_wide$mean_mean_r2_Freshwater),
    mean(gene_wide$mean_max_r2_Freshwater)
  ),
  marine_median = c(
    median(gene_wide$mean_LD_density_Marine),
    median(gene_wide$mean_high_pairs_Marine),
    median(gene_wide$mean_mean_r2_Marine),
    median(gene_wide$mean_max_r2_Marine)
  ),
  freshwater_median = c(
    median(gene_wide$mean_LD_density_Freshwater),
    median(gene_wide$mean_high_pairs_Freshwater),
    median(gene_wide$mean_mean_r2_Freshwater),
    median(gene_wide$mean_max_r2_Freshwater)
  )
)

fwrite(test_out, file.path(outdir, "marine_vs_freshwater_LD_tests.tsv"), sep = "\t")

print(test_out)

# =========================
# Plot 1: LD density per gene
# =========================
p1 <- ggplot(gene_habitat, aes(x = habitat, y = mean_LD_density, fill = habitat)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.55) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Mean LD density per gene",
    title = "Marine vs freshwater LD density",
    subtitle = "LD density = proportion of SNP pairs with R² ≥ 0.2"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold")
  )

ggsave(
  file.path(outdir, "Marine_vs_Freshwater_LD_density_boxplot.png"),
  p1,
  width = 5.5,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, "Marine_vs_Freshwater_LD_density_boxplot.pdf"),
  p1,
  width = 5.5,
  height = 5
)

# =========================
# Plot 2: paired gene-level comparison
# =========================
gene_habitat[, gene_label := gene]

p2 <- ggplot(gene_habitat, aes(x = habitat, y = mean_LD_density, group = gene)) +
  geom_line(alpha = 0.25, linewidth = 0.35) +
  geom_point(aes(color = habitat), size = 1.8, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Mean LD density per gene",
    title = "Gene-level paired comparison of LD density",
    subtitle = "Each line connects the same gene between marine and freshwater groups"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold")
  )

ggsave(
  file.path(outdir, "Marine_vs_Freshwater_LD_density_paired_genes.png"),
  p2,
  width = 5.5,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, "Marine_vs_Freshwater_LD_density_paired_genes.pdf"),
  p2,
  width = 5.5,
  height = 5
)

# =========================
# Plot 3: population-level LD density
# =========================
p3 <- ggplot(pop_sum, aes(x = reorder(pop, mean_LD_density), y = mean_LD_density, fill = habitat)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    x = "Population",
    y = "Mean LD density across genes",
    title = "Population-level LD density"
  ) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold")
  )

ggsave(
  file.path(outdir, "Population_level_LD_density_barplot.png"),
  p3,
  width = 6,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(outdir, "Population_level_LD_density_barplot.pdf"),
  p3,
  width = 6,
  height = 7
)

# =========================
# C3-only comparison:
# Marine RS/SAY vs C3 freshwater
# This is useful because RS/SAY belong to mtCluster C3.
# =========================
c3_fw_pops <- c("ROB", "WB", "LG", "SL", "LAW", "BOOT", "JOE", "FG")
c3_pops <- c(marine_pops, c3_fw_pops)

c3_dat <- dat_full[pop %in% c3_pops]
c3_dat[, group := ifelse(pop %in% marine_pops, "Marine_C3", "Freshwater_C3")]
c3_dat[, group := factor(group, levels = c("Marine_C3", "Freshwater_C3"))]

c3_gene <- c3_dat[, .(
  mean_LD_density = mean(LD_density),
  mean_high_pairs = mean(high_pairs),
  mean_mean_r2 = mean(mean_r2),
  mean_max_r2 = mean(max_r2)
), by = .(gene, group)]

c3_wide <- dcast(
  c3_gene,
  gene ~ group,
  value.var = c("mean_LD_density", "mean_high_pairs", "mean_mean_r2", "mean_max_r2")
)

c3_test_density <- wilcox.test(
  c3_wide$mean_LD_density_Marine_C3,
  c3_wide$mean_LD_density_Freshwater_C3,
  paired = TRUE,
  exact = FALSE
)

cat("\n================ C3-only Marine vs C3 freshwater ================\n")
print(c3_test_density)

fwrite(c3_wide, file.path(outdir, "gene_level_marine_vs_C3freshwater_LD.tsv"), sep = "\t")

p4 <- ggplot(c3_gene, aes(x = group, y = mean_LD_density, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.55) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Mean LD density per gene",
    title = "Marine C3 vs freshwater C3 LD density",
    subtitle = "RS/SAY compared only against C3 freshwater populations"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold")
  )

ggsave(
  file.path(outdir, "MarineC3_vs_FreshwaterC3_LD_density_boxplot.png"),
  p4,
  width = 5.8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, "MarineC3_vs_FreshwaterC3_LD_density_boxplot.pdf"),
  p4,
  width = 5.8,
  height = 5
)

cat("\nDONE\n")