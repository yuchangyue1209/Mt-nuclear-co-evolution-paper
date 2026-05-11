fst based
#/work/cyu/TableS2.csv all genes bed

python3 - <<'PY'
import csv

infile = "/work/cyu/TableS2.csv"
outfile = "/work/cyu/gene_regions.bed"

with open(infile, newline="") as f, open(outfile, "w") as out:
    reader = csv.DictReader(f)
    for row in reader:
        chrom = row["chr_filled"].strip()
        start = row["region_start"].strip()
        end   = row["region_end"].strip()
        gene  = row["stickleback_name"].strip()

        if not chrom or not start or not end or not gene:
            continue

        s = int(start) - 1   # BED start is 0-based
        e = int(end)

        out.write(f"{chrom}\t{s}\t{e}\t{gene}\n")
PY

#step2
mkdir -p /mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work
cd /mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work

awk '{
  print $1 "\t" $2-1 "\t" $2 "\t" $0
}' OFS='\t' /mnt/spareHD_2/nuclear_marked_duplicates/nuclear.sync > nuclear.sync.bed

#
bedtools intersect \
  -a nuclear.sync.bed \
  -b /work/cyu/gene_regions.bed \
  -wa -wb > nuclear_gene.bed.sync

#
mkdir -p gene_sync_raw

awk '{
  gene=$NF
  print $0 >> "gene_sync_raw/" gene ".bed.sync"
}' nuclear_gene.bed.sync

#
mkdir -p gene_sync_fixed

for f in gene_sync_raw/*.bed.sync; do
  base=$(basename "$f" .bed.sync)
  awk '{
    print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33
  }' OFS='\t' "$f" > "gene_sync_fixed/${base}.sync"
done

mkdir -p gene_sync_dedup

for f in gene_sync_fixed/*.sync; do
  base=$(basename "$f" .sync)
  awk '!seen[$1 FS $2]++' "$f" > "gene_sync_dedup/${base}.sync"
done


conda activate /home/cyu/y/envs/grenedalf

mkdir -p /mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/fst_out_all
cd /mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/fst_out_all

for f in /mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/gene_sync_dedup_clean/*.sync; do
  gene=$(basename "$f" .sync)

  awk -F',' -v g="$gene" '{print g "." NR "\t" $1}' \
    /work/cyu/poolseq/PPalign_output/fst/poolsize.txt > "${gene}.rename.txt"

  grenedalf fst \
    --method unbiased-nei \
    --sync-path "$f" \
    --rename-samples-list "${gene}.rename.txt" \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting

  mv fst.csv "${gene}_fst.csv"
done

ls /mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/fst_out_all/*_fst.csv | wc -l



#mt genome fst without dloop
/work/cyu/poolseq/PPalign_output/ann
bedtools intersect \
  -a fish.sync.bed \
  -b dloop.bed \
  -v > fish_noDloop.bed.sync
  cut -f4-33 fish_noDloop.bed.sync > fish_noDloop.sync
head fish_noDloop.sync
awk 'NR<=5{print NF}' fish_noDloop.sync
cd /work/cyu/poolseq/PPalign_output/fst

awk -F',' '{print "fish_noDloop." NR "," $1}' poolsize.txt > rename_mt_noDloop.txt
conda activate /home/cyu/y/envs/grenedalf
mkdir -p /work/cyu/poolseq/PPalign_output/fst_mt_noDloop
cd /work/cyu/poolseq/PPalign_output/fst_mt_noDloop

grenedalf fst \
  --method unbiased-nei \
  --sync-path /work/cyu/poolseq/PPalign_output/ann/fish_noDloop.sync \
  --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename_mt_noDloop.txt \
  --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
  --window-type chromosomes \
  --window-average-policy valid-loci \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --no-extra-columns \
  --allow-file-overwriting

mv fst.csv mtgenome_noDloop_fst.csv

/work/cyu/poolseq/PPalign_output/fst_mt_noDloop/mtgenome_noDloop_fst.csv


#Mantel test：gene vs mt genome 是不是有gene和 mt的分化一致

/work/cyu/mantel_mt_noDloop_vs_nuclear_genes.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
  library(ggplot2)
})

# ============================================================
# Paths
# ============================================================
MT_FST_FILE <- "/work/cyu/poolseq/PPalign_output/fst_mt_noDloop/mtgenome_noDloop_fst.csv"
NUCLEAR_FST_DIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/fst_out_all"
OUTDIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/mantel_mt_vs_nuclear"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper: grenedalf fst csv -> symmetric matrix
# Input format:
# chrom,start,end,THE:JOE.fst,THE:BEA.fst,...
# chrX,1,15049116,0.40,0.30,...
# ============================================================
fst_to_matrix <- function(df) {
  if (nrow(df) < 1) stop("FST file has no rows.")
  if (ncol(df) < 4) stop("FST file has too few columns.")

  row <- df[1]

  pair_cols <- colnames(df)[4:ncol(df)]
  values <- as.numeric(unlist(row[, 4:ncol(df), with = FALSE]))

  # remove trailing ".fst"
  pair_names <- sub("\\.fst$", "", pair_cols)

  # split "THE:JOE"
  sp <- tstrsplit(pair_names, ":", fixed = TRUE)
  pop1 <- sp[[1]]
  pop2 <- sp[[2]]

  pops <- sort(unique(c(pop1, pop2)))

  mat <- matrix(0, nrow = length(pops), ncol = length(pops),
                dimnames = list(pops, pops))

  for (i in seq_along(values)) {
    if (is.na(pop1[i]) || is.na(pop2[i])) next
    mat[pop1[i], pop2[i]] <- values[i]
    mat[pop2[i], pop1[i]] <- values[i]
  }

  diag(mat) <- 0
  mat
}

# ============================================================
# Helper: run Mantel safely
# ============================================================
safe_mantel <- function(mt_mat, g_mat, gene, perms = 999) {
  common <- intersect(rownames(mt_mat), rownames(g_mat))
  common <- sort(common)

  if (length(common) < 4) {
    return(data.table(
      gene = gene,
      n_pops = length(common),
      n_pairs = NA_integer_,
      mantel_r = NA_real_,
      p = NA_real_,
      status = "too_few_pops"
    ))
  }

  mt_sub <- mt_mat[common, common, drop = FALSE]
  g_sub  <- g_mat[common, common, drop = FALSE]

  # vectorized distances
  x <- as.vector(as.dist(mt_sub))
  y <- as.vector(as.dist(g_sub))

  ok <- is.finite(x) & is.finite(y)
  x_ok <- x[ok]
  y_ok <- y[ok]

  if (length(x_ok) < 10) {
    return(data.table(
      gene = gene,
      n_pops = length(common),
      n_pairs = length(x_ok),
      mantel_r = NA_real_,
      p = NA_real_,
      status = "too_few_pairs"
    ))
  }

  if (var(x_ok) == 0 || var(y_ok) == 0) {
    return(data.table(
      gene = gene,
      n_pops = length(common),
      n_pairs = length(x_ok),
      mantel_r = NA_real_,
      p = NA_real_,
      status = "zero_variance"
    ))
  }

  # rebuild cleaned matrices if any NA existed
  # easiest robust route: use correlation on cleaned vectors + permutation manually
  obs_r <- suppressWarnings(cor(x_ok, y_ok, method = "pearson"))

  if (!is.finite(obs_r)) {
    return(data.table(
      gene = gene,
      n_pops = length(common),
      n_pairs = length(x_ok),
      mantel_r = NA_real_,
      p = NA_real_,
      status = "cor_failed"
    ))
  }

  # permutation test on labels of one matrix
  n <- nrow(mt_sub)
  perm_r <- numeric(perms)

  for (b in seq_len(perms)) {
    idx <- sample.int(n)
    g_perm <- g_sub[idx, idx, drop = FALSE]

    xp <- as.vector(as.dist(mt_sub))
    yp <- as.vector(as.dist(g_perm))

    okp <- is.finite(xp) & is.finite(yp)
    xp <- xp[okp]
    yp <- yp[okp]

    if (length(xp) < 10 || var(xp) == 0 || var(yp) == 0) {
      perm_r[b] <- NA_real_
    } else {
      perm_r[b] <- suppressWarnings(cor(xp, yp, method = "pearson"))
    }
  }

  perm_r <- perm_r[is.finite(perm_r)]

  if (length(perm_r) == 0) {
    return(data.table(
      gene = gene,
      n_pops = length(common),
      n_pairs = length(x_ok),
      mantel_r = obs_r,
      p = NA_real_,
      status = "perm_failed"
    ))
  }

  pval <- (sum(abs(perm_r) >= abs(obs_r)) + 1) / (length(perm_r) + 1)

  data.table(
    gene = gene,
    n_pops = length(common),
    n_pairs = length(x_ok),
    mantel_r = obs_r,
    p = pval,
    status = "ok"
  )
}

# ============================================================
# Read mt matrix
# ============================================================
cat("[read] mt no-Dloop FST:\n", MT_FST_FILE, "\n", sep = "")
mt_df <- fread(MT_FST_FILE)
mt_mat <- fst_to_matrix(mt_df)

cat("[info] mt matrix dimensions: ", nrow(mt_mat), " x ", ncol(mt_mat), "\n", sep = "")
cat("[info] mt populations:\n")
print(rownames(mt_mat))

# ============================================================
# Read nuclear files
# ============================================================
files <- list.files(
  NUCLEAR_FST_DIR,
  pattern = "_fst\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No nuclear *_fst.csv files found in: ", NUCLEAR_FST_DIR)
}

cat("[info] nuclear FST files found: ", length(files), "\n", sep = "")

# ============================================================
# Run Mantel for each gene
# ============================================================
results <- vector("list", length(files))

for (i in seq_along(files)) {
  f <- files[i]
  gene <- sub("_fst\\.csv$", "", basename(f))

  cat("[run] ", i, "/", length(files), " : ", gene, "\n", sep = "")

  res_i <- tryCatch({
    df <- fread(f)
    g_mat <- fst_to_matrix(df)
    safe_mantel(mt_mat, g_mat, gene, perms = 999)
  }, error = function(e) {
    data.table(
      gene = gene,
      n_pops = NA_integer_,
      n_pairs = NA_integer_,
      mantel_r = NA_real_,
      p = NA_real_,
      status = paste0("error: ", conditionMessage(e))
    )
  })

  results[[i]] <- res_i
}

res <- rbindlist(results, use.names = TRUE, fill = TRUE)

# ============================================================
# Multiple testing correction
# ============================================================
res[, p_bh := NA_real_]
okp <- is.finite(res$p) & !is.na(res$p)
res[okp, p_bh := p.adjust(p, method = "BH")]

setorder(res, p, -mantel_r)

# ============================================================
# Write results
# ============================================================
OUT_TSV <- file.path(OUTDIR, "mantel_mt_noDloop_vs_nuclear_genes.tsv")
OUT_TOP <- file.path(OUTDIR, "mantel_top20.tsv")
OUT_FAIL <- file.path(OUTDIR, "mantel_failed_or_skipped.tsv")

fwrite(res, OUT_TSV, sep = "\t")
fwrite(res[1:min(20, .N)], OUT_TOP, sep = "\t")
fwrite(res[status != "ok"], OUT_FAIL, sep = "\t")

cat("\n[write] all results: ", OUT_TSV, "\n", sep = "")
cat("[write] top20: ", OUT_TOP, "\n", sep = "")
cat("[write] failed/skipped: ", OUT_FAIL, "\n", sep = "")

cat("\n=== top 20 genes ===\n")
print(res[1:min(20, .N)])

cat("\n=== status summary ===\n")
print(res[, .N, by = status][order(-N)])

# ============================================================
# Plots
# ============================================================
plot_dt <- res[is.finite(mantel_r) & is.finite(p)]

if (nrow(plot_dt) > 0) {
  p1 <- ggplot(plot_dt, aes(x = mantel_r, y = -log10(p))) +
    geom_point() +
    theme_classic(base_size = 14) +
    labs(
      x = "Mantel r",
      y = expression(-log[10](p)),
      title = "Mantel test: mt genome (no D-loop) vs nuclear gene FST"
    )

  ggsave(
    filename = file.path(OUTDIR, "Fig_mantel_volcano.png"),
    plot = p1, width = 7, height = 5, dpi = 300
  )

  top20 <- res[is.finite(mantel_r)][order(p)][1:min(20, sum(is.finite(res$mantel_r)))]
  p2 <- ggplot(top20, aes(x = reorder(gene, mantel_r), y = mantel_r)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      x = "Gene",
      y = "Mantel r",
      title = "Top 20 nuclear genes matching mt divergence"
    )

  ggsave(
    filename = file.path(OUTDIR, "Fig_mantel_top20_barplot.png"),
    plot = p2, width = 7, height = 8, dpi = 300
  )

  cat("[write] plots saved in: ", OUTDIR, "\n", sep = "")
}

cat("\n[OK] done.\n")

#cluster
library(data.table)

# ============================================================
# 0. Input / output paths
# ============================================================
FST_DIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/fst_out_all"
OUTDIR  <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/mtcluster_fst"
PAIR_OUTFILE <- file.path(OUTDIR, "all_pairs_long.tsv")
RES_OUTFILE  <- file.path(OUTDIR, "gene_fst_within_between_mtcluster.tsv")

MT_CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Read all gene-level FST files and convert to long format
# ============================================================
files <- list.files(FST_DIR, pattern = "_fst\\.csv$", full.names = TRUE)

cat("Number of gene FST files:", length(files), "\n")

if (length(files) == 0) {
  stop("No *_fst.csv files found in: ", FST_DIR)
}

pair_list <- vector("list", length(files))

for (i in seq_along(files)) {
  f <- files[i]
  gene <- sub("_fst\\.csv$", "", basename(f))

  df <- fread(f)

  if (nrow(df) < 1 || ncol(df) < 4) {
    next
  }

  row <- df[1]

  # pairwise FST columns start at column 4
  pairs <- colnames(df)[4:ncol(df)]
  vals  <- as.numeric(unlist(row[, 4:ncol(df), with = FALSE]))

  pair_names <- sub("\\.fst$", "", pairs)
  sp <- tstrsplit(pair_names, ":", fixed = TRUE)

  pair_list[[i]] <- data.table(
    gene = gene,
    pop1 = sp[[1]],
    pop2 = sp[[2]],
    fst  = vals
  )
}

dt_pair <- rbindlist(pair_list, use.names = TRUE, fill = TRUE)

# remove missing fst rows if any
dt_pair <- dt_pair[is.finite(fst)]

fwrite(dt_pair, PAIR_OUTFILE, sep = "\t")

cat("DONE. Pairwise long table written to:\n", PAIR_OUTFILE, "\n")
cat("Preview of dt_pair:\n")
print(head(dt_pair))

# ============================================================
# 2. Read and clean mt cluster file
# ============================================================
mt_cluster <- fread(
  MT_CLUSTER_FILE,
  fill = TRUE,
  header = FALSE
)

# remove empty lines
mt_cluster <- mt_cluster[V1 != "" & !is.na(V1)]

# keep first two columns only
mt_cluster <- mt_cluster[, .(pop = V1, cluster = V2)]

# remove accidental header row if present
mt_cluster <- mt_cluster[pop != "pop"]

cat("\nmt_cluster preview:\n")
print(head(mt_cluster))
cat("Unique clusters:\n")
print(unique(mt_cluster$cluster))
cat("Number of populations in mt_cluster:", length(unique(mt_cluster$pop)), "\n")

# ============================================================
# 3. For each gene: compare within-cluster vs between-cluster FST
# ============================================================
gene_list <- unique(dt_pair$gene)

cat("\nNumber of genes in dt_pair:", length(gene_list), "\n")

res_list <- vector("list", length(gene_list))

for (i in seq_along(gene_list)) {
  gene_i <- gene_list[i]

  dt_g <- dt_pair[gene == gene_i]

  # merge cluster for pop1
  dt_g <- merge(
    dt_g,
    mt_cluster,
    by.x = "pop1",
    by.y = "pop",
    all.x = TRUE
  )
  setnames(dt_g, "cluster", "cluster1")

  # merge cluster for pop2
  dt_g <- merge(
    dt_g,
    mt_cluster,
    by.x = "pop2",
    by.y = "pop",
    all.x = TRUE
  )
  setnames(dt_g, "cluster", "cluster2")

  # keep only rows where both populations have mt cluster
  dt_g <- dt_g[!is.na(cluster1) & !is.na(cluster2)]

  # define within vs between
  dt_g[, within := cluster1 == cluster2]

  n_within  <- sum(dt_g$within, na.rm = TRUE)
  n_between <- sum(!dt_g$within, na.rm = TRUE)

  # skip genes without enough comparisons
  if (n_within < 1 || n_between < 1) {
    res_list[[i]] <- data.table(
      gene = gene_i,
      fst_within = NA_real_,
      fst_between = NA_real_,
      diff = NA_real_,
      p = NA_real_,
      n_within = n_within,
      n_between = n_between
    )
    next
  }

  fst_within  <- mean(dt_g$fst[dt_g$within == TRUE], na.rm = TRUE)
  fst_between <- mean(dt_g$fst[dt_g$within == FALSE], na.rm = TRUE)
  diff_i <- fst_between - fst_within

  p_i <- tryCatch(
    wilcox.test(fst ~ within, data = dt_g)$p.value,
    error = function(e) NA_real_
  )

  res_list[[i]] <- data.table(
    gene = gene_i,
    fst_within = fst_within,
    fst_between = fst_between,
    diff = diff_i,
    p = p_i,
    n_within = n_within,
    n_between = n_between
  )
}

res_fst <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

# multiple testing correction
res_fst[, p_bh := p.adjust(p, method = "BH")]

# sort by strongest cluster structuring
setorder(res_fst, -diff)

# write result
fwrite(res_fst, RES_OUTFILE, sep = "\t")

cat("\nDONE. Gene-level within vs between mt-cluster FST results written to:\n",
    RES_OUTFILE, "\n", sep = "")

cat("\nTop 20 genes:\n")
print(res_fst[1:min(20, .N)])

cat("\nSummary:\n")
print(summary(res_fst$diff))







#copen
cat /work/cyu/cophenetic_mt_vs_nuclear_genes.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
})

# ============================================================
# Paths
# ============================================================
MT_FST_FILE <- "/work/cyu/poolseq/PPalign_output/fst_mt_noDloop/mtgenome_noDloop_fst.csv"
NUCLEAR_FST_DIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/fst_out_all"
OUTDIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/cophenetic_mt_vs_nuclear"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper: grenedalf fst csv -> symmetric matrix
# ============================================================
fst_to_matrix <- function(df) {
  if (nrow(df) < 1) stop("FST file has no rows.")
  if (ncol(df) < 4) stop("FST file has too few columns.")

  row <- df[1]
  pair_cols <- colnames(df)[4:ncol(df)]
  values <- as.numeric(unlist(row[, 4:ncol(df), with = FALSE]))

  pair_names <- sub("\\.fst$", "", pair_cols)
  sp <- tstrsplit(pair_names, ":", fixed = TRUE)
  pop1 <- sp[[1]]
  pop2 <- sp[[2]]

  pops <- sort(unique(c(pop1, pop2)))
  mat <- matrix(0, nrow = length(pops), ncol = length(pops),
                dimnames = list(pops, pops))

  for (i in seq_along(values)) {
    if (is.na(pop1[i]) || is.na(pop2[i])) next
    mat[pop1[i], pop2[i]] <- values[i]
    mat[pop2[i], pop1[i]] <- values[i]
  }

  diag(mat) <- 0
  mat
}

# ============================================================
# Helper: safe NJ + cophenetic correlation
# ============================================================
safe_cophenetic_cor <- function(mt_mat, g_mat, gene, perms = 999) {
  common <- intersect(rownames(mt_mat), rownames(g_mat))
  common <- sort(common)

  if (length(common) < 4) {
    return(data.table(
      gene = gene,
      n_pops = length(common),
      n_pairs = NA_integer_,
      cophenetic_r = NA_real_,
      p = NA_real_,
      status = "too_few_pops"
    ))
  }

  mt_sub <- mt_mat[common, common, drop = FALSE]
  g_sub  <- g_mat[common, common, drop = FALSE]

  # remove NA/Inf from original matrices if any
  if (any(!is.finite(mt_sub)) || any(!is.finite(g_sub))) {
    return(data.table(
      gene = gene,
      n_pops = length(common),
      n_pairs = NA_integer_,
      cophenetic_r = NA_real_,
      p = NA_real_,
      status = "non_finite_matrix"
    ))
  }

  # NJ trees
  mt_tree <- tryCatch(nj(as.dist(mt_sub)), error = function(e) NULL)
  g_tree  <- tryCatch(nj(as.dist(g_sub)),  error = function(e) NULL)

  if (is.null(mt_tree) || is.null(g_tree)) {
    return(data.table(
      gene = gene,
      n_pops = length(common),
      n_pairs = NA_integer_,
      cophenetic_r = NA_real_,
      p = NA_real_,
      status = "tree_failed"
    ))
  }

  mt_coph <- cophenetic(mt_tree)
  g_coph  <- cophenetic(g_tree)

  common2 <- intersect(rownames(mt_coph), rownames(g_coph))
  common2 <- sort(common2)

  mt_coph <- mt_coph[common2, common2, drop = FALSE]
  g_coph  <- g_coph[common2, common2, drop = FALSE]

  x <- as.vector(as.dist(mt_coph))
  y <- as.vector(as.dist(g_coph))

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 10) {
    return(data.table(
      gene = gene,
      n_pops = length(common2),
      n_pairs = length(x),
      cophenetic_r = NA_real_,
      p = NA_real_,
      status = "too_few_pairs"
    ))
  }

  if (var(x) == 0 || var(y) == 0) {
    return(data.table(
      gene = gene,
      n_pops = length(common2),
      n_pairs = length(x),
      cophenetic_r = NA_real_,
      p = NA_real_,
      status = "zero_variance"
    ))
  }

  obs_r <- suppressWarnings(cor(x, y, method = "pearson"))

  if (!is.finite(obs_r)) {
    return(data.table(
      gene = gene,
      n_pops = length(common2),
      n_pairs = length(x),
      cophenetic_r = NA_real_,
      p = NA_real_,
      status = "cor_failed"
    ))
  }

  # permutation by relabeling gene tree tips
  n <- nrow(g_coph)
  perm_r <- numeric(perms)

  for (b in seq_len(perms)) {
    idx <- sample.int(n)
    g_perm <- g_coph[idx, idx, drop = FALSE]

    xp <- as.vector(as.dist(mt_coph))
    yp <- as.vector(as.dist(g_perm))

    okp <- is.finite(xp) & is.finite(yp)
    xp <- xp[okp]
    yp <- yp[okp]

    if (length(xp) < 10 || var(xp) == 0 || var(yp) == 0) {
      perm_r[b] <- NA_real_
    } else {
      perm_r[b] <- suppressWarnings(cor(xp, yp, method = "pearson"))
    }
  }

  perm_r <- perm_r[is.finite(perm_r)]

  if (length(perm_r) == 0) {
    return(data.table(
      gene = gene,
      n_pops = length(common2),
      n_pairs = length(x),
      cophenetic_r = obs_r,
      p = NA_real_,
      status = "perm_failed"
    ))
  }

  pval <- (sum(abs(perm_r) >= abs(obs_r)) + 1) / (length(perm_r) + 1)

  data.table(
    gene = gene,
    n_pops = length(common2),
    n_pairs = length(x),
    cophenetic_r = obs_r,
    p = pval,
    status = "ok"
  )
}

# ============================================================
# Read mt matrix
# ============================================================
cat("[read] mt no-Dloop FST:\n", MT_FST_FILE, "\n", sep = "")
mt_df <- fread(MT_FST_FILE)
mt_mat <- fst_to_matrix(mt_df)

cat("[info] mt matrix dimensions: ", nrow(mt_mat), " x ", ncol(mt_mat), "\n", sep = "")

# ============================================================
# Read all nuclear gene FST files
# ============================================================
files <- list.files(
  NUCLEAR_FST_DIR,
  pattern = "_fst\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No nuclear *_fst.csv files found in: ", NUCLEAR_FST_DIR)
}

cat("[info] nuclear FST files found: ", length(files), "\n", sep = "")

# ============================================================
# Run cophenetic correlation per gene
# ============================================================
results <- vector("list", length(files))

for (i in seq_along(files)) {
  f <- files[i]
  gene <- sub("_fst\\.csv$", "", basename(f))

  cat("[run] ", i, "/", length(files), " : ", gene, "\n", sep = "")

  res_i <- tryCatch({
    df <- fread(f)
    g_mat <- fst_to_matrix(df)
    safe_cophenetic_cor(mt_mat, g_mat, gene, perms = 999)
  }, error = function(e) {
    data.table(
      gene = gene,
      n_pops = NA_integer_,
      n_pairs = NA_integer_,
      cophenetic_r = NA_real_,
      p = NA_real_,
      status = paste0("error: ", conditionMessage(e))
    )
  })

  results[[i]] <- res_i
}

res_tree <- rbindlist(results, use.names = TRUE, fill = TRUE)

# multiple testing
res_tree[, p_bh := NA_real_]
okp <- is.finite(res_tree$p) & !is.na(res_tree$p)
res_tree[okp, p_bh := p.adjust(p, method = "BH")]

setorder(res_tree, p, -cophenetic_r)

# ============================================================
# Write results
# ============================================================
OUT_TSV  <- file.path(OUTDIR, "cophenetic_mt_noDloop_vs_nuclear_genes.tsv")
OUT_TOP  <- file.path(OUTDIR, "cophenetic_top20.tsv")
OUT_FAIL <- file.path(OUTDIR, "cophenetic_failed_or_skipped.tsv")

fwrite(res_tree, OUT_TSV, sep = "\t")
fwrite(res_tree[1:min(20, .N)], OUT_TOP, sep = "\t")
fwrite(res_tree[status != "ok"], OUT_FAIL, sep = "\t")

cat("\n[write] all results: ", OUT_TSV, "\n", sep = "")
cat("[write] top20: ", OUT_TOP, "\n", sep = "")
cat("[write] failed/skipped: ", OUT_FAIL, "\n", sep = "")

cat("\n=== top 20 genes ===\n")
print(res_tree[1:min(20, .N)])

cat("\n=== status summary ===\n")
print(res_tree[, .N, by = status][order(-N)])

# ============================================================
# Plots
# ============================================================
plot_dt <- res_tree[is.finite(cophenetic_r) & is.finite(p)]

if (nrow(plot_dt) > 0) {
  p1 <- ggplot(plot_dt, aes(x = cophenetic_r, y = -log10(p))) +
    geom_point() +
    theme_classic(base_size = 14) +
    labs(
      x = "Cophenetic correlation (r)",
      y = expression(-log[10](p)),
      title = "Cophenetic correlation: mt tree vs nuclear gene tree"
    )

  ggsave(
    filename = file.path(OUTDIR, "Fig_cophenetic_volcano.png"),
    plot = p1, width = 7, height = 5, dpi = 300
  )

  top20 <- res_tree[is.finite(cophenetic_r)][order(p)][1:min(20, sum(is.finite(res_tree$cophenetic_r)))]
  p2 <- ggplot(top20, aes(x = reorder(gene, cophenetic_r), y = cophenetic_r)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      x = "Gene",
      y = "Cophenetic r",
      title = "Top 20 nuclear genes whose trees most resemble mt tree"
    )

  ggsave(
    filename = file.path(OUTDIR, "Fig_cophenetic_top20_barplot.png"),
    plot = p2, width = 7, height = 8, dpi = 300
  )

  cat("[write] plots saved in: ", OUTDIR, "\n", sep = "")
}

cat("\n[OK] done.\n")







####pbs gene level
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# paths
# =========================
MT_FST_FILE <- "/work/cyu/poolseq/PPalign_output/fst_mt_noDloop/mtgenome_noDloop_fst.csv"
NU_FST_DIR  <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/fst_out_all"

OUTDIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/pbs_mt_vs_nuclear"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# =========================
# population design
# =========================
AK_FW <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_FW <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")

# marine references
# AK freshwater: marine = RS, outgroup = SAY
# BC freshwater: marine = SAY, outgroup = RS
design <- rbindlist(list(
  data.table(focal = AK_FW, region = "AK", marine = "RS",  outgroup = "SAY"),
  data.table(focal = BC_FW, region = "BC", marine = "SAY", outgroup = "RS")
))

# remove populations not present if needed later
design <- unique(design)

# =========================
# helper functions
# =========================
fst_to_matrix <- function(file) {
  df <- fread(file)

  if (nrow(df) < 1 || ncol(df) < 4) {
    stop("Bad FST file: ", file)
  }

  row <- df[1]
  pair_cols <- colnames(df)[4:ncol(df)]
  vals <- as.numeric(unlist(row[, 4:ncol(df), with = FALSE]))

  pair_names <- sub("\\.fst$", "", pair_cols)
  sp <- tstrsplit(pair_names, ":", fixed = TRUE)

  pop1 <- sp[[1]]
  pop2 <- sp[[2]]

  pops <- sort(unique(c(pop1, pop2)))

  mat <- matrix(NA_real_, nrow = length(pops), ncol = length(pops),
                dimnames = list(pops, pops))

  diag(mat) <- 0

  for (i in seq_along(vals)) {
    mat[pop1[i], pop2[i]] <- vals[i]
    mat[pop2[i], pop1[i]] <- vals[i]
  }

  mat
}

get_fst <- function(mat, a, b) {
  if (!(a %in% rownames(mat)) || !(b %in% colnames(mat))) return(NA_real_)
  mat[a, b]
}

calc_pbs_one <- function(mat, focal, marine, outgroup) {
  fst_AB <- get_fst(mat, focal, marine)
  fst_AC <- get_fst(mat, focal, outgroup)
  fst_BC <- get_fst(mat, marine, outgroup)

  # valid range protection
  fst_AB <- pmin(pmax(fst_AB, 0), 0.999999)
  fst_AC <- pmin(pmax(fst_AC, 0), 0.999999)
  fst_BC <- pmin(pmax(fst_BC, 0), 0.999999)

  if (!all(is.finite(c(fst_AB, fst_AC, fst_BC)))) {
    return(data.table(
      fst_focal_marine = fst_AB,
      fst_focal_outgroup = fst_AC,
      fst_marine_outgroup = fst_BC,
      T_focal_marine = NA_real_,
      T_focal_outgroup = NA_real_,
      T_marine_outgroup = NA_real_,
      PBS = NA_real_
    ))
  }

  T_AB <- -log(1 - fst_AB)
  T_AC <- -log(1 - fst_AC)
  T_BC <- -log(1 - fst_BC)

  PBS <- (T_AB + T_AC - T_BC) / 2

  data.table(
    fst_focal_marine = fst_AB,
    fst_focal_outgroup = fst_AC,
    fst_marine_outgroup = fst_BC,
    T_focal_marine = T_AB,
    T_focal_outgroup = T_AC,
    T_marine_outgroup = T_BC,
    PBS = PBS
  )
}

calc_pbs_table <- function(mat, design_dt) {
  rbindlist(lapply(seq_len(nrow(design_dt)), function(i) {
    d <- design_dt[i]
    out <- calc_pbs_one(mat, d$focal, d$marine, d$outgroup)
    cbind(d, out)
  }))
}

# =========================
# 1. mt PBS
# =========================
cat("[read] mt FST\n")
mt_mat <- fst_to_matrix(MT_FST_FILE)

design_use <- design[
  focal %in% rownames(mt_mat) &
    marine %in% rownames(mt_mat) &
    outgroup %in% rownames(mt_mat)
]

mt_pbs <- calc_pbs_table(mt_mat, design_use)
setnames(mt_pbs, "PBS", "mt_PBS")

fwrite(mt_pbs, file.path(OUTDIR, "mt_noDloop_PBS_by_population.tsv"), sep = "\t")

cat("[write] mt PBS table\n")

# =========================
# 2. nuclear gene PBS
# =========================
files <- list.files(NU_FST_DIR, pattern = "_fst\\.csv$", full.names = TRUE)

cat("[info] nuclear gene FST files:", length(files), "\n")

nu_list <- vector("list", length(files))

for (i in seq_along(files)) {
  f <- files[i]
  gene <- sub("_fst\\.csv$", "", basename(f))

  cat("[run] ", i, "/", length(files), " ", gene, "\n", sep = "")

  nu_list[[i]] <- tryCatch({
    mat <- fst_to_matrix(f)
    tmp <- calc_pbs_table(mat, design_use)
    tmp[, gene := gene]
    tmp
  }, error = function(e) {
    data.table(
      focal = NA_character_,
      region = NA_character_,
      marine = NA_character_,
      outgroup = NA_character_,
      gene = gene,
      fst_focal_marine = NA_real_,
      fst_focal_outgroup = NA_real_,
      fst_marine_outgroup = NA_real_,
      T_focal_marine = NA_real_,
      T_focal_outgroup = NA_real_,
      T_marine_outgroup = NA_real_,
      PBS = NA_real_
    )
  })
}

nu_pbs <- rbindlist(nu_list, use.names = TRUE, fill = TRUE)
setnames(nu_pbs, "PBS", "nu_PBS")

fwrite(nu_pbs, file.path(OUTDIR, "nuclear_gene_PBS_by_population.tsv"), sep = "\t")

# =========================
# 3. merge mt and nuclear PBS
# =========================
merged <- merge(
  nu_pbs,
  mt_pbs[, .(focal, region, marine, outgroup, mt_PBS)],
  by = c("focal", "region", "marine", "outgroup"),
  all.x = TRUE
)

fwrite(merged, file.path(OUTDIR, "merged_mt_nuclear_PBS_by_gene_population.tsv"), sep = "\t")

# =========================
# 4. gene-wise correlation: nu PBS pattern vs mt PBS pattern
# =========================
gene_res <- merged[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))
      data.table(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    } else {
      data.table(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    }
  },
  by = gene
]

gene_res[, p_bh := p.adjust(p, method = "BH")]
setorder(gene_res, p, -cor_r)

fwrite(gene_res, file.path(OUTDIR, "gene_correlation_nuPBS_vs_mtPBS.tsv"), sep = "\t")

# =========================
# 5. optional: region-specific correlation
# =========================
gene_region_res <- merged[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 4 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))
      data.table(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    } else {
      data.table(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    }
  },
  by = .(gene, region)
]

gene_region_res[, p_bh := p.adjust(p, method = "BH"), by = region]
setorder(gene_region_res, region, p, -cor_r)

fwrite(gene_region_res, file.path(OUTDIR, "gene_region_correlation_nuPBS_vs_mtPBS.tsv"), sep = "\t")

# =========================
# 6. plots
# =========================
plot_dt <- gene_res[is.finite(cor_r) & is.finite(p)]

if (nrow(plot_dt) > 0) {
  p1 <- ggplot(plot_dt, aes(x = cor_r, y = -log10(p))) +
    geom_point() +
    theme_classic(base_size = 14) +
    labs(
      x = "Correlation between nuclear gene PBS and mt PBS",
      y = expression(-log[10](p)),
      title = "Genes whose PBS trajectories match mitochondrial PBS"
    )

  ggsave(
    file.path(OUTDIR, "Fig_gene_nuPBS_vs_mtPBS_volcano.png"),
    p1, width = 7, height = 5, dpi = 300
  )

  top20 <- gene_res[is.finite(cor_r)][order(p)][1:min(20, .N)]

  p2 <- ggplot(top20, aes(x = reorder(gene, cor_r), y = cor_r)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      x = "Gene",
      y = "Correlation r",
      title = "Top genes matching mt PBS trajectory"
    )

  ggsave(
    file.path(OUTDIR, "Fig_top20_gene_nuPBS_vs_mtPBS.png"),
    p2, width = 7, height = 8, dpi = 300
  )
}

cat("\n[OK] done\n")
cat("[outdir] ", OUTDIR, "\n", sep = "")







#pbsn1 gene level
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# paths
# =========================
MT_FST_FILE <- "/work/cyu/poolseq/PPalign_output/fst_mt_noDloop/mtgenome_noDloop_fst.csv"
NU_FST_DIR  <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/fst_out_all"

OUTDIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/pbsn1_gene_level"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# =========================
# population design
# =========================
AK_FW <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_FW <- c("SWA","THE","JOE","BEA","MUC","PYE","BOOT","ECHO","LAW","GOS","ROB")

design <- rbindlist(list(
  data.table(focal = AK_FW, region = "AK", marine = "RS",  outgroup = "SAY"),
  data.table(focal = BC_FW, region = "BC", marine = "SAY", outgroup = "RS")
))

# =========================
# helper functions
# =========================
fst_to_matrix <- function(file) {
  df <- fread(file)

  row <- df[1]
  pair_cols <- colnames(df)[4:ncol(df)]
  vals <- as.numeric(unlist(row[, 4:ncol(df), with = FALSE]))

  pair_names <- sub("\\.fst$", "", pair_cols)
  sp <- tstrsplit(pair_names, ":", fixed = TRUE)

  pop1 <- sp[[1]]
  pop2 <- sp[[2]]
  pops <- sort(unique(c(pop1, pop2)))

  mat <- matrix(NA_real_, nrow = length(pops), ncol = length(pops),
                dimnames = list(pops, pops))
  diag(mat) <- 0

  for (i in seq_along(vals)) {
    mat[pop1[i], pop2[i]] <- vals[i]
    mat[pop2[i], pop1[i]] <- vals[i]
  }

  mat
}

get_fst <- function(mat, a, b) {
  if (!(a %in% rownames(mat)) || !(b %in% colnames(mat))) return(NA_real_)
  mat[a, b]
}

calc_all_pbs <- function(mat, A, B, C) {
  fst_AB <- get_fst(mat, A, B)
  fst_AC <- get_fst(mat, A, C)
  fst_BC <- get_fst(mat, B, C)

  fst_AB <- pmin(pmax(fst_AB, 0), 0.999999)
  fst_AC <- pmin(pmax(fst_AC, 0), 0.999999)
  fst_BC <- pmin(pmax(fst_BC, 0), 0.999999)

  if (!all(is.finite(c(fst_AB, fst_AC, fst_BC)))) {
    return(data.table(
      fst_AB = fst_AB, fst_AC = fst_AC, fst_BC = fst_BC,
      PBS_A = NA_real_, PBS_B = NA_real_, PBS_C = NA_real_,
      PBSn1_A = NA_real_
    ))
  }

  T_AB <- -log(1 - fst_AB)
  T_AC <- -log(1 - fst_AC)
  T_BC <- -log(1 - fst_BC)

  PBS_A <- (T_AB + T_AC - T_BC) / 2
  PBS_B <- (T_AB + T_BC - T_AC) / 2
  PBS_C <- (T_AC + T_BC - T_AB) / 2

  denom <- PBS_A + PBS_B + PBS_C

  PBSn1_A <- ifelse(is.finite(denom) && denom > 0, PBS_A / denom, NA_real_)

  data.table(
    fst_AB = fst_AB,
    fst_AC = fst_AC,
    fst_BC = fst_BC,
    T_AB = T_AB,
    T_AC = T_AC,
    T_BC = T_BC,
    PBS_A = PBS_A,
    PBS_B = PBS_B,
    PBS_C = PBS_C,
    PBSn1_A = PBSn1_A
  )
}

calc_pbsn1_table <- function(mat, design_dt) {
  rbindlist(lapply(seq_len(nrow(design_dt)), function(i) {
    d <- design_dt[i]
    out <- calc_all_pbs(mat, d$focal, d$marine, d$outgroup)
    cbind(d, out)
  }))
}

# =========================
# 1. mt PBSn1
# =========================
cat("[read] mt FST\n")
mt_mat <- fst_to_matrix(MT_FST_FILE)

design_use <- design[
  focal %in% rownames(mt_mat) &
    marine %in% rownames(mt_mat) &
    outgroup %in% rownames(mt_mat)
]

mt_pbsn1 <- calc_pbsn1_table(mt_mat, design_use)

setnames(mt_pbsn1, "PBS_A", "mt_PBS")
setnames(mt_pbsn1, "PBSn1_A", "mt_PBSn1")

fwrite(mt_pbsn1, file.path(OUTDIR, "mt_noDloop_PBSn1_by_population.tsv"), sep = "\t")

# =========================
# 2. nuclear gene PBSn1
# =========================
files <- list.files(NU_FST_DIR, pattern = "_fst\\.csv$", full.names = TRUE)
cat("[info] nuclear gene FST files:", length(files), "\n")

nu_list <- vector("list", length(files))

for (i in seq_along(files)) {
  f <- files[i]
  gene <- sub("_fst\\.csv$", "", basename(f))

  cat("[run] ", i, "/", length(files), " ", gene, "\n", sep = "")

  nu_list[[i]] <- tryCatch({
    mat <- fst_to_matrix(f)
    tmp <- calc_pbsn1_table(mat, design_use)
    tmp[, gene := gene]
    tmp
  }, error = function(e) {
    data.table(gene = gene, error = conditionMessage(e))
  })
}

nu_pbsn1 <- rbindlist(nu_list, use.names = TRUE, fill = TRUE)

setnames(nu_pbsn1, "PBS_A", "nu_PBS")
setnames(nu_pbsn1, "PBSn1_A", "nu_PBSn1")

fwrite(nu_pbsn1, file.path(OUTDIR, "nuclear_gene_PBSn1_by_population.tsv"), sep = "\t")

# =========================
# 3. merge mt and nuclear PBSn1
# =========================
merged <- merge(
  nu_pbsn1,
  mt_pbsn1[, .(focal, region, marine, outgroup, mt_PBS, mt_PBSn1)],
  by = c("focal", "region", "marine", "outgroup"),
  all.x = TRUE
)

fwrite(merged, file.path(OUTDIR, "merged_mt_nuclear_PBSn1_by_gene_population.tsv"), sep = "\t")

# =========================
# 4. gene-wise correlation: nu PBSn1 vs mt PBSn1
# =========================
gene_res <- merged[
  is.finite(nu_PBSn1) & is.finite(mt_PBSn1),
  {
    if (.N >= 5 && var(nu_PBSn1) > 0 && var(mt_PBSn1) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBSn1, mt_PBSn1, method = "pearson"))
      data.table(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBSn1 = mean(nu_PBSn1, na.rm = TRUE),
        mean_mt_PBSn1 = mean(mt_PBSn1, na.rm = TRUE)
      )
    } else {
      data.table(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBSn1 = mean(nu_PBSn1, na.rm = TRUE),
        mean_mt_PBSn1 = mean(mt_PBSn1, na.rm = TRUE)
      )
    }
  },
  by = gene
]

gene_res[, p_bh := p.adjust(p, method = "BH")]
setorder(gene_res, p, -cor_r)

fwrite(gene_res, file.path(OUTDIR, "gene_correlation_nuPBSn1_vs_mtPBSn1.tsv"), sep = "\t")

# =========================
# 5. plot
# =========================
plot_dt <- gene_res[is.finite(cor_r) & is.finite(p)]

if (nrow(plot_dt) > 0) {
  p1 <- ggplot(plot_dt, aes(x = cor_r, y = -log10(p))) +
    geom_point() +
    theme_classic(base_size = 14) +
    labs(
      x = "Correlation between nuclear gene PBSn1 and mt PBSn1",
      y = expression(-log[10](p)),
      title = "Gene-level PBSn1 concordance with mitochondrial PBSn1"
    )

  ggsave(
    file.path(OUTDIR, "Fig_gene_nuPBSn1_vs_mtPBSn1_volcano.png"),
    p1, width = 7, height = 5, dpi = 300
  )

  top20 <- gene_res[is.finite(cor_r)][order(p)][1:min(20, .N)]

  p2 <- ggplot(top20, aes(x = reorder(gene, cor_r), y = cor_r)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      x = "Gene",
      y = "Correlation r",
      title = "Top genes matching mt PBSn1 trajectory"
    )

  ggsave(
    file.path(OUTDIR, "Fig_top20_gene_nuPBSn1_vs_mtPBSn1.png"),
    p2, width = 7, height = 8, dpi = 300
  )
}

cat("\n[OK] done\n")
cat("[outdir] ", OUTDIR, "\n", sep = "")




#max pbs gene level
/work/cyu/maxSNP_PBS_nuclear_genes.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# input paths
# =========================
SYNC_FILE <- "/mnt/spareHD_2/nuclear_marked_duplicates/nuclear.sync"
GENE_BED  <- "/work/cyu/nuOXPHOS_genes_final.clean4.bed"
POOL_FILE <- "/work/cyu/poolseq/PPalign_output/fst/poolsize.txt"

OUTDIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/maxSNP_PBS"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# =========================
# population design
# =========================
AK_FW <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_FW <- c("SWA","THE","JOE","BEA","MUC","PYE","BOOT","ECHO","LAW","GOS","ROB")

design <- rbindlist(list(
  data.table(focal = AK_FW, region = "AK", marine = "RS",  outgroup = "SAY"),
  data.table(focal = BC_FW, region = "BC", marine = "SAY", outgroup = "RS")
))

# =========================
# read population names
# =========================
pool <- fread(POOL_FILE, header = FALSE)
setnames(pool, c("pop", "pool_size"))

pops <- pool$pop
cat("[info] populations from poolsize:\n")
print(pops)

# =========================
# sync parser
# sync columns:
# chr pos ref sample1 sample2 ...
# sample format: A:T:C:G:N:del
# =========================
parse_sync_counts <- function(x) {
  sp <- strsplit(x, ":", fixed = TRUE)[[1]]
  vals <- as.numeric(sp[1:4]) # A T C G
  vals[is.na(vals)] <- 0
  vals
}

get_major_minor_freq <- function(count1, count2) {
  total <- count1 + count2
  major <- which.max(total)
  minor <- order(total, decreasing = TRUE)[2]

  n1 <- sum(count1)
  n2 <- sum(count2)

  if (n1 < 4 || n2 < 4) return(c(NA_real_, NA_real_))

  p1 <- count1[minor] / n1
  p2 <- count2[minor] / n2

  c(p1, p2)
}

fst_pair_simple <- function(p1, p2) {
  if (!is.finite(p1) || !is.finite(p2)) return(NA_real_)

  pbar <- (p1 + p2) / 2
  ht <- 2 * pbar * (1 - pbar)
  hs <- (2 * p1 * (1 - p1) + 2 * p2 * (1 - p2)) / 2

  if (!is.finite(ht) || ht <= 0) return(NA_real_)

  fst <- (ht - hs) / ht
  fst <- max(fst, 0)
  fst <- min(fst, 0.999999)
  fst
}

calc_pair_fst_from_counts <- function(counts_list, popA, popB, pop_index) {
  i <- pop_index[[popA]]
  j <- pop_index[[popB]]

  if (is.null(i) || is.null(j)) return(NA_real_)

  f <- get_major_minor_freq(counts_list[[i]], counts_list[[j]])
  fst_pair_simple(f[1], f[2])
}

calc_pbs <- function(fst_AB, fst_AC, fst_BC) {
  if (!all(is.finite(c(fst_AB, fst_AC, fst_BC)))) return(NA_real_)

  fst_AB <- min(max(fst_AB, 0), 0.999999)
  fst_AC <- min(max(fst_AC, 0), 0.999999)
  fst_BC <- min(max(fst_BC, 0), 0.999999)

  T_AB <- -log(1 - fst_AB)
  T_AC <- -log(1 - fst_AC)
  T_BC <- -log(1 - fst_BC)

  (T_AB + T_AC - T_BC) / 2
}

# =========================
# prepare gene regions
# =========================
gene_bed <- fread(GENE_BED, header = FALSE)
gene_bed <- gene_bed[, .(
  chr = V1,
  start = as.integer(V2) + 1,  # BED 0-based to 1-based
  end = as.integer(V3),
  gene = V4
)]

setkey(gene_bed, chr, start, end)

# =========================
# read sync
# =========================
cat("[read] sync file\n")
sync <- fread(SYNC_FILE, header = FALSE)

# expected columns: chr pos ref + samples
setnames(sync, c("chr", "pos", "ref", pops))

sync[, start := as.integer(pos)]
sync[, end := as.integer(pos)]

setkey(sync, chr, start, end)

# =========================
# intersect sync with gene BED
# =========================
cat("[intersect] SNPs with genes\n")
snps_gene <- foverlaps(
  sync,
  gene_bed,
  by.x = c("chr", "start", "end"),
  by.y = c("chr", "start", "end"),
  type = "within",
  nomatch = 0
)

cat("[info] SNP-gene rows:", nrow(snps_gene), "\n")

pop_index <- as.list(seq_along(pops))
names(pop_index) <- pops

# =========================
# calculate SNP-level PBS
# =========================
res_list <- vector("list", nrow(snps_gene))

cat("[run] SNP-level PBS\n")

for (r in seq_len(nrow(snps_gene))) {
  if (r %% 10000 == 0) cat("[progress]", r, "/", nrow(snps_gene), "\n")

  row <- snps_gene[r]

  counts_list <- lapply(pops, function(p) parse_sync_counts(row[[p]]))

  tmp <- rbindlist(lapply(seq_len(nrow(design)), function(i) {
    d <- design[i]

    fst_AB <- calc_pair_fst_from_counts(counts_list, d$focal, d$marine, pop_index)
    fst_AC <- calc_pair_fst_from_counts(counts_list, d$focal, d$outgroup, pop_index)
    fst_BC <- calc_pair_fst_from_counts(counts_list, d$marine, d$outgroup, pop_index)

    pbs <- calc_pbs(fst_AB, fst_AC, fst_BC)

    data.table(
      gene = row$gene,
      chr = row$chr,
      pos = row$pos,
      focal = d$focal,
      region = d$region,
      marine = d$marine,
      outgroup = d$outgroup,
      fst_focal_marine = fst_AB,
      fst_focal_outgroup = fst_AC,
      fst_marine_outgroup = fst_BC,
      snp_PBS = pbs
    )
  }))

  res_list[[r]] <- tmp
}

snp_pbs <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

snp_pbs <- snp_pbs[is.finite(snp_PBS)]

fwrite(
  snp_pbs,
  file.path(OUTDIR, "nuclear_gene_snp_level_PBS.tsv.gz"),
  sep = "\t"
)

# =========================
# max SNP PBS per gene x population
# =========================
max_pbs <- snp_pbs[
  ,
  .SD[which.max(snp_PBS)],
  by = .(gene, focal, region, marine, outgroup)
]

setnames(max_pbs, "snp_PBS", "max_snp_PBS")

fwrite(
  max_pbs,
  file.path(OUTDIR, "nuclear_gene_maxSNP_PBS_by_population.tsv"),
  sep = "\t"
)

# =========================
# mean and max summary per gene
# =========================
gene_summary <- max_pbs[
  ,
  .(
    n_pop = .N,
    mean_maxSNP_PBS = mean(max_snp_PBS, na.rm = TRUE),
    median_maxSNP_PBS = median(max_snp_PBS, na.rm = TRUE),
    max_observed_PBS = max(max_snp_PBS, na.rm = TRUE)
  ),
  by = gene
][order(-mean_maxSNP_PBS)]

fwrite(
  gene_summary,
  file.path(OUTDIR, "gene_summary_maxSNP_PBS.tsv"),
  sep = "\t"
)

cat("\n[OK] done\n")
cat("[outdir] ", OUTDIR, "\n", sep = "")


#perm corrrelation
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

MAX_FILE <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/maxSNP_PBS/nuclear_gene_maxSNP_PBS_by_population.tsv"
MT_FILE  <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/pbsn1_gene_level/mt_noDloop_PBSn1_by_population.tsv"
OUTDIR   <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/maxSNP_PBS"

nu <- fread(MAX_FILE)
mt <- fread(MT_FILE)

mt <- mt[, .(focal, region, marine, outgroup, mt_PBS, mt_PBSn1)]

dt <- merge(
  nu,
  mt,
  by = c("focal", "region", "marine", "outgroup"),
  all.x = TRUE
)

set.seed(123)

perm_test <- function(x, y, nperm = 2000) {
  r_obs <- cor(x, y)
  r_perm <- replicate(nperm, cor(x, sample(y)))
  p <- mean(abs(r_perm) >= abs(r_obs))
  list(r = r_obs, p = p)
}

# ============================================================
# max SNP PBS vs mt PBS
# ============================================================
res_mtPBS <- dt[
  is.finite(max_snp_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 5 && var(max_snp_PBS) > 0 && var(mt_PBS) > 0) {
      pt <- perm_test(max_snp_PBS, mt_PBS, nperm = 2000)
      .(
        n_pop = .N,
        cor_r = pt$r,
        p_perm = pt$p,
        mean_maxSNP_PBS = mean(max_snp_PBS),
        mean_mt_PBS = mean(mt_PBS)
      )
    } else {
      .(
        n_pop = .N,
        cor_r = NA_real_,
        p_perm = NA_real_,
        mean_maxSNP_PBS = mean(max_snp_PBS),
        mean_mt_PBS = mean(mt_PBS)
      )
    }
  },
  by = gene
]

res_mtPBS[, p_bh := p.adjust(p_perm, method = "BH")]
setorder(res_mtPBS, p_perm, -cor_r)

fwrite(
  res_mtPBS,
  file.path(OUTDIR, "gene_correlation_maxSNP_PBS_vs_mtPBS_perm.tsv"),
  sep = "\t"
)

# ============================================================
# max SNP PBS vs mt PBSn1
# ============================================================
res_mtPBSn1 <- dt[
  is.finite(max_snp_PBS) & is.finite(mt_PBSn1),
  {
    if (.N >= 5 && var(max_snp_PBS) > 0 && var(mt_PBSn1) > 0) {
      pt <- perm_test(max_snp_PBS, mt_PBSn1, nperm = 2000)
      .(
        n_pop = .N,
        cor_r = pt$r,
        p_perm = pt$p,
        mean_maxSNP_PBS = mean(max_snp_PBS),
        mean_mt_PBSn1 = mean(mt_PBSn1)
      )
    } else {
      .(
        n_pop = .N,
        cor_r = NA_real_,
        p_perm = NA_real_,
        mean_maxSNP_PBS = mean(max_snp_PBS),
        mean_mt_PBSn1 = mean(mt_PBSn1)
      )
    }
  },
  by = gene
]

res_mtPBSn1[, p_bh := p.adjust(p_perm, method = "BH")]
setorder(res_mtPBSn1, p_perm, -cor_r)

fwrite(
  res_mtPBSn1,
  file.path(OUTDIR, "gene_correlation_maxSNP_PBS_vs_mtPBSn1_perm.tsv"),
  sep = "\t"
)

cat("\n=== max SNP PBS vs mt PBS permutation top ===\n")
print(head(res_mtPBS, 20))

cat("\n=== max SNP PBS vs mt PBSn1 permutation top ===\n")
print(head(res_mtPBSn1, 20))

cat("\n[write]\n")
cat(file.path(OUTDIR, "gene_correlation_maxSNP_PBS_vs_mtPBS_perm.tsv"), "\n")
cat(file.path(OUTDIR, "gene_correlation_maxSNP_PBS_vs_mtPBSn1_perm.tsv"), "\n")
cd /mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/maxSNP_PBS
awk -F'\t' 'NR>1 && $3 > 0 && $7 < 0.05' gene_correlation_maxSNP_PBS_vs_mtPBS_perm.tsv
ndufb1	20	0.646795681979088	0.001	0.566149764433928	0.133429270045148	0.0415
ndufa1	20	0.635082151053412	0.001	0.382942279674176	0.133429270045148	0.0415
cox7a2l	20	0.607201679558001	0.0015	1.7204054734779	0.133429270045148	0.0415







#pbsn1 shared ak/bc
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ==============================================================================
# 1. Paths
# ==============================================================================
INPUT_FILE <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/pbsn1_gene_level/merged_mt_nuclear_PBSn1_by_gene_population.tsv"
ANNOT_FILE <- "/work/cyu/codeml_sites_summary.merged.tsv"

OUTDIR <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/validation_results_subunit_only"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("[Loading] Reading merged PBS data...\n")
dt <- fread(INPUT_FILE)

cat("[Loading] Reading gene annotation...\n")
annot <- fread(ANNOT_FILE)
annot_gene <- unique(annot[, .(gene, role, complex)])

# ==============================================================================
# 2. Calculate region-wise mean PBS
# ==============================================================================
gene_summary <- dt[is.finite(nu_PBS), .(
  mean_PBS = mean(nu_PBS, na.rm = TRUE),
  n_pop = .N
), by = .(gene, region)]

ak_dt <- gene_summary[region == "AK", .(
  gene,
  mean_PBS_AK = mean_PBS,
  n_pop_AK = n_pop
)]

bc_dt <- gene_summary[region == "BC", .(
  gene,
  mean_PBS_BC = mean_PBS,
  n_pop_BC = n_pop
)]

compare_dt <- merge(ak_dt, bc_dt, by = "gene")

# ==============================================================================
# 3. Add annotation and keep subunits only
# ==============================================================================
compare_dt <- merge(compare_dt, annot_gene, by = "gene", all.x = TRUE)
compare_dt <- compare_dt[role == "subunit"]

cat("[Info] Number of subunit genes:", nrow(compare_dt), "\n")

# ==============================================================================
# 4. Define Shared / Region-specific status within subunits
# ==============================================================================
cutoff_ak <- quantile(compare_dt$mean_PBS_AK, 0.95, na.rm = TRUE)
cutoff_bc <- quantile(compare_dt$mean_PBS_BC, 0.95, na.rm = TRUE)

cat("[Info] AK 95% cutoff:", cutoff_ak, "\n")
cat("[Info] BC 95% cutoff:", cutoff_bc, "\n")

compare_dt[, status := "Background"]

compare_dt[
  mean_PBS_AK > cutoff_ak & mean_PBS_BC > cutoff_bc,
  status := "Shared_Selection"
]

compare_dt[
  mean_PBS_AK > cutoff_ak & mean_PBS_BC <= cutoff_bc,
  status := "AK_Specific"
]

compare_dt[
  mean_PBS_AK <= cutoff_ak & mean_PBS_BC > cutoff_bc,
  status := "BC_Specific"
]

cat("[Info] Status counts:\n")
print(compare_dt[, .N, by = status])

# ==============================================================================
# 5. Save comparison table
# ==============================================================================
OUT_TABLE <- file.path(OUTDIR, "subunit_only_shared_vs_specific_comparison.tsv")
fwrite(compare_dt, OUT_TABLE, sep = "\t")
cat("[Write] Table:", OUT_TABLE, "\n")

# ==============================================================================
# 6. Plot: subunit only; label all non-background genes
# ==============================================================================
cat("[Plotting] Generating subunit-only comparison scatter plot...\n")

label_dt <- compare_dt[status != "Background"]

p1 <- ggplot(compare_dt, aes(x = mean_PBS_BC, y = mean_PBS_AK, color = status)) +
  geom_point(size = 2.2, alpha = 0.7) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_text_repel(
    data = label_dt,
    aes(label = gene),
    size = 3,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "AK_Specific" = "red",
    "BC_Specific" = "blue",
    "Shared_Selection" = "purple",
    "Background" = "grey80"
  )) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black")
  ) +
  labs(
    x = "Mean PBS (British Columbia)",
    y = "Mean PBS (Alaska)",
    color = "Status"
  )

OUT_PNG <- file.path(OUTDIR, "Fig_subunit_only_AK_vs_BC_PBS_Comparison.png")
OUT_PDF <- file.path(OUTDIR, "Fig_subunit_only_AK_vs_BC_PBS_Comparison.pdf")

ggsave(OUT_PNG, p1, width = 8, height = 7, dpi = 300)
ggsave(OUT_PDF, p1, width = 8, height = 7)

cat("[Write] Figure PNG:", OUT_PNG, "\n")
cat("[Write] Figure PDF:", OUT_PDF, "\n")

# ==============================================================================
# 7. Merge mitochondrial concordance validation
# ==============================================================================
cat("[Validation] Checking concordance with mitochondrial trajectory...\n")

COR_FILE <- "/mnt/spareHD_2/nuclear_marked_duplicates/gene_fst_work/pbsn1_gene_level/gene_correlation_nuPBSn1_vs_mtPBSn1.tsv"
cor_dt <- fread(COR_FILE)

final_dt <- merge(
  compare_dt,
  cor_dt[, .(gene, cor_r, p_bh)],
  by = "gene",
  all.x = TRUE
)

true_signals <- final_dt[
  status == "AK_Specific" & (is.na(cor_r) | abs(cor_r) < 0.5)
]

OUT_TRUE <- file.path(OUTDIR, "subunit_only_verified_AK_specific_genes.tsv")
fwrite(true_signals, OUT_TRUE, sep = "\t")

cat("[Write] Verified AK genes:", OUT_TRUE, "\n")

cat("\n[Done] Analysis complete.\n")
print(p1)