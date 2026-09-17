
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# Input
# ============================================================
DELTA_FILE <- "/path/to/workspace/ldx_all_subunits/ld/deltaAF_long.noAMO_noLB.ldPruned_kept.keepMarine.R2_0.2.tsv.gz"

# temporary lineage file; replace later with UPGMA-defined mt lineage
CLUSTER_FILE <- "/path/to/workspace/mt_lineage_for_glm.tsv"

OUTDIR <- "/path/to/workspace/ldx_all_subunits/ld/pairwise_deltaAF_similarity_ldpruned_effectFiltered_oneSNPperGene_10000perm"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parameters
# ============================================================
N_PERM <- 10000
set.seed(123)

MIN_DEPTH <- 20
EFFECT_THRESHOLD <- 0.20
MIN_SNPS_MATRIX <- 10

# ============================================================
# Functions
# ============================================================
normalize_pop <- function(x){
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl = TRUE)
}

perm_test <- function(dist_obj, cl, nperm = 10000){

  M <- as.matrix(dist_obj)
  idx <- which(upper.tri(M), arr.ind = TRUE)

  d <- M[upper.tri(M)]
  same <- cl[idx[,1]] == cl[idx[,2]]

  if(sum(same, na.rm = TRUE) == 0 || sum(!same, na.rm = TRUE) == 0){
    return(list(
      obs = NA_real_,
      p_perm = NA_real_,
      within = NA_real_,
      between = NA_real_,
      n_within_pairs = sum(same, na.rm = TRUE),
      n_between_pairs = sum(!same, na.rm = TRUE),
      n_perm_used = 0L,
      note = "Not enough within/between pairs"
    ))
  }

  within_obs  <- mean(d[same], na.rm = TRUE)
  between_obs <- mean(d[!same], na.rm = TRUE)

  # obs < 0 means within-lineage populations are more similar
  obs <- within_obs - between_obs

  ge <- 0L
  used <- 0L

  for(b in seq_len(nperm)){
    clp <- sample(cl)
    samep <- clp[idx[,1]] == clp[idx[,2]]

    if(sum(samep, na.rm = TRUE) == 0 || sum(!samep, na.rm = TRUE) == 0) next

    within_p  <- mean(d[samep], na.rm = TRUE)
    between_p <- mean(d[!samep], na.rm = TRUE)

    statp <- within_p - between_p

    if(!is.finite(statp)) next

    used <- used + 1L

    # one-sided: is observed within-between smaller than random?
    if(statp <= obs) ge <- ge + 1L
  }

  p <- if(used > 0) (1 + ge) / (1 + used) else NA_real_

  list(
    obs = obs,
    p_perm = p,
    within = within_obs,
    between = between_obs,
    n_within_pairs = sum(same, na.rm = TRUE),
    n_between_pairs = sum(!same, na.rm = TRUE),
    n_perm_used = used,
    note = "obs = within - between; negative means within mt lineage more similar"
  )
}

make_pairwise_table <- function(mat, cl, region_name, metric_name){

  cor_mat <- cor(t(mat), use = "pairwise.complete.obs")
  dist_mat <- 1 - cor_mat

  M <- as.matrix(dist_mat)
  idx <- which(upper.tri(M), arr.ind = TRUE)

  data.table(
    region = region_name,
    metric = metric_name,
    pop1 = rownames(M)[idx[,1]],
    pop2 = colnames(M)[idx[,2]],
    lineage1 = as.character(cl[idx[,1]]),
    lineage2 = as.character(cl[idx[,2]]),
    dist = M[upper.tri(M)],
    type = ifelse(cl[idx[,1]] == cl[idx[,2]], "Within", "Between")
  )
}

run_one_metric <- function(mat, cl, region_name, metric_name){

  if(metric_name == "Signed_deltaAF"){
    use_mat <- mat
  } else if(metric_name == "Abs_deltaAF"){
    use_mat <- abs(mat)
  } else {
    stop("Unknown metric")
  }

  cor_mat <- cor(t(use_mat), use = "pairwise.complete.obs")
  dist_obj <- as.dist(1 - cor_mat)

  out <- perm_test(dist_obj, cl, nperm = N_PERM)

  summary <- data.table(
    region = region_name,
    metric = metric_name,
    obs = out$obs,
    p_perm = out$p_perm,
    within = out$within,
    between = out$between,
    n_within_pairs = out$n_within_pairs,
    n_between_pairs = out$n_between_pairs,
    n_perm_used = out$n_perm_used,
    note = out$note
  )

  pairwise <- make_pairwise_table(
    use_mat,
    cl,
    region_name,
    metric_name
  )

  list(summary = summary, pairwise = pairwise)
}

# ============================================================
# Read deltaAF
# ============================================================
cat("[INFO] Reading deltaAF file...\n")
DEL <- fread(DELTA_FILE)

if ("Pop" %in% names(DEL)) setnames(DEL, "Pop", "pop")
if ("Gene" %in% names(DEL)) setnames(DEL, "Gene", "gene")

DEL[, pop := normalize_pop(pop)]
DEL[, gene := tolower(as.character(gene))]

if (!"deltaAF" %in% names(DEL)) {
  stop("Cannot find deltaAF column.")
}

if (!"snp" %in% names(DEL)) {
  if (all(c("chr", "pos") %in% names(DEL))) {
    DEL[, snp := paste(chr, pos, sep = ":")]
  } else {
    stop("Need either snp column or chr/pos columns.")
  }
}

if (!"region" %in% names(DEL)) {
  DEL[, region := "ALL"]
}

if (!all(c("af", "depth") %in% names(DEL))) {
  stop("Need af and depth columns.")
}

DEL[, af := as.numeric(af)]
DEL[, depth := as.numeric(depth)]
DEL[, deltaAF := as.numeric(deltaAF)]

DEL <- DEL[
  is.finite(deltaAF) &
  is.finite(af) &
  is.finite(depth) &
  af >= 0 & af <= 1 &
  depth >= MIN_DEPTH
]

cat("[INFO] rows after depth filter:", nrow(DEL), "\n")
cat("[INFO] populations:", uniqueN(DEL$pop), "\n")
cat("[INFO] genes:", uniqueN(DEL$gene), "\n")
cat("[INFO] SNPs:", uniqueN(DEL$snp), "\n")

# ============================================================
# Read mt lineage
# ============================================================
CL <- fread(CLUSTER_FILE)
setnames(CL, names(CL), tolower(names(CL)))

if (!all(c("pop", "mt_lineage") %in% names(CL))) {
  stop("CLUSTER_FILE must contain columns: pop, mt_lineage")
}

CL[, pop := normalize_pop(pop)]
CL[, mt_lineage := as.factor(mt_lineage)]
CL <- unique(CL, by = "pop")

cat("[INFO] mt lineages loaded:\n")
print(CL[, .N, by = mt_lineage][order(mt_lineage)])

# ============================================================
# Merge
# ============================================================
DEL <- merge(DEL, CL, by = "pop", all.x = TRUE)
DEL <- DEL[!is.na(mt_lineage)]

cat("[INFO] rows after lineage merge:", nrow(DEL), "\n")
cat("[INFO] populations after merge:", uniqueN(DEL$pop), "\n")
print(DEL[, .N, by = .(pop, mt_lineage)][order(mt_lineage, pop)])

# ============================================================
# Calculate max lineage AF difference per region x gene x SNP
# ============================================================
cat("\n[INFO] Calculating max lineage AF difference per SNP...\n")

DEL[, focal_count := round(af * depth)]
DEL[, focal_count := pmax(0, pmin(focal_count, depth))]
DEL[, other_count := depth - focal_count]

lineage_af <- DEL[, .(
  lineage_focal_count = sum(focal_count, na.rm = TRUE),
  lineage_depth = sum(depth, na.rm = TRUE),
  n_pops_lineage = uniqueN(pop)
), by = .(region, gene, snp, mt_lineage)]

lineage_af <- lineage_af[lineage_depth > 0]
lineage_af[, lineage_AF := lineage_focal_count / lineage_depth]

snp_effect <- lineage_af[, .(
  n_lineages = uniqueN(mt_lineage),
  max_lineage_AF = max(lineage_AF, na.rm = TRUE),
  min_lineage_AF = min(lineage_AF, na.rm = TRUE),
  max_lineage_AF_diff = max(lineage_AF, na.rm = TRUE) - min(lineage_AF, na.rm = TRUE),
  max_lineage = as.character(mt_lineage[which.max(lineage_AF)][1]),
  min_lineage = as.character(mt_lineage[which.min(lineage_AF)][1])
), by = .(region, gene, snp)]

snp_effect[, pass_effect_filter := n_lineages >= 2 & max_lineage_AF_diff >= EFFECT_THRESHOLD]

# ============================================================
# One SNP per gene:
# within each region x gene, keep SNP with largest max_lineage_AF_diff
# tie-breaker: largest mean abs deltaAF, then SNP name
# ============================================================
mean_abs_daf <- DEL[, .(
  mean_abs_deltaAF = mean(abs(deltaAF), na.rm = TRUE),
  max_abs_deltaAF = max(abs(deltaAF), na.rm = TRUE)
), by = .(region, gene, snp)]

snp_effect <- merge(
  snp_effect,
  mean_abs_daf,
  by = c("region", "gene", "snp"),
  all.x = TRUE
)

candidate_snps <- snp_effect[pass_effect_filter == TRUE]

one_snp_per_gene <- candidate_snps[
  order(region, gene, -max_lineage_AF_diff, -mean_abs_deltaAF, snp)
][
  , .SD[1], by = .(region, gene)
]

selected_out <- file.path(
  OUTDIR,
  paste0("selected_oneSNPperGene_maxLineageAFdiff_effect", EFFECT_THRESHOLD, ".tsv")
)

fwrite(one_snp_per_gene, selected_out, sep = "\t")
cat("[OK] wrote:", selected_out, "\n")

cat("\n[INFO] Selected SNP summary:\n")
print(one_snp_per_gene[, .(
  n_genes_selected = uniqueN(gene),
  n_snps_selected = uniqueN(snp),
  median_max_lineage_AF_diff = median(max_lineage_AF_diff, na.rm = TRUE),
  max_max_lineage_AF_diff = max(max_lineage_AF_diff, na.rm = TRUE)
), by = region])

# Keep only selected SNPs
KEEP <- one_snp_per_gene[, .(region, gene, snp)]

DEL_filt <- merge(
  DEL,
  KEEP,
  by = c("region", "gene", "snp"),
  all = FALSE
)

cat("\n[INFO] rows after one-SNP-per-gene filter:", nrow(DEL_filt), "\n")
cat("[INFO] selected SNPs by region:\n")
print(DEL_filt[, .(
  n_genes = uniqueN(gene),
  n_snps = uniqueN(snp),
  n_rows = .N
), by = region])

# ============================================================
# Pairwise profile similarity within each region
# ============================================================
ALL_PAIRWISE <- list()
SUMMARY_LIST <- list()

regions <- sort(unique(DEL_filt$region))

for (r in regions) {

  cat("\n================ Region:", r, "================\n")

  X <- DEL_filt[region == r]

  W <- dcast(
    X,
    pop + mt_lineage ~ snp,
    value.var = "deltaAF",
    fun.aggregate = mean,
    fill = NA_real_
  )

  if (nrow(W) < 4) {
    cat("[WARN] too few populations; skip region:", r, "\n")
    next
  }

  mat <- as.matrix(W[, -(1:2)])
  rownames(mat) <- W$pop
  cl <- W$mt_lineage

  # Remove SNP columns with too few finite values or zero variance
  keep_cols <- apply(mat, 2, function(z) {
    sum(is.finite(z)) >= 3 && sd(z, na.rm = TRUE) > 0
  })

  mat <- mat[, keep_cols, drop = FALSE]

  cat("[INFO] populations:", nrow(mat), "\n")
  cat("[INFO] selected one-SNP-per-gene SNPs retained in matrix:", ncol(mat), "\n")
  print(data.table(pop = rownames(mat), mt_lineage = cl)[order(mt_lineage, pop)])

  if (ncol(mat) < MIN_SNPS_MATRIX) {
    cat("[WARN] too few SNPs after filtering; skip region:", r, "\n")
    next
  }

  signed_res <- run_one_metric(
    mat = mat,
    cl = cl,
    region_name = r,
    metric_name = "Signed_deltaAF"
  )

  abs_res <- run_one_metric(
    mat = mat,
    cl = cl,
    region_name = r,
    metric_name = "Abs_deltaAF"
  )

  SUMMARY_LIST[[paste(r, "Signed_deltaAF", sep = "__")]] <- signed_res$summary
  SUMMARY_LIST[[paste(r, "Abs_deltaAF", sep = "__")]] <- abs_res$summary

  ALL_PAIRWISE[[paste(r, "Signed_deltaAF", sep = "__")]] <- signed_res$pairwise
  ALL_PAIRWISE[[paste(r, "Abs_deltaAF", sep = "__")]] <- abs_res$pairwise
}

PAIRWISE <- rbindlist(ALL_PAIRWISE, fill = TRUE)
SUMMARY <- rbindlist(SUMMARY_LIST, fill = TRUE)

pairwise_out <- file.path(
  OUTDIR,
  paste0("pairwise_deltaAF_similarity_within_between_effect", EFFECT_THRESHOLD, "_oneSNPperGene_10000perm.tsv")
)

summary_out <- file.path(
  OUTDIR,
  paste0("pairwise_deltaAF_similarity_summary_effect", EFFECT_THRESHOLD, "_oneSNPperGene_10000perm.tsv")
)

fwrite(PAIRWISE, pairwise_out, sep = "\t")
fwrite(SUMMARY, summary_out, sep = "\t")

cat("\n[OK] wrote:", pairwise_out, "\n")
cat("[OK] wrote:", summary_out, "\n")

cat("\nSummary:\n")
print(SUMMARY)

# ============================================================
# Plot functions
# ============================================================
plot_bar <- function(metric_name, ylab, outfile){

  P <- PAIRWISE[metric == metric_name]
  P[, type := factor(type, levels = c("Within", "Between"))]

  BAR <- P[, .(
    mean_dist = mean(dist, na.rm = TRUE),
    se_dist = sd(dist, na.rm = TRUE) / sqrt(.N),
    n_pairs = .N
  ), by = .(region, type)]

  LAB <- SUMMARY[metric == metric_name]
  LAB[, label := ifelse(
    p_perm < 0.05,
    paste0("p = ", signif(p_perm, 2), " *"),
    paste0("p = ", signif(p_perm, 2))
  )]

  ypos <- BAR[, .(
    y = max(mean_dist + se_dist, na.rm = TRUE) * 1.10
  ), by = region]

  LAB <- merge(LAB, ypos, by = "region", all.x = TRUE)

  p <- ggplot(BAR, aes(x = type, y = mean_dist, fill = type)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(
      aes(ymin = mean_dist - se_dist, ymax = mean_dist + se_dist),
      width = 0.18,
      linewidth = 0.7
    ) +
    geom_text(
      data = LAB,
      aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      size = 4
    ) +
    facet_wrap(~region, scales = "free_y") +
    theme_classic(base_size = 14) +
    labs(
      x = NULL,
      y = ylab,
      title = paste0("One SNP per gene; max lineage AF difference ≥ ", EFFECT_THRESHOLD)
    ) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(outfile, p, width = 6.8, height = 4.6, dpi = 300)
  cat("[OK] saved:", outfile, "\n")
}

plot_box <- function(metric_name, ylab, outfile){

  P <- PAIRWISE[metric == metric_name]
  P[, type := factor(type, levels = c("Within", "Between"))]

  p <- ggplot(P, aes(x = type, y = dist, fill = type)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.65) +
    geom_jitter(width = 0.12, alpha = 0.4, size = 1) +
    facet_wrap(~region, scales = "free_y") +
    theme_classic(base_size = 14) +
    labs(
      x = NULL,
      y = ylab,
      title = paste0("One SNP per gene; max lineage AF difference ≥ ", EFFECT_THRESHOLD)
    ) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(outfile, p, width = 6.8, height = 4.6, dpi = 300)
  cat("[OK] saved:", outfile, "\n")
}

plot_bar(
  "Signed_deltaAF",
  "1 - correlation of ΔAF profiles",
  file.path(OUTDIR, paste0("Fig_pairwise_similarity_barplot_signed_effect", EFFECT_THRESHOLD, "_oneSNPperGene_10000perm.png"))
)

plot_box(
  "Signed_deltaAF",
  "1 - correlation of ΔAF profiles",
  file.path(OUTDIR, paste0("Fig_pairwise_similarity_boxplot_signed_effect", EFFECT_THRESHOLD, "_oneSNPperGene_10000perm.png"))
)

plot_bar(
  "Abs_deltaAF",
  "1 - correlation of |ΔAF| profiles",
  file.path(OUTDIR, paste0("Fig_pairwise_similarity_barplot_abs_effect", EFFECT_THRESHOLD, "_oneSNPperGene_10000perm.png"))
)

plot_box(
  "Abs_deltaAF",
  "1 - correlation of |ΔAF| profiles",
  file.path(OUTDIR, paste0("Fig_pairwise_similarity_boxplot_abs_effect", EFFECT_THRESHOLD, "_oneSNPperGene_10000perm.png"))
)

cat("\nDONE\n")

