#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

set.seed(20260810)

erc_dir <- paste0(
  "/mnt/spareHD_2/genomewide_codeml_kuster/",
  "10_erc_results/mtPCG_composite_t_spearman"
)

nuclear_file <- file.path(
  erc_dir,
  "nuclear_relative_rates.tsv.gz"
)

mt_file <- file.path(
  erc_dir,
  "mt13_relative_rates.tsv"
)

annotation_file <- file.path(
  erc_dir,
  "erc_mtPCG_composite_annotated.tsv"
)

n_permutations <- 100000L
n_bootstrap <- 10000L

# ----------------------------------------------------------
# Read matrices
# ----------------------------------------------------------

read_rate_matrix <- function(path) {

  dat <- fread(path)

  id_candidates <- c(
    "nuclear_gene",
    "gene_id",
    "gene",
    "mt_gene",
    "Gene",
    names(dat)[1]
  )

  id_column <- id_candidates[
    id_candidates %in% names(dat)
  ][1]

  if (is.na(id_column)) {
    stop("Cannot identify ID column: ", path)
  }

  ids <- as.character(dat[[id_column]])
  rate_columns <- setdiff(names(dat), id_column)

  for (column_name in rate_columns) {
    set(
      dat,
      j = column_name,
      value = suppressWarnings(
        as.numeric(dat[[column_name]])
      )
    )
  }

  result <- as.matrix(dat[, ..rate_columns])
  storage.mode(result) <- "double"
  rownames(result) <- ids

  result
}

nuclear <- read_rate_matrix(nuclear_file)
mt <- read_rate_matrix(mt_file)

rownames(mt) <- toupper(
  sub("^MT[-_]", "", rownames(mt))
)

common_branches <- intersect(
  colnames(nuclear),
  colnames(mt)
)

nuclear <- nuclear[, common_branches, drop = FALSE]
mt <- mt[, common_branches, drop = FALSE]

cat(
  "[load] nuclear =", nrow(nuclear),
  "mt =", nrow(mt),
  "branches =", length(common_branches),
  "\n"
)

# ----------------------------------------------------------
# Annotation
# ----------------------------------------------------------

annotation <- fread(annotation_file)

id_candidates <- c(
  "nuclear_gene",
  "gene_id",
  "stickleback_gene",
  "stickleback_gene_id"
)

id_column <- id_candidates[
  id_candidates %in% names(annotation)
][1]

if (is.na(id_column)) {
  stop("Cannot identify annotation gene-ID column")
}

setnames(annotation, id_column, "analysis_gene_id")
annotation <- unique(annotation, by = "analysis_gene_id")

for (column_name in c(
  "own_role",
  "own_complex",
  "new_symbol",
  "old_name"
)) {
  if (!column_name %in% names(annotation)) {
    annotation[, (column_name) := NA_character_]
  }
}

# ----------------------------------------------------------
# Helper functions
# ----------------------------------------------------------

safe_cor <- function(x, y) {

  keep <- is.finite(x) & is.finite(y)

  if (
    sum(keep) < 5 ||
    sd(x[keep]) == 0 ||
    sd(y[keep]) == 0
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x[keep],
      y[keep],
      method = "spearman"
    )
  )
}

composite <- function(rate_matrix, genes) {

  genes <- intersect(
    genes,
    rownames(rate_matrix)
  )

  if (length(genes) == 0) {
    stop("Composite requested for an empty gene set")
  }

  colMeans(
    rate_matrix[genes, , drop = FALSE],
    na.rm = TRUE
  )
}

clean <- function(x) {
  trimws(as.character(x))
}

# ----------------------------------------------------------
# Define CIV genes
# ----------------------------------------------------------

civ_nuclear_genes <- annotation[
  clean(own_role) == "subunit" &
    clean(own_complex) == "CIV" &
    analysis_gene_id %in% rownames(nuclear),
  analysis_gene_id
]

civ_mt_genes <- intersect(
  c("COX1", "COX2", "COX3"),
  rownames(mt)
)

if (length(civ_nuclear_genes) < 2) {
  stop("Fewer than two nuclear CIV genes found")
}

if (length(civ_mt_genes) != 3) {
  stop(
    "Expected COX1, COX2 and COX3; found: ",
    paste(civ_mt_genes, collapse = ", ")
  )
}

cat("\n===== Nuclear CIV genes =====\n")

civ_annotation <- annotation[
  analysis_gene_id %in% civ_nuclear_genes,
  .(
    analysis_gene_id,
    new_symbol,
    old_name,
    own_role,
    own_complex
  )
]

print(civ_annotation)

cat(
  "\nMitochondrial CIV genes:",
  paste(civ_mt_genes, collapse = ", "),
  "\n"
)

# ----------------------------------------------------------
# Full CIV ERC
# ----------------------------------------------------------

civ_nuclear_composite <- composite(
  nuclear,
  civ_nuclear_genes
)

civ_mt_composite <- composite(
  mt,
  civ_mt_genes
)

observed_rs <- safe_cor(
  civ_nuclear_composite,
  civ_mt_composite
)

cat("\n===== Full CIV ERC =====\n")
cat("Spearman rs:", observed_rs, "\n")

# ----------------------------------------------------------
# Leave one mitochondrial gene out
# ----------------------------------------------------------

loo_mt <- rbindlist(
  lapply(civ_mt_genes, function(excluded_gene) {

    retained <- setdiff(
      civ_mt_genes,
      excluded_gene
    )

    data.table(
      excluded_mt_gene = excluded_gene,
      retained_mt_genes = paste(
        retained,
        collapse = ","
      ),
      rs = safe_cor(
        civ_nuclear_composite,
        composite(mt, retained)
      )
    )
  })
)

cat("\n===== Leave-one-mt-gene-out =====\n")
print(loo_mt)

# Also inspect each mt CIV gene separately
single_mt <- rbindlist(
  lapply(civ_mt_genes, function(mt_gene) {

    data.table(
      mt_gene = mt_gene,
      rs = safe_cor(
        civ_nuclear_composite,
        as.numeric(mt[mt_gene, ])
      )
    )
  })
)

cat("\n===== Individual mt CIV genes =====\n")
print(single_mt)

# ----------------------------------------------------------
# Leave one nuclear gene out
# ----------------------------------------------------------

loo_nuclear <- rbindlist(
  lapply(civ_nuclear_genes, function(excluded_gene) {

    retained <- setdiff(
      civ_nuclear_genes,
      excluded_gene
    )

    gene_meta <- annotation[
      analysis_gene_id == excluded_gene
    ][1]

    data.table(
      excluded_nuclear_gene = excluded_gene,
      excluded_symbol = gene_meta$new_symbol,
      excluded_old_name = gene_meta$old_name,
      retained_nuclear_n = length(retained),
      rs = safe_cor(
        composite(nuclear, retained),
        civ_mt_composite
      )
    )
  })
)

cat("\n===== Leave-one-nuclear-gene-out =====\n")
print(loo_nuclear)

# Individual nuclear genes
single_nuclear <- rbindlist(
  lapply(civ_nuclear_genes, function(nuclear_gene) {

    gene_meta <- annotation[
      analysis_gene_id == nuclear_gene
    ][1]

    data.table(
      nuclear_gene = nuclear_gene,
      new_symbol = gene_meta$new_symbol,
      old_name = gene_meta$old_name,
      rs = safe_cor(
        as.numeric(nuclear[nuclear_gene, ]),
        civ_mt_composite
      )
    )
  })
)

setorder(single_nuclear, -rs)

cat("\n===== Individual nuclear CIV genes =====\n")
print(single_nuclear)

# ----------------------------------------------------------
# Leave one branch out
# ----------------------------------------------------------

loo_branch <- rbindlist(
  lapply(seq_along(common_branches), function(branch_index) {

    retained <- setdiff(
      seq_along(common_branches),
      branch_index
    )

    data.table(
      excluded_branch = common_branches[branch_index],
      rs = safe_cor(
        civ_nuclear_composite[retained],
        civ_mt_composite[retained]
      )
    )
  })
)

cat("\n===== Leave-one-branch-out summary =====\n")

print(
  loo_branch[
    ,
    .(
      n = .N,
      minimum_rs = min(rs, na.rm = TRUE),
      median_rs = median(rs, na.rm = TRUE),
      maximum_rs = max(rs, na.rm = TRUE),
      all_positive = all(rs > 0, na.rm = TRUE)
    )
  ]
)

cat("\nLowest leave-one-branch ERC values:\n")
print(loo_branch[order(rs)][1:min(10, .N)])

# ----------------------------------------------------------
# Branch type: terminal versus internal
# ----------------------------------------------------------

branch_is_terminal <- function(branch_label) {

  numbers <- as.integer(
    unlist(
      regmatches(
        branch_label,
        gregexpr("[0-9]+", branch_label)
      )
    )
  )

  if (length(numbers) < 2) {
    return(NA)
  }

  # Tips are numbered 1 through 27.
  any(numbers <= 27)
}

branch_type <- vapply(
  common_branches,
  branch_is_terminal,
  logical(1)
)

if (anyNA(branch_type)) {
  warning(
    "Could not classify some branch labels; ",
    "stratified permutation will omit classification"
  )
}

cat("\n===== Branch classes =====\n")
cat("Terminal branches:", sum(branch_type, na.rm = TRUE), "\n")
cat("Internal branches:", sum(!branch_type, na.rm = TRUE), "\n")

# ----------------------------------------------------------
# Permutation tests
# ----------------------------------------------------------

unrestricted_null <- numeric(n_permutations)
stratified_null <- numeric(n_permutations)

terminal_indices <- which(branch_type)
internal_indices <- which(!branch_type)

cat(
  "\n[permutation] Running",
  n_permutations,
  "branch-label permutations\n"
)

for (iteration in seq_len(n_permutations)) {

  unrestricted_null[iteration] <- safe_cor(
    civ_nuclear_composite,
    sample(civ_mt_composite)
  )

  permuted_mt <- civ_mt_composite

  permuted_mt[terminal_indices] <- sample(
    civ_mt_composite[terminal_indices]
  )

  permuted_mt[internal_indices] <- sample(
    civ_mt_composite[internal_indices]
  )

  stratified_null[iteration] <- safe_cor(
    civ_nuclear_composite,
    permuted_mt
  )
}

permutation_summary <- data.table(
  permutation_type = c(
    "unrestricted",
    "terminal_internal_stratified"
  ),
  observed_rs = observed_rs,
  null_mean = c(
    mean(unrestricted_null, na.rm = TRUE),
    mean(stratified_null, na.rm = TRUE)
  ),
  null_sd = c(
    sd(unrestricted_null, na.rm = TRUE),
    sd(stratified_null, na.rm = TRUE)
  ),
  null_95_quantile = c(
    quantile(
      unrestricted_null,
      0.95,
      na.rm = TRUE,
      names = FALSE
    ),
    quantile(
      stratified_null,
      0.95,
      na.rm = TRUE,
      names = FALSE
    )
  ),
  empirical_p_upper = c(
    (
      1 + sum(
        unrestricted_null >= observed_rs,
        na.rm = TRUE
      )
    ) / (
      1 + sum(is.finite(unrestricted_null))
    ),
    (
      1 + sum(
        stratified_null >= observed_rs,
        na.rm = TRUE
      )
    ) / (
      1 + sum(is.finite(stratified_null))
    )
  )
)

cat("\n===== CIV branch-permutation tests =====\n")
print(permutation_summary)

# ----------------------------------------------------------
# Nmt-ARS versus cyto-ARS paired bootstrap
# ----------------------------------------------------------

nmt_ars <- annotation[
  clean(own_role) == "Nmt-ARS" &
    analysis_gene_id %in% rownames(nuclear),
  analysis_gene_id
]

cyto_ars <- annotation[
  clean(own_role) == "cyto-ARS" &
    analysis_gene_id %in% rownames(nuclear),
  analysis_gene_id
]

mt_all_genes <- rownames(mt)

observed_nmt_ars <- safe_cor(
  composite(nuclear, nmt_ars),
  composite(mt, mt_all_genes)
)

observed_cyto_ars <- safe_cor(
  composite(nuclear, cyto_ars),
  composite(mt, mt_all_genes)
)

observed_ars_difference <-
  observed_nmt_ars - observed_cyto_ars

ars_bootstrap <- data.table(
  iteration = seq_len(n_bootstrap),
  nmt_ars_rs = NA_real_,
  cyto_ars_rs = NA_real_,
  difference_rs = NA_real_
)

cat(
  "\n[bootstrap] Nmt-ARS vs cyto-ARS:",
  n_bootstrap,
  "iterations\n"
)

for (iteration in seq_len(n_bootstrap)) {

  # Use the same resampled mt genes in both correlations,
  # making the functional comparison paired.
  sampled_mt <- sample(
    mt_all_genes,
    length(mt_all_genes),
    replace = TRUE
  )

  sampled_nmt_ars <- sample(
    nmt_ars,
    length(nmt_ars),
    replace = TRUE
  )

  sampled_cyto_ars <- sample(
    cyto_ars,
    length(cyto_ars),
    replace = TRUE
  )

  mt_rate <- composite(
    mt,
    sampled_mt
  )

  nmt_value <- safe_cor(
    composite(nuclear, sampled_nmt_ars),
    mt_rate
  )

  cyto_value <- safe_cor(
    composite(nuclear, sampled_cyto_ars),
    mt_rate
  )

  ars_bootstrap[iteration, `:=`(
    nmt_ars_rs = nmt_value,
    cyto_ars_rs = cyto_value,
    difference_rs = nmt_value - cyto_value
  )]
}

ars_differences <- ars_bootstrap[
  is.finite(difference_rs),
  difference_rs
]

ars_comparison <- data.table(
  nmt_ars_n = length(nmt_ars),
  cyto_ars_n = length(cyto_ars),
  observed_nmt_ars_rs = observed_nmt_ars,
  observed_cyto_ars_rs = observed_cyto_ars,
  observed_difference_rs = observed_ars_difference,
  bootstrap_mean_difference = mean(ars_differences),
  difference_ci_lower = quantile(
    ars_differences,
    0.025,
    names = FALSE
  ),
  difference_ci_upper = quantile(
    ars_differences,
    0.975,
    names = FALSE
  ),
  bootstrap_p_nmt_greater =
    (
      1 + sum(ars_differences <= 0)
    ) /
    (
      1 + length(ars_differences)
    )
)

cat("\n===== Nmt-ARS versus cyto-ARS =====\n")
print(ars_comparison)

# ----------------------------------------------------------
# Save results
# ----------------------------------------------------------

fwrite(
  civ_annotation,
  file.path(erc_dir, "CIV_nuclear_gene_members.tsv"),
  sep = "\t"
)

fwrite(
  loo_mt,
  file.path(erc_dir, "CIV_leave_one_mt_gene_out.tsv"),
  sep = "\t"
)

fwrite(
  single_mt,
  file.path(erc_dir, "CIV_individual_mt_gene_ERC.tsv"),
  sep = "\t"
)

fwrite(
  loo_nuclear,
  file.path(erc_dir, "CIV_leave_one_nuclear_gene_out.tsv"),
  sep = "\t"
)

fwrite(
  single_nuclear,
  file.path(erc_dir, "CIV_individual_nuclear_gene_ERC.tsv"),
  sep = "\t"
)

fwrite(
  loo_branch,
  file.path(erc_dir, "CIV_leave_one_branch_out.tsv"),
  sep = "\t"
)

fwrite(
  permutation_summary,
  file.path(erc_dir, "CIV_branch_permutation_summary.tsv"),
  sep = "\t"
)

fwrite(
  data.table(
    iteration = seq_len(n_permutations),
    unrestricted_rs = unrestricted_null,
    stratified_rs = stratified_null
  ),
  file.path(
    erc_dir,
    "CIV_branch_permutation_null_100000.tsv.gz"
  ),
  sep = "\t"
)

fwrite(
  ars_comparison,
  file.path(erc_dir, "Nmt_ARS_vs_cyto_ARS_summary.tsv"),
  sep = "\t"
)

fwrite(
  ars_bootstrap,
  file.path(
    erc_dir,
    "Nmt_ARS_vs_cyto_ARS_bootstrap_10000.tsv.gz"
  ),
  sep = "\t"
)

cat("\n[OK] CIV robustness analyses completed.\n")
