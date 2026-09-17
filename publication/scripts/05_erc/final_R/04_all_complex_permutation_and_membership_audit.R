#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

set.seed(20260810)

erc_dir <- paste0(
  "/path/to/data/genomewide_codeml_kuster/",
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

erc_annotation_file <- file.path(
  erc_dir,
  "erc_mtPCG_composite_annotated.tsv"
)

master_annotation_file <- paste0(
  "/path/to/data/genomewide_codeml_kuster/",
  "07_codeml_genomewide/codeml_master_analysis.tsv"
)

n_permutations <- 100000L

# ----------------------------------------------------------
# Read relative-rate matrices
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
    stop("Cannot identify ID column in: ", path)
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

  keep <- rowSums(is.finite(result)) >= 5
  result[keep, , drop = FALSE]
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
# General helper functions
# ----------------------------------------------------------

clean <- function(x) {
  result <- trimws(as.character(x))
  result[
    is.na(result) |
    result == "" |
    result == "NA"
  ] <- NA_character_
  result
}

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
    stop("Cannot calculate composite for empty gene set")
  }

  colMeans(
    rate_matrix[genes, , drop = FALSE],
    na.rm = TRUE
  )
}

identify_id_column <- function(dat) {

  candidates <- c(
    "nuclear_gene",
    "gene_id",
    "stickleback_gene",
    "stickleback_gene_id",
    "Gene stable ID",
    "gene"
  )

  result <- candidates[
    candidates %in% names(dat)
  ][1]

  if (is.na(result)) {
    stop(
      "Cannot identify gene-ID column. Columns: ",
      paste(names(dat), collapse = ", ")
    )
  }

  result
}

# ----------------------------------------------------------
# Read ERC annotations
# ----------------------------------------------------------

annotation <- fread(erc_annotation_file)

annotation_id <- identify_id_column(annotation)
setnames(annotation, annotation_id, "analysis_gene_id")
annotation <- unique(annotation, by = "analysis_gene_id")

for (column_name in c(
  "own_role",
  "own_complex",
  "new_symbol",
  "old_name",
  "core_status"
)) {
  if (!column_name %in% names(annotation)) {
    annotation[, (column_name) := NA_character_]
  }

  annotation[
    ,
    (column_name) := clean(get(column_name))
  ]
}

# ----------------------------------------------------------
# Define nuclear and mitochondrial complex memberships
# ----------------------------------------------------------

complex_names <- c(
  "CI", "CII", "CIII", "CIV", "CV"
)

nuclear_complex_genes <- setNames(
  lapply(complex_names, function(complex_name) {

    annotation[
      own_role == "subunit" &
        own_complex == complex_name &
        analysis_gene_id %in% rownames(nuclear),
      analysis_gene_id
    ]
  }),
  complex_names
)

mt_complex_genes <- list(
  CI = c(
    "ND1", "ND2", "ND3", "ND4",
    "ND4L", "ND5", "ND6"
  ),

  # Complex II has no mitochondrially encoded subunit.
  # Following Weaver, compare it with all 13 mtPCGs.
  CII = rownames(mt),

  CIII = "CYTB",

  CIV = c(
    "COX1", "COX2", "COX3"
  ),

  CV = c(
    "ATP6", "ATP8"
  )
)

for (complex_name in complex_names) {

  mt_complex_genes[[complex_name]] <- intersect(
    mt_complex_genes[[complex_name]],
    rownames(mt)
  )

  cat(
    "[membership]",
    complex_name,
    "nuclear =", length(
      nuclear_complex_genes[[complex_name]]
    ),
    "mt =", length(
      mt_complex_genes[[complex_name]]
    ),
    "\n"
  )
}

# ----------------------------------------------------------
# Classify terminal and internal branches
# ----------------------------------------------------------

branch_is_terminal <- function(branch_label) {

  matches <- regmatches(
    branch_label,
    gregexpr("[0-9]+", branch_label)
  )[[1]]

  numbers <- suppressWarnings(
    as.integer(matches)
  )

  if (length(numbers) < 2) {
    return(NA)
  }

  # The 27 population tips are numbered 1 through 27.
  any(numbers <= 27)
}

terminal_status <- vapply(
  common_branches,
  branch_is_terminal,
  logical(1)
)

if (anyNA(terminal_status)) {
  stop(
    "Unable to classify terminal/internal status for: ",
    paste(
      common_branches[is.na(terminal_status)],
      collapse = ", "
    )
  )
}

terminal_indices <- which(terminal_status)
internal_indices <- which(!terminal_status)

cat(
  "[branches] terminal =", length(terminal_indices),
  "internal =", length(internal_indices),
  "\n"
)

# ----------------------------------------------------------
# Run the same permutation test for every complex
# ----------------------------------------------------------

permutation_results <- list()
permutation_nulls <- list()

for (complex_index in seq_along(complex_names)) {

  complex_name <- complex_names[complex_index]

  nuclear_genes <- nuclear_complex_genes[[complex_name]]

  mitochondrial_genes <- mt_complex_genes[[complex_name]]

  if (
    length(nuclear_genes) == 0 ||
    length(mitochondrial_genes) == 0
  ) {
    warning(
      "Skipping ", complex_name,
      " because one gene set is empty"
    )
    next
  }

  nuclear_rate <- composite(
    nuclear,
    nuclear_genes
  )

  mt_rate <- composite(
    mt,
    mitochondrial_genes
  )

  observed_rs <- safe_cor(
    nuclear_rate,
    mt_rate
  )

  unrestricted_null <- numeric(
    n_permutations
  )

  stratified_null <- numeric(
    n_permutations
  )

  cat(
    "\n[permutation]",
    complex_name,
    "observed rs =", observed_rs,
    "iterations =", n_permutations,
    "\n"
  )

  for (iteration in seq_len(n_permutations)) {

    unrestricted_null[iteration] <- safe_cor(
      nuclear_rate,
      sample(mt_rate)
    )

    permuted_mt <- mt_rate

    permuted_mt[terminal_indices] <- sample(
      mt_rate[terminal_indices]
    )

    permuted_mt[internal_indices] <- sample(
      mt_rate[internal_indices]
    )

    stratified_null[iteration] <- safe_cor(
      nuclear_rate,
      permuted_mt
    )
  }

  unrestricted_p <- (
    1 + sum(
      unrestricted_null >= observed_rs,
      na.rm = TRUE
    )
  ) / (
    1 + sum(is.finite(unrestricted_null))
  )

  stratified_p <- (
    1 + sum(
      stratified_null >= observed_rs,
      na.rm = TRUE
    )
  ) / (
    1 + sum(is.finite(stratified_null))
  )

  permutation_results[[complex_name]] <- data.table(
    complex = complex_name,
    nuclear_n = length(nuclear_genes),
    mt_n = length(mitochondrial_genes),
    nuclear_genes = paste(
      nuclear_genes,
      collapse = ","
    ),
    mt_genes = paste(
      mitochondrial_genes,
      collapse = ","
    ),
    observed_rs = observed_rs,

    unrestricted_null_mean = mean(
      unrestricted_null,
      na.rm = TRUE
    ),

    unrestricted_null_95 = quantile(
      unrestricted_null,
      0.95,
      na.rm = TRUE,
      names = FALSE
    ),

    unrestricted_p = unrestricted_p,

    stratified_null_mean = mean(
      stratified_null,
      na.rm = TRUE
    ),

    stratified_null_95 = quantile(
      stratified_null,
      0.95,
      na.rm = TRUE,
      names = FALSE
    ),

    stratified_p = stratified_p
  )

  permutation_nulls[[complex_name]] <- data.table(
    iteration = seq_len(n_permutations),
    complex = complex_name,
    unrestricted_rs = unrestricted_null,
    stratified_rs = stratified_null
  )
}

permutation_results <- rbindlist(
  permutation_results,
  fill = TRUE
)

permutation_results[
  ,
  unrestricted_padj :=
    p.adjust(
      unrestricted_p,
      method = "BH"
    )
]

permutation_results[
  ,
  stratified_padj :=
    p.adjust(
      stratified_p,
      method = "BH"
    )
]

permutation_results[
  ,
  unrestricted_bonferroni :=
    p.adjust(
      unrestricted_p,
      method = "bonferroni"
    )
]

permutation_results[
  ,
  stratified_bonferroni :=
    p.adjust(
      stratified_p,
      method = "bonferroni"
    )
]

setorder(
  permutation_results,
  stratified_p
)

permutation_nulls <- rbindlist(
  permutation_nulls,
  fill = TRUE
)

cat("\n===== All-complex permutation results =====\n")

print(
  permutation_results[
    ,
    .(
      complex,
      nuclear_n,
      mt_n,
      observed_rs,
      unrestricted_p,
      unrestricted_padj,
      stratified_p,
      stratified_padj,
      stratified_bonferroni
    )
  ]
)

# ----------------------------------------------------------
# Membership audit using the master annotation table
# ----------------------------------------------------------

if (!file.exists(master_annotation_file)) {
  stop(
    "Master annotation file not found: ",
    master_annotation_file
  )
}

master <- fread(master_annotation_file)

master_id <- identify_id_column(master)
setnames(master, master_id, "analysis_gene_id")

master <- unique(
  master,
  by = "analysis_gene_id"
)

for (column_name in c(
  "own_role",
  "own_complex",
  "new_symbol",
  "old_name",
  "core_status"
)) {

  if (!column_name %in% names(master)) {
    master[, (column_name) := NA_character_]
  }

  master[
    ,
    (column_name) := clean(get(column_name))
  ]
}

all_oxphos_structural <- master[
  own_role == "subunit" &
    own_complex %in% complex_names
]

all_oxphos_structural[
  ,
  in_nuclear_rer_matrix :=
    analysis_gene_id %in% rownames(nuclear)
]

all_oxphos_structural[
  ,
  in_erc_annotation :=
    analysis_gene_id %in% annotation$analysis_gene_id
]

all_oxphos_structural[
  ,
  used_in_complex_erc :=
    in_nuclear_rer_matrix &
    in_erc_annotation
]

all_oxphos_structural[
  ,
  exclusion_status := fifelse(
    used_in_complex_erc,
    "included",
    fifelse(
      !in_nuclear_rer_matrix,
      "absent_from_nuclear_RER_matrix",
      "absent_from_ERC_annotation"
    )
  )
]

# If useful filter/status columns exist, retain them in output.
preferred_columns <- c(
  "analysis_gene_id",
  "new_symbol",
  "old_name",
  "own_role",
  "own_complex",
  "core_status",
  "in_nuclear_rer_matrix",
  "in_erc_annotation",
  "used_in_complex_erc",
  "exclusion_status",
  "omega",
  "dN",
  "dS",
  "N",
  "S",
  "quality_status",
  "filter_status",
  "codeml_status"
)

audit_columns <- intersect(
  preferred_columns,
  names(all_oxphos_structural)
)

membership_audit <- all_oxphos_structural[
  ,
  ..audit_columns
]

setorder(
  membership_audit,
  own_complex,
  -used_in_complex_erc,
  new_symbol,
  analysis_gene_id
)

membership_summary <- membership_audit[
  ,
  .(
    curated_total = .N,
    included_in_erc = sum(
      used_in_complex_erc,
      na.rm = TRUE
    ),
    excluded_from_erc = sum(
      !used_in_complex_erc,
      na.rm = TRUE
    )
  ),
  by = own_complex
]

membership_summary[
  ,
  percent_retained :=
    100 * included_in_erc / curated_total
]

cat("\n===== OXPHOS structural membership audit =====\n")
print(membership_summary[order(own_complex)])

cat("\n===== CIV membership audit =====\n")

print(
  membership_audit[
    own_complex == "CIV"
  ]
)

cat("\n===== Included CIV genes =====\n")

print(
  membership_audit[
    own_complex == "CIV" &
      used_in_complex_erc == TRUE
  ]
)

cat("\n===== Excluded CIV genes =====\n")

print(
  membership_audit[
    own_complex == "CIV" &
      used_in_complex_erc == FALSE
  ]
)

# ----------------------------------------------------------
# Investigate unnamed CIV locus
# ----------------------------------------------------------

unnamed_civ <- membership_audit[
  own_complex == "CIV" &
    (
      is.na(new_symbol) |
      new_symbol == analysis_gene_id |
      grepl("^ENSGACG", new_symbol)
    )
]

cat("\n===== Unnamed CIV loci =====\n")
print(unnamed_civ)

# Output every master-table column for unnamed loci so that
# possible annotation fields are not accidentally omitted.
unnamed_civ_full <- master[
  analysis_gene_id %in%
    unnamed_civ$analysis_gene_id
]

# ----------------------------------------------------------
# Save outputs
# ----------------------------------------------------------

fwrite(
  permutation_results,
  file.path(
    erc_dir,
    "OXPHOS_all_complex_permutation_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  permutation_nulls,
  file.path(
    erc_dir,
    "OXPHOS_all_complex_permutation_null_100000.tsv.gz"
  ),
  sep = "\t"
)

fwrite(
  membership_summary,
  file.path(
    erc_dir,
    "OXPHOS_complex_membership_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  membership_audit,
  file.path(
    erc_dir,
    "OXPHOS_structural_membership_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  membership_audit[
    own_complex == "CIV"
  ],
  file.path(
    erc_dir,
    "CIV_structural_membership_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  unnamed_civ_full,
  file.path(
    erc_dir,
    "CIV_unnamed_loci_full_annotation.tsv"
  ),
  sep = "\t"
)

cat(
  "\n[OK] All-complex permutation and ",
  "membership audit completed.\n"
)
