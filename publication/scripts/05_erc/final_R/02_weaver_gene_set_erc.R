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

annotation_file <- file.path(
  erc_dir,
  "erc_mtPCG_composite_annotated.tsv"
)

n_boot <- 10000L

# ----------------------------------------------------------
# Read relative-rate matrices
# ----------------------------------------------------------

read_rate_matrix <- function(path, object_name) {

  dat <- fread(path)

  if (ncol(dat) < 3) {
    stop(object_name, " matrix has fewer than three columns")
  }

  id_candidates <- c(
    "nuclear_gene",
    "gene_id",
    "gene",
    "mt_gene",
    "Gene",
    names(dat)[1]
  )

  id_column <- id_candidates[id_candidates %in% names(dat)][1]

  if (is.na(id_column)) {
    stop("Cannot identify ID column in ", path)
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

  matrix_data <- as.matrix(dat[, ..rate_columns])
  storage.mode(matrix_data) <- "double"

  rownames(matrix_data) <- ids

  keep_rows <- rowSums(is.finite(matrix_data)) > 1
  matrix_data <- matrix_data[keep_rows, , drop = FALSE]

  cat(
    "[load]", object_name,
    "genes =", nrow(matrix_data),
    "branches =", ncol(matrix_data), "\n"
  )

  matrix_data
}

nuclear <- read_rate_matrix(
  nuclear_file,
  "nuclear"
)

mt <- read_rate_matrix(
  mt_file,
  "mitochondrial"
)

common_branches <- intersect(
  colnames(nuclear),
  colnames(mt)
)

if (length(common_branches) < 5) {
  stop(
    "Too few shared branch columns: ",
    paste(common_branches, collapse = ", ")
  )
}

nuclear <- nuclear[, common_branches, drop = FALSE]
mt <- mt[, common_branches, drop = FALSE]

cat("[branches] Shared branches =", length(common_branches), "\n")

# Normalize mitochondrial gene labels
rownames(mt) <- toupper(
  sub("^MT[-_]", "", rownames(mt))
)

# ----------------------------------------------------------
# Read annotations
# ----------------------------------------------------------

annotation <- fread(annotation_file)

annotation_id_candidates <- c(
  "nuclear_gene",
  "gene_id",
  "stickleback_gene",
  "stickleback_gene_id"
)

annotation_id <- annotation_id_candidates[
  annotation_id_candidates %in% names(annotation)
][1]

if (is.na(annotation_id)) {
  stop("Cannot identify annotation gene-ID column")
}

setnames(annotation, annotation_id, "analysis_gene_id")

annotation <- unique(
  annotation,
  by = "analysis_gene_id"
)

annotation <- annotation[
  analysis_gene_id %in% rownames(nuclear)
]

clean_value <- function(values) {
  values <- trimws(as.character(values))
  values[
    is.na(values) |
    values == "" |
    values == "NA"
  ] <- NA_character_
  values
}

for (column_name in c(
  "primary_class",
  "own_role",
  "own_complex",
  "core_status"
)) {
  if (!column_name %in% names(annotation)) {
    annotation[, (column_name) := NA_character_]
  }

  annotation[
    ,
    (column_name) := clean_value(get(column_name))
  ]
}

# Harmonize Kuster labels
annotation[
  primary_class %in% c(
    "direct-n-mt",
    "direct_n-mt",
    "direct n-mt"
  ),
  primary_class := "direct_n-mt"
]

annotation[
  primary_class %in% c(
    "indirect-n-mt",
    "indirect_n-mt",
    "indirect n-mt"
  ),
  primary_class := "indirect_n-mt"
]

annotation[
  primary_class %in% c(
    "non_n_mt",
    "non_n-mt",
    "non-n-mt",
    "non n-mt"
  ),
  primary_class := "non-n-mt"
]

# ----------------------------------------------------------
# Define mitochondrial complexes
# ----------------------------------------------------------

mt_complexes <- list(
  CI = c(
    "ND1", "ND2", "ND3", "ND4",
    "ND4L", "ND5", "ND6"
  ),
  CIII = "CYTB",
  CIV = c("COX1", "COX2", "COX3"),
  CV = c("ATP6", "ATP8")
)

missing_mt <- setdiff(
  unique(unlist(mt_complexes)),
  rownames(mt)
)

if (length(missing_mt) > 0) {
  stop(
    "Missing mitochondrial genes: ",
    paste(missing_mt, collapse = ", ")
  )
}

# ----------------------------------------------------------
# Define nuclear gene sets
# ----------------------------------------------------------

genes_where <- function(mask) {

  mask_expression <- substitute(mask)

  selected <- eval(
    mask_expression,
    envir = as.list(annotation),
    enclos = parent.frame()
  )

  selected[is.na(selected)] <- FALSE

  intersect(
    annotation$analysis_gene_id[selected],
    rownames(nuclear)
  )
}

gene_sets <- list(
  Nmt_OXPHOS_structural = genes_where(
    own_role == "subunit"
  ),

  OXPHOS_assembly_factors = genes_where(
    own_role == "assembly_factor"
  ),

  Nmt_ribosomal = genes_where(
    own_role == "Nmt-ribo"
  ),

  Nmt_ARS = genes_where(
    own_role == "Nmt-ARS"
  ),

  Cyto_ribosomal_control = genes_where(
    own_role == "cyto-ribo"
  ),

  Cyto_ARS_control = genes_where(
    own_role == "cyto-ARS"
  ),

  Direct_nmt = genes_where(
    primary_class == "direct_n-mt"
  ),

  Indirect_nmt = genes_where(
    primary_class == "indirect_n-mt"
  ),

  OXPHOS_core = genes_where(
    core_status == "nu_core"
  ),

  OXPHOS_noncore = genes_where(
    core_status == "nu_noncore"
  )
)

# Complex-specific nuclear structural subunits
for (complex_name in c("CI", "CII", "CIII", "CIV", "CV")) {

  gene_sets[[paste0("OXPHOS_", complex_name)]] <- genes_where(
    own_role == "subunit" &
      own_complex == complex_name
  )
}

cat("\n===== Gene-set sizes =====\n")

for (set_name in names(gene_sets)) {
  cat(set_name, ":", length(gene_sets[[set_name]]), "\n")
}

# ----------------------------------------------------------
# Composite and bootstrap functions
# ----------------------------------------------------------

composite_rate <- function(rate_matrix, genes) {

  genes <- intersect(genes, rownames(rate_matrix))

  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(rate_matrix)))
  }

  colMeans(
    rate_matrix[genes, , drop = FALSE],
    na.rm = TRUE
  )
}

safe_spearman <- function(x, y) {

  keep <- is.finite(x) & is.finite(y)

  if (sum(keep) < 5) {
    return(NA_real_)
  }

  if (
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

bootstrap_erc <- function(
  nuclear_genes,
  mt_genes,
  iterations = 10000L
) {

  nuclear_genes <- intersect(
    nuclear_genes,
    rownames(nuclear)
  )

  mt_genes <- intersect(
    mt_genes,
    rownames(mt)
  )

  if (
    length(nuclear_genes) == 0 ||
    length(mt_genes) == 0
  ) {
    return(rep(NA_real_, iterations))
  }

  result <- numeric(iterations)

  for (iteration in seq_len(iterations)) {

    sampled_nuclear <- sample(
      nuclear_genes,
      size = length(nuclear_genes),
      replace = TRUE
    )

    sampled_mt <- sample(
      mt_genes,
      size = length(mt_genes),
      replace = TRUE
    )

    nuclear_composite <- composite_rate(
      nuclear,
      sampled_nuclear
    )

    mt_composite <- composite_rate(
      mt,
      sampled_mt
    )

    result[iteration] <- safe_spearman(
      nuclear_composite,
      mt_composite
    )
  }

  result
}

summarize_bootstrap <- function(
  analysis,
  nuclear_set,
  mt_set,
  observed,
  bootstrap_values
) {

  finite_bootstrap <- bootstrap_values[
    is.finite(bootstrap_values)
  ]

  if (length(finite_bootstrap) == 0) {
    return(
      data.table(
        analysis = analysis,
        nuclear_set = nuclear_set,
        mt_set = mt_set,
        nuclear_n = length(gene_sets[[nuclear_set]]),
        mt_n = length(mt_set),
        observed_rs = observed,
        bootstrap_mean_rs = NA_real_,
        bootstrap_median_rs = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        bootstrap_p_positive = NA_real_
      )
    )
  }

  data.table(
    analysis = analysis,
    nuclear_set = nuclear_set,
    mt_set = paste(mt_set, collapse = ","),
    nuclear_n = length(gene_sets[[nuclear_set]]),
    mt_n = length(mt_set),
    observed_rs = observed,
    bootstrap_mean_rs = mean(finite_bootstrap),
    bootstrap_median_rs = median(finite_bootstrap),
    ci_lower = quantile(
      finite_bootstrap,
      0.025,
      names = FALSE
    ),
    ci_upper = quantile(
      finite_bootstrap,
      0.975,
      names = FALSE
    ),
    bootstrap_p_positive =
      (1 + sum(finite_bootstrap <= 0)) /
      (length(finite_bootstrap) + 1)
  )
}

# ----------------------------------------------------------
# Run Weaver-style analyses
# ----------------------------------------------------------

analysis_definitions <- list(
  list(
    analysis = "functional_set",
    nuclear_set = "Nmt_OXPHOS_structural",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "functional_set",
    nuclear_set = "OXPHOS_assembly_factors",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "functional_set",
    nuclear_set = "Nmt_ribosomal",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "functional_set",
    nuclear_set = "Nmt_ARS",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "functional_set",
    nuclear_set = "Cyto_ribosomal_control",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "functional_set",
    nuclear_set = "Cyto_ARS_control",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "Kuster_class",
    nuclear_set = "Direct_nmt",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "Kuster_class",
    nuclear_set = "Indirect_nmt",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "core_status",
    nuclear_set = "OXPHOS_core",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "core_status",
    nuclear_set = "OXPHOS_noncore",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "within_complex",
    nuclear_set = "OXPHOS_CI",
    mt_genes = mt_complexes$CI
  ),
  list(
    analysis = "within_complex_control",
    nuclear_set = "OXPHOS_CII",
    mt_genes = rownames(mt)
  ),
  list(
    analysis = "within_complex",
    nuclear_set = "OXPHOS_CIII",
    mt_genes = mt_complexes$CIII
  ),
  list(
    analysis = "within_complex",
    nuclear_set = "OXPHOS_CIV",
    mt_genes = mt_complexes$CIV
  ),
  list(
    analysis = "within_complex",
    nuclear_set = "OXPHOS_CV",
    mt_genes = mt_complexes$CV
  )
)

summary_results <- list()
bootstrap_table <- list()

for (index in seq_along(analysis_definitions)) {

  definition <- analysis_definitions[[index]]

  set_name <- definition$nuclear_set
  nuclear_genes <- gene_sets[[set_name]]
  mitochondrial_genes <- definition$mt_genes

  cat(
    "\n[bootstrap]",
    set_name,
    "nuclear =", length(nuclear_genes),
    "mt =", length(mitochondrial_genes),
    "\n"
  )

  observed_nuclear <- composite_rate(
    nuclear,
    nuclear_genes
  )

  observed_mt <- composite_rate(
    mt,
    mitochondrial_genes
  )

  observed_rs <- safe_spearman(
    observed_nuclear,
    observed_mt
  )

  bootstrap_values <- bootstrap_erc(
    nuclear_genes,
    mitochondrial_genes,
    n_boot
  )

  summary_results[[index]] <- summarize_bootstrap(
    definition$analysis,
    set_name,
    mitochondrial_genes,
    observed_rs,
    bootstrap_values
  )

  bootstrap_table[[index]] <- data.table(
    iteration = seq_len(n_boot),
    analysis = definition$analysis,
    nuclear_set = set_name,
    rs = bootstrap_values
  )
}

summary_results <- rbindlist(
  summary_results,
  fill = TRUE
)

summary_results[
  ,
  bootstrap_padj_positive :=
    p.adjust(
      bootstrap_p_positive,
      method = "BH"
    )
]

bootstrap_table <- rbindlist(
  bootstrap_table,
  fill = TRUE
)

# ----------------------------------------------------------
# Compare each focal distribution against cyto-ribo control
# ----------------------------------------------------------

control_values <- bootstrap_table[
  nuclear_set == "Cyto_ribosomal_control",
  rs
]

comparison_results <- list()
comparison_counter <- 0L

focal_sets <- c(
  "Nmt_OXPHOS_structural",
  "OXPHOS_assembly_factors",
  "Nmt_ribosomal",
  "Nmt_ARS",
  "Direct_nmt",
  "Indirect_nmt",
  "OXPHOS_core",
  "OXPHOS_noncore"
)

for (set_name in focal_sets) {

  focal_values <- bootstrap_table[
    nuclear_set == set_name,
    rs
  ]

  comparison_length <- min(
    length(focal_values),
    length(control_values)
  )

  difference <- focal_values[
    seq_len(comparison_length)
  ] - control_values[
    seq_len(comparison_length)
  ]

  difference <- difference[
    is.finite(difference)
  ]

  if (length(difference) == 0) {
    next
  }

  comparison_counter <- comparison_counter + 1L

  comparison_results[[comparison_counter]] <- data.table(
    focal_set = set_name,
    control_set = "Cyto_ribosomal_control",
    mean_difference_rs = mean(difference),
    median_difference_rs = median(difference),
    difference_ci_lower = quantile(
      difference,
      0.025,
      names = FALSE
    ),
    difference_ci_upper = quantile(
      difference,
      0.975,
      names = FALSE
    ),
    bootstrap_p_focal_greater =
      (1 + sum(difference <= 0)) /
      (length(difference) + 1)
  )
}

comparison_results <- rbindlist(
  comparison_results,
  fill = TRUE
)

comparison_results[
  ,
  bootstrap_padj_focal_greater :=
    p.adjust(
      bootstrap_p_focal_greater,
      method = "BH"
    )
]

# ----------------------------------------------------------
# Save branch-level composite values
# ----------------------------------------------------------

composite_rows <- list()
composite_counter <- 0L

for (definition in analysis_definitions) {

  set_name <- definition$nuclear_set

  nuclear_composite <- composite_rate(
    nuclear,
    gene_sets[[set_name]]
  )

  mt_composite <- composite_rate(
    mt,
    definition$mt_genes
  )

  composite_counter <- composite_counter + 1L

  composite_rows[[composite_counter]] <- data.table(
    analysis = definition$analysis,
    nuclear_set = set_name,
    branch = common_branches,
    nuclear_composite_rer = nuclear_composite,
    mt_composite_rer = mt_composite
  )
}

composite_rows <- rbindlist(
  composite_rows,
  fill = TRUE
)

# ----------------------------------------------------------
# Output
# ----------------------------------------------------------

fwrite(
  summary_results,
  file.path(
    erc_dir,
    "weaver_gene_set_erc_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  bootstrap_table,
  file.path(
    erc_dir,
    "weaver_gene_set_erc_bootstrap_10000.tsv.gz"
  ),
  sep = "\t"
)

fwrite(
  comparison_results,
  file.path(
    erc_dir,
    "weaver_gene_set_erc_vs_cytoribo.tsv"
  ),
  sep = "\t"
)

fwrite(
  composite_rows,
  file.path(
    erc_dir,
    "weaver_gene_set_branch_composites.tsv"
  ),
  sep = "\t"
)

cat("\n\n===== Weaver-style gene-set ERC =====\n")

print(
  summary_results[
    ,
    .(
      analysis,
      nuclear_set,
      nuclear_n,
      mt_n,
      observed_rs,
      bootstrap_mean_rs,
      ci_lower,
      ci_upper,
      bootstrap_p_positive,
      bootstrap_padj_positive
    )
  ]
)

cat("\n===== Difference from cyto-ribosomal control =====\n")

print(comparison_results)

cat("\n[OK] Weaver-style ERC analysis completed.\n")
