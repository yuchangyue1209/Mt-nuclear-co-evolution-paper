#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

erc_dir <- "/path/to/data/genomewide_codeml_kuster/10_erc_results/mtPCG_composite_t_spearman"

input_file <- file.path(
  erc_dir,
  "erc_mtPCG_composite_annotated.tsv"
)

x <- fread(input_file)

required <- c(
  "nuclear_gene", "r", "p", "padj", "primary_class",
  "own_role", "own_complex", "core_status"
)

missing_columns <- setdiff(required, names(x))

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

x[, r := as.numeric(r)]
x <- x[is.finite(r)]

# Harmonize class labels
x[, primary_class_clean := gsub("_", "-", primary_class)]
x[
  primary_class_clean %in% c("direct-n-mt", "direct n-mt"),
  primary_class_clean := "direct_n-mt"
]
x[
  primary_class_clean %in% c("indirect-n-mt", "indirect n-mt"),
  primary_class_clean := "indirect_n-mt"
]
x[
  primary_class_clean %in% c("non-n-mt", "non n-mt"),
  primary_class_clean := "non-n-mt"
]

background <- x[
  primary_class_clean == "non-n-mt" & is.finite(r),
  r
]

if (length(background) == 0) {
  stop("No finite non-n-mt background genes found")
}

cat("Non-n-mt background N:", length(background), "\n")
cat("Background median:", median(background), "\n\n")

# ----------------------------------------------------------
# Empirical percentile and empirical upper-tail P per gene
# ----------------------------------------------------------

sorted_background <- sort(background)
n_background <- length(sorted_background)

empirical_percentile <- function(value) {
  findInterval(value, sorted_background) / n_background
}

empirical_upper_p <- function(value) {
  (1 + sum(sorted_background >= value)) /
    (n_background + 1)
}

x[, background_percentile := vapply(
  r,
  empirical_percentile,
  numeric(1)
)]

x[, empirical_p_upper := vapply(
  r,
  empirical_upper_p,
  numeric(1)
)]

x[, empirical_padj := p.adjust(
  empirical_p_upper,
  method = "BH"
)]

# Background thresholds
proportions <- c(
  top10 = 0.90,
  top05 = 0.95,
  top025 = 0.975,
  top01 = 0.99
)

thresholds <- quantile(
  background,
  probs = proportions,
  na.rm = TRUE,
  names = TRUE
)

for (threshold_name in names(thresholds)) {
  x[
    ,
    (paste0("above_", threshold_name)) :=
      r >= thresholds[[threshold_name]]
  ]
}

threshold_table <- data.table(
  threshold = names(thresholds),
  background_quantile = as.numeric(proportions),
  r_cutoff = as.numeric(thresholds)
)

fwrite(
  threshold_table,
  file.path(erc_dir, "erc_background_thresholds.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------
# Construct category membership
# ----------------------------------------------------------

membership <- list()

add_category <- function(group, category, mask) {
  membership[[length(membership) + 1]] <<- data.table(
    nuclear_gene = x$nuclear_gene[mask],
    group = group,
    category = category
  )
}

# Kuster categories
for (category_name in c(
  "direct_n-mt",
  "indirect_n-mt",
  "non-n-mt"
)) {
  add_category(
    "Kuster_class",
    category_name,
    x$primary_class_clean == category_name
  )
}

# Functional roles
roles <- sort(unique(na.omit(x$own_role)))

for (category_name in roles) {
  if (
    nzchar(category_name) &&
    category_name != "NA"
  ) {
    add_category(
      "functional_role",
      category_name,
      !is.na(x$own_role) &
        x$own_role == category_name
    )
  }
}

# OXPHOS complexes
complexes <- sort(unique(na.omit(x$own_complex)))

for (category_name in complexes) {
  if (
    nzchar(category_name) &&
    category_name != "NA"
  ) {
    add_category(
      "OXPHOS_complex",
      category_name,
      !is.na(x$own_complex) &
        x$own_complex == category_name
    )
  }
}

# Core/noncore
core_categories <- sort(unique(na.omit(x$core_status)))

for (category_name in core_categories) {
  if (
    nzchar(category_name) &&
    category_name != "NA"
  ) {
    add_category(
      "core_status",
      category_name,
      !is.na(x$core_status) &
        x$core_status == category_name
    )
  }
}

membership_table <- unique(rbindlist(membership, fill = TRUE))

# ----------------------------------------------------------
# Distribution tests: Wilcoxon + rank-biserial/Cliff's delta
# ----------------------------------------------------------

distribution_results <- list()
counter <- 0L

category_rows <- unique(
  membership_table[, .(group, category)]
)

for (row_index in seq_len(nrow(category_rows))) {

  current_group <- category_rows$group[row_index]
  current_category <- category_rows$category[row_index]

  # Do not compare the background category against itself
  if (
    current_group == "Kuster_class" &&
    current_category == "non-n-mt"
  ) {
    next
  }

  genes_in_category <- membership_table[
    group == current_group &
      category == current_category,
    nuclear_gene
  ]

  focal <- x[
    nuclear_gene %in% genes_in_category,
    r
  ]

  focal <- focal[is.finite(focal)]

  if (length(focal) < 2) {
    next
  }

  test <- suppressWarnings(
    wilcox.test(
      focal,
      background,
      alternative = "greater",
      exact = FALSE
    )
  )

  # For a two-sample Wilcoxon test, W is the Mann-Whitney U.
  auc <- as.numeric(test$statistic) /
    (length(focal) * length(background))

  cliffs_delta <- 2 * auc - 1

  counter <- counter + 1L

  distribution_results[[counter]] <- data.table(
    group = current_group,
    category = current_category,
    focal_n = length(focal),
    background_n = length(background),
    focal_median_r = median(focal),
    background_median_r = median(background),
    median_difference = median(focal) - median(background),
    auc_probability = auc,
    cliffs_delta = cliffs_delta,
    wilcoxon_p_greater = test$p.value
  )
}

distribution_results <- rbindlist(
  distribution_results,
  fill = TRUE
)

distribution_results[
  ,
  wilcoxon_padj :=
    p.adjust(wilcoxon_p_greater, method = "BH")
]

setorder(
  distribution_results,
  wilcoxon_padj,
  -cliffs_delta
)

fwrite(
  distribution_results,
  file.path(erc_dir, "erc_category_distribution_tests.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------
# Fisher enrichment at top 10%, 5%, 2.5%, and 1%
# ----------------------------------------------------------

enrichment_results <- list()
counter <- 0L

for (row_index in seq_len(nrow(category_rows))) {

  current_group <- category_rows$group[row_index]
  current_category <- category_rows$category[row_index]

  if (
    current_group == "Kuster_class" &&
    current_category == "non-n-mt"
  ) {
    next
  }

  genes_in_category <- membership_table[
    group == current_group &
      category == current_category,
    nuclear_gene
  ]

  focal_rows <- x[nuclear_gene %in% genes_in_category]
  background_rows <- x[primary_class_clean == "non-n-mt"]

  if (nrow(focal_rows) == 0) {
    next
  }

  for (threshold_name in names(thresholds)) {

    cutoff <- thresholds[[threshold_name]]

    focal_above <- sum(focal_rows$r >= cutoff)
    focal_below <- sum(focal_rows$r < cutoff)

    background_above <- sum(background_rows$r >= cutoff)
    background_below <- sum(background_rows$r < cutoff)

    contingency <- matrix(
      c(
        focal_above,
        focal_below,
        background_above,
        background_below
      ),
      nrow = 2,
      byrow = TRUE
    )

    fisher_result <- fisher.test(
      contingency,
      alternative = "greater"
    )

    counter <- counter + 1L

    enrichment_results[[counter]] <- data.table(
      group = current_group,
      category = current_category,
      threshold = threshold_name,
      r_cutoff = cutoff,
      focal_total = nrow(focal_rows),
      focal_above = focal_above,
      focal_percent_above =
        100 * focal_above / nrow(focal_rows),
      background_total = nrow(background_rows),
      background_above = background_above,
      background_percent_above =
        100 * background_above / nrow(background_rows),
      odds_ratio = unname(fisher_result$estimate),
      fisher_p_greater = fisher_result$p.value
    )
  }
}

enrichment_results <- rbindlist(
  enrichment_results,
  fill = TRUE
)

# FDR across all category × threshold tests
enrichment_results[
  ,
  fisher_padj_global :=
    p.adjust(fisher_p_greater, method = "BH")
]

# Also correct within each threshold
enrichment_results[
  ,
  fisher_padj_within_threshold :=
    p.adjust(fisher_p_greater, method = "BH"),
  by = threshold
]

setorder(
  enrichment_results,
  threshold,
  fisher_padj_within_threshold,
  -odds_ratio
)

fwrite(
  enrichment_results,
  file.path(erc_dir, "erc_category_threshold_enrichment.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------
# Candidate tables
# ----------------------------------------------------------

setorder(x, -r)

fwrite(
  x,
  file.path(erc_dir, "erc_all_genes_empirical_statistics.tsv.gz"),
  sep = "\t"
)

candidate_top05 <- x[
  r >= thresholds[["top05"]]
]

candidate_top01 <- x[
  r >= thresholds[["top01"]]
]

fwrite(
  candidate_top05,
  file.path(erc_dir, "erc_candidates_background_top05.tsv"),
  sep = "\t"
)

fwrite(
  candidate_top01,
  file.path(erc_dir, "erc_candidates_background_top01.tsv"),
  sep = "\t"
)

mitochondrial_top05 <- candidate_top05[
  primary_class_clean %in%
    c("direct_n-mt", "indirect_n-mt") |
  (!is.na(own_role) & nzchar(own_role))
]

mitochondrial_top01 <- candidate_top01[
  primary_class_clean %in%
    c("direct_n-mt", "indirect_n-mt") |
  (!is.na(own_role) & nzchar(own_role))
]

fwrite(
  mitochondrial_top05,
  file.path(
    erc_dir,
    "erc_mitochondrial_candidates_background_top05.tsv"
  ),
  sep = "\t"
)

fwrite(
  mitochondrial_top01,
  file.path(
    erc_dir,
    "erc_mitochondrial_candidates_background_top01.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------
# Console summary
# ----------------------------------------------------------

cat("\n===== Background thresholds =====\n")
print(threshold_table)

cat("\n===== Distribution tests =====\n")
print(
  distribution_results[
    ,
    .(
      group,
      category,
      focal_n,
      focal_median_r,
      median_difference,
      cliffs_delta,
      wilcoxon_p_greater,
      wilcoxon_padj
    )
  ]
)

cat("\n===== Threshold enrichment =====\n")
print(
  enrichment_results[
    ,
    .(
      group,
      category,
      threshold,
      focal_total,
      focal_above,
      focal_percent_above,
      odds_ratio,
      fisher_p_greater,
      fisher_padj_within_threshold
    )
  ]
)

cat("\n===== Mitochondrial-associated top 5% candidates =====\n")
cat("N =", nrow(mitochondrial_top05), "\n")

print(
  mitochondrial_top05[
    ,
    .(
      nuclear_gene,
      new_symbol,
      old_name,
      r,
      background_percentile,
      empirical_p_upper,
      padj,
      primary_class,
      own_role,
      own_complex,
      core_status
    )
  ]
)

cat("\n===== Mitochondrial-associated top 1% candidates =====\n")
cat("N =", nrow(mitochondrial_top01), "\n")

print(
  mitochondrial_top01[
    ,
    .(
      nuclear_gene,
      new_symbol,
      old_name,
      r,
      background_percentile,
      empirical_p_upper,
      padj,
      primary_class,
      own_role,
      own_complex,
      core_status
    )
  ]
)

cat("\n[OK] ERC statistical screening completed.\n")
