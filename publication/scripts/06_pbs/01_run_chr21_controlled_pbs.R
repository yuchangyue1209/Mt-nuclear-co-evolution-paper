#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# Configuration
# ============================================================

set.seed(123)

N_PERM <- 10000
MIN_POP <- 10

BASE_DIR <- paste0(
  "/path/to/workspace/gene_fst_work_withNorway/",
  "pbs_mt_vs_nuclear_NorwayOutgroup"
)

# Try these files in order
GENE_PBS_CANDIDATES <- c(
  file.path(
    BASE_DIR,
    "OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_noAMO.tsv"
  ),
  file.path(
    BASE_DIR,
    "OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_by_gene_population.tsv"
  )
)

CHR21_PBS_FILE <- file.path(
  BASE_DIR,
  "chr21_neutral_PBS_by_population.tsv"
)

OUTDIR <- file.path(
  BASE_DIR,
  "OXPHOS72_chr21_controlled"
)

OUT_TABLE <- file.path(
  OUTDIR,
  "OXPHOS72_chr21_controlled_mitonuclear_PBS.tsv"
)

OUT_CANDIDATES <- file.path(
  OUTDIR,
  "OXPHOS72_chr21_controlled_candidates.tsv"
)

OUT_PNG <- file.path(
  OUTDIR,
  "OXPHOS72_chr21_controlled_mitonuclear_PBS.png"
)

OUT_PDF <- file.path(
  OUTDIR,
  "OXPHOS72_chr21_controlled_mitonuclear_PBS.pdf"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Locate input file
# ============================================================

existing_gene_files <- GENE_PBS_CANDIDATES[
  file.exists(GENE_PBS_CANDIDATES)
]

if (length(existing_gene_files) == 0) {
  stop(
    "Could not find an OXPHOS72 PBS input file.\nChecked:\n",
    paste(GENE_PBS_CANDIDATES, collapse = "\n")
  )
}

GENE_PBS_FILE <- existing_gene_files[1]

if (!file.exists(CHR21_PBS_FILE)) {
  stop(
    "chr21 PBS file does not exist:\n",
    CHR21_PBS_FILE
  )
}

cat("[Input] Gene PBS:\n", GENE_PBS_FILE, "\n\n")
cat("[Input] chr21 PBS:\n", CHR21_PBS_FILE, "\n\n")

# ============================================================
# Read data
# ============================================================

gene_dt <- fread(GENE_PBS_FILE)
chr21_dt <- fread(CHR21_PBS_FILE)

cat("[Info] Gene PBS columns:\n")
print(names(gene_dt))

cat("\n[Info] chr21 PBS columns:\n")
print(names(chr21_dt))

required_gene_columns <- c(
  "gene",
  "focal",
  "region",
  "nu_PBS",
  "mt_PBS"
)

required_chr21_columns <- c(
  "focal",
  "region",
  "chr21_PBS"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  names(gene_dt)
)

missing_chr21_columns <- setdiff(
  required_chr21_columns,
  names(chr21_dt)
)

if (length(missing_gene_columns) > 0) {
  stop(
    "Missing columns in gene PBS table: ",
    paste(missing_gene_columns, collapse = ", ")
  )
}

if (length(missing_chr21_columns) > 0) {
  stop(
    "Missing columns in chr21 PBS table: ",
    paste(missing_chr21_columns, collapse = ", ")
  )
}

# ============================================================
# Prepare and merge data
# ============================================================

# Match the main PBS analysis
gene_dt <- gene_dt[
  !focal %in% c("AMO", "LB")
]

chr21_dt <- chr21_dt[
  !focal %in% c("AMO", "LB")
]

chr21_use <- unique(
  chr21_dt[
    ,
    .(
      focal,
      region,
      chr21_PBS
    )
  ]
)

dat <- merge(
  gene_dt,
  chr21_use,
  by = c("focal", "region"),
  all.x = TRUE
)

dat <- dat[
  is.finite(nu_PBS) &
  is.finite(mt_PBS) &
  is.finite(chr21_PBS)
]

dat[, region := factor(region)]

# Check duplicated gene-population records
duplicate_rows <- dat[
  ,
  .N,
  by = .(
    gene,
    focal,
    region
  )
][N > 1]

if (nrow(duplicate_rows) > 0) {
  cat(
    "\n[Warning] Duplicate gene-population rows detected.\n",
    "Averaging PBS values within duplicated records.\n"
  )

  dat <- dat[
    ,
    .(
      nu_PBS = mean(nu_PBS, na.rm = TRUE),
      mt_PBS = mean(mt_PBS, na.rm = TRUE),
      chr21_PBS = mean(chr21_PBS, na.rm = TRUE)
    ),
    by = .(
      gene,
      focal,
      region
    )
  ]

  dat[, region := factor(region)]
}

# ============================================================
# Dataset checks
# ============================================================

population_table <- unique(
  dat[
    ,
    .(
      focal,
      region
    )
  ]
)

setorder(
  population_table,
  region,
  focal
)

cat("\n[Info] Populations retained:\n")
print(population_table)

cat("\n[Info] Population counts by region:\n")

region_counts <- population_table[
  ,
  .N,
  by = region
]

print(region_counts)

if (
  region_counts[region == "AK", N] != 8 ||
  region_counts[region == "BC", N] != 11
) {
  warning(
    paste0(
      "Expected 8 AK and 11 BC populations, but observed:\n",
      paste(
        region_counts$region,
        region_counts$N,
        collapse = "; "
      )
    )
  )
}

cat(
  "\n[Info] Number of genes:",
  uniqueN(dat$gene),
  "\n"
)

# ============================================================
# Overall covariate diagnostics
# ============================================================

covariate_dt <- unique(
  dat[
    ,
    .(
      focal,
      region,
      mt_PBS,
      chr21_PBS
    )
  ]
)

raw_mt_chr21_cor <- cor(
  covariate_dt$mt_PBS,
  covariate_dt$chr21_PBS,
  use = "complete.obs"
)

mt_covariate_model <- lm(
  mt_PBS ~ chr21_PBS + region,
  data = covariate_dt
)

mt_covariate_R2 <- summary(
  mt_covariate_model
)$r.squared

mt_vif <- 1 / (1 - mt_covariate_R2)

cat(
  "\n[Diagnostic] Raw correlation between mtPBS and chr21 PBS:",
  sprintf("%.4f", raw_mt_chr21_cor),
  "\n"
)

cat(
  "[Diagnostic] R2 for mtPBS ~ chr21PBS + region:",
  sprintf("%.4f", mt_covariate_R2),
  "\n"
)

cat(
  "[Diagnostic] Approximate VIF for mtPBS:",
  sprintf("%.4f", mt_vif),
  "\n"
)

if (mt_vif > 5) {
  warning(
    paste0(
      "mtPBS is strongly associated with chr21PBS and/or region. ",
      "Interpret individual coefficients cautiously."
    )
  )
}

# ============================================================
# Helper: permute residuals within region
# ============================================================

permute_within_region <- function(
  residual_vector,
  region_vector
) {

  output <- residual_vector

  for (reg in unique(region_vector)) {

    index <- which(
      region_vector == reg
    )

    output[index] <- sample(
      residual_vector[index],
      length(index),
      replace = FALSE
    )
  }

  output
}

# ============================================================
# Per-gene Freedman-Lane test
# ============================================================

test_one_gene <- function(
  d,
  nperm = 10000
) {

  d <- copy(d)

  d <- d[
    is.finite(nu_PBS) &
    is.finite(mt_PBS) &
    is.finite(chr21_PBS)
  ]

  d[, region := droplevels(factor(region))]

  n_population <- nrow(d)

  empty_result <- function() {
    data.table(
      n_pop = n_population,
      n_AK = sum(d$region == "AK"),
      n_BC = sum(d$region == "BC"),
      beta_mt = NA_real_,
      beta_chr21 = NA_real_,
      partial_R2_mt = NA_real_,
      cor_raw = NA_real_,
      cor_residual = NA_real_,
      p_parametric = NA_real_,
      p_perm = NA_real_,
      mt_vif = NA_real_
    )
  }

  if (n_population < MIN_POP) {
    return(empty_result())
  }

  if (
    length(unique(d$region)) < 2 ||
    var(d$nu_PBS) <= 0 ||
    var(d$mt_PBS) <= 0 ||
    var(d$chr21_PBS) <= 0
  ) {
    return(empty_result())
  }

  # Full and reduced design matrices
  X_full <- model.matrix(
    ~ mt_PBS + chr21_PBS + region,
    data = d
  )

  X_reduced <- model.matrix(
    ~ chr21_PBS + region,
    data = d
  )

  y <- d$nu_PBS

  # Ensure matrices have full rank
  if (
    qr(X_full)$rank < ncol(X_full) ||
    qr(X_reduced)$rank < ncol(X_reduced)
  ) {
    return(empty_result())
  }

  # Fit observed models
  fit_full <- lm.fit(
    x = X_full,
    y = y
  )

  fit_reduced <- lm.fit(
    x = X_reduced,
    y = y
  )

  full_coefficients <- fit_full$coefficients

  beta_mt <- unname(
    full_coefficients["mt_PBS"]
  )

  beta_chr21 <- unname(
    full_coefficients["chr21_PBS"]
  )

  if (!is.finite(beta_mt)) {
    return(empty_result())
  }

  rss_full <- sum(
    fit_full$residuals^2
  )

  rss_reduced <- sum(
    fit_reduced$residuals^2
  )

  if (
    !is.finite(rss_full) ||
    !is.finite(rss_reduced) ||
    rss_reduced <= 0
  ) {
    return(empty_result())
  }

  partial_R2_mt <- (
    rss_reduced - rss_full
  ) / rss_reduced

  # Raw correlation
  cor_raw <- cor(
    d$nu_PBS,
    d$mt_PBS,
    use = "complete.obs"
  )

  # Residual correlation after controlling chr21 and region
  nu_neutral_fit <- lm(
    nu_PBS ~ chr21_PBS + region,
    data = d
  )

  mt_neutral_fit <- lm(
    mt_PBS ~ chr21_PBS + region,
    data = d
  )

  cor_residual <- cor(
    residuals(nu_neutral_fit),
    residuals(mt_neutral_fit),
    use = "complete.obs"
  )

  # Approximate VIF for mtPBS within this gene's population set
  mt_model_R2 <- summary(
    mt_neutral_fit
  )$r.squared

  gene_mt_vif <- 1 / (
    1 - mt_model_R2
  )

  # Parametric P value for reference
  full_lm <- lm(
    nu_PBS ~ mt_PBS + chr21_PBS + region,
    data = d
  )

  coefficient_table <- summary(
    full_lm
  )$coefficients

  p_parametric <- if (
    "mt_PBS" %in% rownames(coefficient_table)
  ) {
    coefficient_table[
      "mt_PBS",
      "Pr(>|t|)"
    ]
  } else {
    NA_real_
  }

  # Freedman-Lane permutation
  fitted_reduced <- fit_reduced$fitted.values
  residual_reduced <- fit_reduced$residuals

  # Precompute coefficient transformation:
  # beta = solve(X'X) X'y
  coefficient_operator <- solve(
    crossprod(X_full),
    t(X_full)
  )

  mt_row <- which(
    colnames(X_full) == "mt_PBS"
  )

  beta_permuted <- numeric(
    nperm
  )

  for (iteration in seq_len(nperm)) {

    permuted_residuals <- permute_within_region(
      residual_vector = residual_reduced,
      region_vector = d$region
    )

    y_permuted <- (
      fitted_reduced +
      permuted_residuals
    )

    permuted_coefficients <- (
      coefficient_operator %*%
      y_permuted
    )

    beta_permuted[iteration] <- permuted_coefficients[
      mt_row
    ]
  }

  p_perm <- (
    sum(
      abs(beta_permuted) >=
      abs(beta_mt)
    ) + 1
  ) / (
    nperm + 1
  )

  data.table(
    n_pop = n_population,
    n_AK = sum(d$region == "AK"),
    n_BC = sum(d$region == "BC"),
    beta_mt = beta_mt,
    beta_chr21 = beta_chr21,
    partial_R2_mt = partial_R2_mt,
    cor_raw = cor_raw,
    cor_residual = cor_residual,
    p_parametric = p_parametric,
    p_perm = p_perm,
    mt_vif = gene_mt_vif
  )
}

# ============================================================
# Run all genes
# ============================================================

cat(
  "\n[Run] Testing genes with",
  N_PERM,
  "permutations per gene\n"
)

result <- dat[
  ,
  test_one_gene(
    .SD,
    nperm = N_PERM
  ),
  by = gene
]

# Keep BH value for transparency, but candidate selection uses p_perm
result[, p_bh := p.adjust(
  p_perm,
  method = "BH"
)]

result[, candidate := (
  is.finite(beta_mt) &
  is.finite(p_perm) &
  beta_mt > 0 &
  p_perm < 0.05
)]

result[, direction := fifelse(
  !is.finite(beta_mt),
  "unresolved",
  fifelse(
    beta_mt > 0,
    "positive",
    "negative"
  )
)]

setorder(
  result,
  -candidate,
  p_perm,
  -partial_R2_mt,
  -beta_mt
)

# ============================================================
# Write results
# ============================================================

fwrite(
  result,
  OUT_TABLE,
  sep = "\t"
)

candidate_result <- result[
  candidate == TRUE
]

fwrite(
  candidate_result,
  OUT_CANDIDATES,
  sep = "\t"
)

cat(
  "\n[Result] Tested genes:",
  nrow(result),
  "\n"
)

cat(
  "[Result] Positive permutation candidates:",
  nrow(candidate_result),
  "\n"
)

cat("\n[Result] Candidate genes:\n")

print(
  candidate_result[
    ,
    .(
      gene,
      n_pop,
      beta_mt,
      beta_chr21,
      partial_R2_mt,
      cor_raw,
      cor_residual,
      p_parametric,
      p_perm,
      mt_vif
    )
  ]
)

# ============================================================
# Plot
# ============================================================

plot_dt <- result[
  is.finite(beta_mt) &
  is.finite(p_perm)
]

plot_dt[
  p_perm <= 0,
  p_perm := 1 / (
    N_PERM + 1
  )
]

plot_dt[, label_gene := ifelse(
  candidate |
  p_perm < 0.10,
  gene,
  NA_character_
)]

p <- ggplot(
  plot_dt,
  aes(
    x = beta_mt,
    y = -log10(p_perm)
  )
) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    color = "grey45"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "grey45"
  ) +
  geom_point(
    aes(
      color = candidate
    ),
    size = 2.8,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "grey65",
      "TRUE" = "#E66101"
    ),
    labels = c(
      "FALSE" = "Other genes",
      "TRUE" = "Candidate"
    ),
    name = NULL
  ) +
  theme_classic(
    base_size = 14
  ) +
  labs(
    x = "mtPBS effect after controlling for chr21 PBS and region",
    y = expression(-log[10](P[permutation])),
    title = "Chr21-controlled mitonuclear PBS associations"
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

# Label without requiring ggrepel
candidate_labels <- plot_dt[
  !is.na(label_gene)
]

if (nrow(candidate_labels) > 0) {
  p <- p +
    geom_text(
      data = candidate_labels,
      aes(label = label_gene),
      size = 3,
      vjust = -0.7,
      check_overlap = TRUE,
      show.legend = FALSE
    )
}

ggsave(
  OUT_PNG,
  p,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  OUT_PDF,
  p,
  width = 8,
  height = 6
)

# ============================================================
# Final summary
# ============================================================

cat("\n[Written]\n")
cat("All results: ", OUT_TABLE, "\n", sep = "")
cat("Candidates:  ", OUT_CANDIDATES, "\n", sep = "")
cat("Plot PNG:    ", OUT_PNG, "\n", sep = "")
cat("Plot PDF:    ", OUT_PDF, "\n", sep = "")
