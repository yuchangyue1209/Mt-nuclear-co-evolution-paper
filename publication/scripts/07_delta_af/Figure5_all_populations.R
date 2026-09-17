#!/usr/bin/env Rscript

# ============================================================
# Manuscript Figure 5: all-population primary analysis
#
# Primary analysis:
#   * 25 freshwater populations (11 Alaska + 14 British Columbia)
#   * 2 regional marine references (RS and SAY)
#   * one shared representative SNP per nuclear OXPHOS gene
#   * the same focal allele is used in every population
#
# Historical follow-up:
#   * five recently colonized populations are identified only after
#     the all-population SNP selection and candidate definition
#   * recent versus established summaries are written separately
# ============================================================

SCRIPT_VERSION <- "1.0.0"

# ------------------------------------------------------------
# Command-line interface, paths, and settings
# ------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) y else x
}

usage <- paste(
  "Usage:",
  "  Rscript Figure5_all_populations.R \\",
  "    --af <allele_frequency_table.tsv[.gz]> \\",
  "    --genes <oxphos72_gene_list.txt> \\",
  "    --outdir <output_directory> \\",
  "    [--min-depth 20] [--threshold 0.5] [--dpi 500]",
  "",
  "Required AF columns: chr, pos, gene, pop, af, depth.",
  "An optional focal_allele column is retained in output provenance.",
  sep = "\n"
)

parse_cli <- function(args) {
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat(usage, "\n")
    quit(save = "no", status = if (length(args)) 0L else 1L)
  }

  allowed <- c(
    "--af", "--genes", "--outdir",
    "--min-depth", "--threshold", "--dpi"
  )

  unknown <- args[grepl("^--", args) & !args %in% allowed]
  if (length(unknown)) {
    stop("Unknown option(s): ", paste(unique(unknown), collapse = ", "))
  }

  if (length(args) %% 2L != 0L) {
    stop("Every option must be followed by a value.\n\n", usage)
  }

  keys <- args[seq(1L, length(args), by = 2L)]
  values <- args[seq(2L, length(args), by = 2L)]

  if (anyDuplicated(keys)) {
    stop("Each command-line option may be supplied only once.")
  }

  out <- as.list(values)
  names(out) <- sub("^--", "", keys)
  out
}

CLI <- parse_cli(commandArgs(trailingOnly = TRUE))

required_options <- c("af", "genes", "outdir")
missing_options <- required_options[!required_options %in% names(CLI)]

if (length(missing_options)) {
  stop(
    "Missing required option(s): ",
    paste0("--", missing_options, collapse = ", "),
    "\n\n",
    usage
  )
}

AF_FILE <- normalizePath(CLI$af, mustWork = FALSE)
GENE_LIST <- normalizePath(CLI$genes, mustWork = FALSE)
OUTDIR <- normalizePath(CLI$outdir, mustWork = FALSE)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MIN_DEPTH <- as.numeric(CLI[["min-depth"]] %||% "20")
THRESHOLD <- as.numeric(CLI$threshold %||% "0.5")
DPI <- as.numeric(CLI$dpi %||% "500")

if (!is.finite(MIN_DEPTH) || MIN_DEPTH < 0) {
  stop("--min-depth must be a nonnegative number.")
}

if (!is.finite(THRESHOLD) || THRESHOLD < 0 || THRESHOLD > 1) {
  stop("--threshold must be between 0 and 1.")
}

if (!is.finite(DPI) || DPI <= 0) {
  stop("--dpi must be a positive number.")
}

required_packages <- c("data.table", "ggplot2", "ggrepel", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

AK_ESTABLISHED <- c(
  "FG", "LG", "SR", "SL", "TL", "WB", "WT", "WK"
)

AK_RECENT <- c("SC", "CH", "LB")

BC_ESTABLISHED <- c(
  "SWA", "THE", "JOE", "BEA", "MUC", "PYE",
  "AMO", "BOOT", "ECHO", "LAW", "GOS", "ROB"
)

BC_RECENT <- c("PACH", "FRED")

AK_FRESH <- c(AK_ESTABLISHED, AK_RECENT)
BC_FRESH <- c(BC_ESTABLISHED, BC_RECENT)

AK_MARINE <- "RS"
BC_MARINE <- "SAY"

EXPECTED_N <- c(
  AK = length(AK_FRESH),
  BC = length(BC_FRESH)
)

REGION_LEVELS <- c("AK", "BC")
REGION_LABELS <- c(AK = "Alaska", BC = "British Columbia")
REGION_COLORS <- c(AK = "#1796C4", BC = "#08A77B")
REGION_SHAPES <- c(AK = 16, BC = 17)
RECENT_OUTLINE <- "#E31A1C"
MARINE_COLOR <- "#7B3294"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

normalize_pop <- function(x) {
  x <- basename(as.character(x))
  x <- sub("_subset\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x <- sub(
    "^(?:[0-9]+_)?([A-Za-z]+)(?:_S[0-9]+)?$",
    "\\1",
    x,
    perl = TRUE
  )
  toupper(x)
}

normalize_gene <- function(x) {
  tolower(trimws(as.character(x)))
}

safe_cor_test <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) {
    return(list(n = length(x), estimate = NA_real_, p.value = NA_real_))
  }

  test <- cor.test(
    x,
    y,
    method = method,
    exact = FALSE
  )

  list(
    n = length(x),
    estimate = unname(test$estimate),
    p.value = test$p.value
  )
}

required_files <- c(AF_FILE, GENE_LIST)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files)) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

RUN_PARAMETERS <- data.table(
  parameter = c(
    "script_version", "af_file", "gene_list", "output_directory",
    "minimum_depth", "candidate_threshold", "figure_dpi"
  ),
  value = c(
    SCRIPT_VERSION, AF_FILE, GENE_LIST, OUTDIR,
    MIN_DEPTH, THRESHOLD, DPI
  )
)

fwrite(
  RUN_PARAMETERS,
  file.path(OUTDIR, "Figure5_run_parameters.tsv"),
  sep = "\t"
)

# ------------------------------------------------------------
# Read and validate the 72-gene AF table
# ------------------------------------------------------------

genes72 <- unique(normalize_gene(readLines(GENE_LIST, warn = FALSE)))
genes72 <- genes72[nzchar(genes72)]

if (length(genes72) != 72L) {
  stop("Expected 72 nuclear OXPHOS genes; found ", length(genes72))
}

AF <- fread(AF_FILE)

required_columns <- c(
  "chr", "pos", "gene", "pop", "af", "depth"
)

missing_columns <- setdiff(required_columns, names(AF))

if (length(missing_columns)) {
  stop(
    "AF table lacks required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

# Some finalized AF tables retain only the consistently oriented frequency
# column and do not repeat the focal-allele nucleotide. The nucleotide label is
# useful for provenance but is not required to calculate freshwater-minus-
# marine delta AF. In that case, retain an explicit missing-value column and
# report the limitation in the run log.
if (!"focal_allele" %in% names(AF)) {
  AF[, focal_allele := NA_character_]
  warning(
    "AF table has no focal_allele column. Using the existing, consistently ",
    "oriented af values; nucleotide identities will be recorded as NA."
  )
}

AF[, gene := normalize_gene(gene)]
AF[, pop := normalize_pop(pop)]
AF[, af := as.numeric(af)]
AF[, depth := as.numeric(depth)]
AF[, focal_allele := toupper(trimws(as.character(focal_allele)))]
AF[, chr := as.character(chr)]
AF[, pos := as.integer(pos)]

KEEP_POPS <- c(
  AK_FRESH,
  BC_FRESH,
  AK_MARINE,
  BC_MARINE
)

AF <- AF[
  gene %chin% genes72 &
    pop %chin% KEEP_POPS &
    is.finite(af) &
    is.finite(depth) &
    depth >= MIN_DEPTH
]

AF[, snp := paste(chr, pos, gene, sep = ":")]

duplicates <- AF[, .N, by = .(gene, snp, pop)][N > 1L]

if (nrow(duplicates)) {
  fwrite(
    duplicates,
    file.path(OUTDIR, "ERROR_duplicate_gene_snp_population_rows.tsv"),
    sep = "\t"
  )
  stop("Duplicate gene-SNP-population rows detected.")
}

if (any(!is.na(AF$focal_allele) & nzchar(AF$focal_allele))) {
  focal_audit <- AF[
    !is.na(focal_allele) & nzchar(focal_allele),
    .(n_focal_alleles = uniqueN(focal_allele)),
    by = snp
  ]

  if (focal_audit[n_focal_alleles != 1L, .N] > 0L) {
    fwrite(
      focal_audit[n_focal_alleles != 1L],
      file.path(OUTDIR, "ERROR_inconsistent_focal_alleles.tsv"),
      sep = "\t"
    )
    stop("At least one SNP has inconsistent focal-allele definitions.")
  }
}

missing_populations <- setdiff(KEEP_POPS, unique(AF$pop))

if (length(missing_populations)) {
  stop(
    "The following required populations are absent after filtering: ",
    paste(missing_populations, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Population metadata
# ------------------------------------------------------------

POP_META <- rbindlist(list(
  data.table(
    pop = AK_ESTABLISHED,
    region = "AK",
    habitat = "freshwater",
    population_history = "established"
  ),
  data.table(
    pop = AK_RECENT,
    region = "AK",
    habitat = "freshwater",
    population_history = "recent"
  ),
  data.table(
    pop = BC_ESTABLISHED,
    region = "BC",
    habitat = "freshwater",
    population_history = "established"
  ),
  data.table(
    pop = BC_RECENT,
    region = "BC",
    habitat = "freshwater",
    population_history = "recent"
  ),
  data.table(
    pop = AK_MARINE,
    region = "AK",
    habitat = "marine",
    population_history = "marine"
  ),
  data.table(
    pop = BC_MARINE,
    region = "BC",
    habitat = "marine",
    population_history = "marine"
  )
))

AF <- merge(AF, POP_META, by = "pop", all.x = TRUE)

if (AF[is.na(region), .N] > 0L) {
  stop("Population metadata could not be assigned to all AF rows.")
}

# ------------------------------------------------------------
# Calculate freshwater-minus-marine delta AF for all populations
# ------------------------------------------------------------

MARINE <- AF[habitat == "marine", .(
  region,
  snp,
  marine_pop = pop,
  marine_af = af,
  marine_depth = depth
)]

if (MARINE[, anyDuplicated(paste(region, snp))] != 0L) {
  stop("More than one regional marine reference row exists for a SNP.")
}

FRESH <- AF[habitat == "freshwater", .(
  region,
  population_history,
  snp,
  chr,
  pos,
  gene,
  pop,
  af,
  depth,
  focal_allele
)]

DELTA <- merge(
  FRESH,
  MARINE,
  by = c("region", "snp"),
  all = FALSE
)

DELTA[, deltaAF := af - marine_af]
DELTA[, abs_deltaAF := abs(deltaAF)]

setorder(DELTA, gene, snp, region, pop)

fwrite(
  DELTA,
  file.path(OUTDIR, "deltaAF_all_freshwater_populations.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

# ------------------------------------------------------------
# Score every SNP using all freshwater populations
# ------------------------------------------------------------

REGIONAL <- DELTA[, {
  med <- median(deltaAF, na.rm = TRUE)
  n_positive <- sum(deltaAF > 0, na.rm = TRUE)
  n_negative <- sum(deltaAF < 0, na.rm = TRUE)
  expected_n <- unname(EXPECTED_N[region[1]])

  list(
    n_pop = uniqueN(pop),
    expected_n_pop = expected_n,
    complete_coverage = uniqueN(pop) == expected_n,
    median_deltaAF = med,
    abs_median_deltaAF = abs(med),
    mean_deltaAF = mean(deltaAF, na.rm = TRUE),
    abs_mean_deltaAF = abs(mean(deltaAF, na.rm = TRUE)),
    directional_concordance = max(n_positive, n_negative) / .N,
    mean_depth = mean(depth, na.rm = TRUE),
    min_depth = min(depth, na.rm = TRUE),
    median_freshwater_af = median(af, na.rm = TRUE),
    marine_af = marine_af[1],
    marine_pop = marine_pop[1],
    marine_depth = marine_depth[1],
    focal_allele = focal_allele[1],
    chr = chr[1],
    pos = pos[1]
  )
}, by = .(region, gene, snp)]

WIDE <- dcast(
  REGIONAL,
  gene + snp + chr + pos ~ region,
  value.var = c(
    "n_pop",
    "complete_coverage",
    "median_deltaAF",
    "abs_median_deltaAF",
    "mean_deltaAF",
    "abs_mean_deltaAF",
    "directional_concordance",
    "mean_depth",
    "min_depth",
    "median_freshwater_af",
    "marine_af",
    "marine_pop",
    "marine_depth"
  )
)

WIDE_COMPLETE <- WIDE[
  complete_coverage_AK == TRUE &
    complete_coverage_BC == TRUE &
    is.finite(median_deltaAF_AK) &
    is.finite(median_deltaAF_BC)
]

WIDE_COMPLETE[, same_nonzero_direction := (
  sign(median_deltaAF_AK) == sign(median_deltaAF_BC) &
    sign(median_deltaAF_AK) != 0
)]

WIDE_COMPLETE[, shared_score := pmin(
  abs_median_deltaAF_AK,
  abs_median_deltaAF_BC
)]

WIDE_COMPLETE[, joint_mean_abs_median := (
  abs_median_deltaAF_AK + abs_median_deltaAF_BC
) / 2]

WIDE_COMPLETE[, minimum_directional_concordance := pmin(
  directional_concordance_AK,
  directional_concordance_BC
)]

WIDE_COMPLETE[, minimum_mean_depth := pmin(
  mean_depth_AK,
  mean_depth_BC
)]

WIDE_COMPLETE[, has_concordant_candidate := any(
  same_nonzero_direction
), by = gene]

WIDE_COMPLETE[, eligible_for_selection := fifelse(
  has_concordant_candidate,
  same_nonzero_direction,
  TRUE
)]

WIDE_COMPLETE[, selection_class := fifelse(
  has_concordant_candidate,
  "concordant_preferred",
  "opposite_or_zero_fallback"
)]

fwrite(
  WIDE_COMPLETE,
  file.path(OUTDIR, "all_complete_shared_SNP_scores.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

# ------------------------------------------------------------
# Coverage audit
# ------------------------------------------------------------

COVERAGE <- WIDE_COMPLETE[, .(
  n_shared_complete_snps = .N,
  n_concordant_shared_snps = sum(same_nonzero_direction)
), by = gene]

COVERAGE <- merge(
  data.table(gene = genes72),
  COVERAGE,
  by = "gene",
  all.x = TRUE
)

COVERAGE[is.na(n_shared_complete_snps), `:=`(
  n_shared_complete_snps = 0L,
  n_concordant_shared_snps = 0L
)]

setorder(COVERAGE, gene)

fwrite(
  COVERAGE,
  file.path(OUTDIR, "shared_SNP_complete_coverage_audit.tsv"),
  sep = "\t"
)

missing_genes <- COVERAGE[n_shared_complete_snps == 0L]

if (nrow(missing_genes)) {
  print(missing_genes)
  stop(
    "At least one gene has no SNP with complete coverage across all ",
    "11 AK and 14 BC freshwater populations. See the coverage audit."
  )
}

# ------------------------------------------------------------
# Select one shared SNP per gene
# ------------------------------------------------------------

ELIGIBLE <- WIDE_COMPLETE[eligible_for_selection == TRUE]

setorder(
  ELIGIBLE,
  gene,
  -shared_score,
  -joint_mean_abs_median,
  -minimum_directional_concordance,
  -minimum_mean_depth,
  chr,
  pos
)

TOP <- ELIGIBLE[, .SD[1], by = gene]

TOP[, selection_method := paste0(
  "all_25_freshwater_populations;prefer_same_direction_then_",
  "maximize_minimum_absolute_regional_median_deltaAF"
)]

setorder(TOP, gene)

if (nrow(TOP) != 72L) {
  stop("Expected 72 representative SNPs; found ", nrow(TOP))
}

KEY <- TOP[, .(gene, snp)]

LONG <- merge(
  DELTA,
  KEY,
  by = c("gene", "snp"),
  all = FALSE
)

setorder(LONG, gene, region, population_history, pop)

SUMMARY <- LONG[, {
  med <- median(deltaAF, na.rm = TRUE)
  n_positive <- sum(deltaAF > 0, na.rm = TRUE)
  n_negative <- sum(deltaAF < 0, na.rm = TRUE)

  list(
    n_pop = uniqueN(pop),
    median_deltaAF = med,
    abs_median_deltaAF = abs(med),
    mean_deltaAF = mean(deltaAF, na.rm = TRUE),
    directional_concordance = max(n_positive, n_negative) / .N,
    median_freshwater_af = median(af, na.rm = TRUE),
    marine_af = marine_af[1],
    marine_pop = marine_pop[1],
    focal_allele = focal_allele[1],
    shared_snp = snp[1],
    chr = chr[1],
    pos = pos[1]
  )
}, by = .(gene, region)]

SUMMARY <- merge(
  SUMMARY,
  TOP[, .(
    gene,
    shared_score,
    same_nonzero_direction,
    selection_class
  )],
  by = "gene",
  all.x = TRUE
)

setorder(SUMMARY, gene, region)

fwrite(
  TOP,
  file.path(OUTDIR, "Figure5_representative_SNP_all_populations.tsv"),
  sep = "\t"
)

fwrite(
  LONG,
  file.path(OUTDIR, "Figure5_population_data_all_populations.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  SUMMARY,
  file.path(OUTDIR, "Figure5_region_summary_all_populations.tsv"),
  sep = "\t"
)

# ------------------------------------------------------------
# Main statistical results using all freshwater populations
# ------------------------------------------------------------

A_DATA <- dcast(
  SUMMARY,
  gene + selection_class + shared_score ~ region,
  value.var = "median_deltaAF"
)

setnames(A_DATA, c("AK", "BC"), c(
  "AK_median_deltaAF",
  "BC_median_deltaAF"
))

A_DATA <- A_DATA[
  is.finite(AK_median_deltaAF) &
    is.finite(BC_median_deltaAF)
]

A_DATA[, same_direction := (
  sign(AK_median_deltaAF) == sign(BC_median_deltaAF) &
    sign(AK_median_deltaAF) != 0
)]

CANDIDATE_GENES <- sort(unique(
  SUMMARY[abs_median_deltaAF >= THRESHOLD, gene]
))

if (!length(CANDIDATE_GENES)) {
  stop("No genes passed the all-population candidate threshold.")
}

A_DATA[, candidate := gene %chin% CANDIDATE_GENES]

PEARSON <- cor.test(
  A_DATA$AK_median_deltaAF,
  A_DATA$BC_median_deltaAF,
  method = "pearson"
)

SPEARMAN <- cor.test(
  A_DATA$AK_median_deltaAF,
  A_DATA$BC_median_deltaAF,
  method = "spearman",
  exact = FALSE
)

DIRECTION_TEST <- binom.test(
  x = sum(A_DATA$same_direction),
  n = nrow(A_DATA),
  p = 0.5,
  alternative = "two.sided"
)

COR_STATS <- data.table(
  test = c("Pearson", "Spearman"),
  estimate = c(
    unname(PEARSON$estimate),
    unname(SPEARMAN$estimate)
  ),
  p_value = c(PEARSON$p.value, SPEARMAN$p.value),
  n_gene = nrow(A_DATA),
  populations = "all 25 freshwater populations"
)

DIRECTION_STATS <- data.table(
  n_gene = nrow(A_DATA),
  same_direction_genes = sum(A_DATA$same_direction),
  proportion_same_direction = mean(A_DATA$same_direction),
  null_probability = 0.5,
  p_value = DIRECTION_TEST$p.value,
  confidence_lower = DIRECTION_TEST$conf.int[1],
  confidence_upper = DIRECTION_TEST$conf.int[2]
)

fwrite(
  COR_STATS,
  file.path(OUTDIR, "Figure5A_correlation_statistics_all_populations.tsv"),
  sep = "\t"
)

fwrite(
  DIRECTION_STATS,
  file.path(OUTDIR, "Figure5A_directional_concordance_all_populations.tsv"),
  sep = "\t"
)

fwrite(
  A_DATA,
  file.path(OUTDIR, "Figure5A_plot_data_all_populations.tsv"),
  sep = "\t"
)

# ------------------------------------------------------------
# Separate historical analysis: recent versus established
# ------------------------------------------------------------

HISTORY_SUMMARY <- LONG[, .(
  n_pop = uniqueN(pop),
  median_af = median(af, na.rm = TRUE),
  median_deltaAF = median(deltaAF, na.rm = TRUE),
  mean_deltaAF = mean(deltaAF, na.rm = TRUE)
), by = .(
  gene,
  region,
  population_history,
  snp,
  marine_af,
  marine_pop
)]

HISTORY_CANDIDATES <- HISTORY_SUMMARY[gene %chin% CANDIDATE_GENES]

EST_MEDIAN <- HISTORY_CANDIDATES[
  population_history == "established",
  .(
    gene,
    region,
    established_median_deltaAF = median_deltaAF,
    established_median_af = median_af
  )
]

RECENT_LONG <- merge(
  LONG[
    gene %chin% CANDIDATE_GENES &
      population_history == "recent"
  ],
  EST_MEDIAN,
  by = c("gene", "region"),
  all.x = TRUE
)

RECENT_LONG[, direction_matches_established := (
  sign(deltaAF) == sign(established_median_deltaAF) &
    sign(established_median_deltaAF) != 0
)]

RECENT_POP_STATS <- RECENT_LONG[, {
  test <- safe_cor_test(
    deltaAF,
    established_median_deltaAF,
    method = "spearman"
  )

  list(
    n_candidate_genes = uniqueN(gene),
    spearman_rho = test$estimate,
    p_value = test$p.value,
    concordant_genes = sum(direction_matches_established, na.rm = TRUE),
    proportion_directionally_concordant = mean(
      direction_matches_established,
      na.rm = TRUE
    ),
    median_recent_deltaAF = median(deltaAF, na.rm = TRUE)
  )
}, by = .(pop, region)]

HISTORY_WIDE <- dcast(
  HISTORY_CANDIDATES,
  gene + region + snp + marine_af + marine_pop ~ population_history,
  value.var = "median_deltaAF"
)

if (all(c("established", "recent") %in% names(HISTORY_WIDE))) {
  HISTORY_WIDE[, recent_minus_established := recent - established]
}

fwrite(
  HISTORY_SUMMARY,
  file.path(OUTDIR, "Figure5_history_group_summary.tsv"),
  sep = "\t"
)

fwrite(
  RECENT_LONG,
  file.path(OUTDIR, "Figure5_recent_population_candidate_data.tsv"),
  sep = "\t"
)

fwrite(
  RECENT_POP_STATS,
  file.path(OUTDIR, "Figure5_recent_vs_established_statistics.tsv"),
  sep = "\t"
)

fwrite(
  HISTORY_WIDE,
  file.path(OUTDIR, "Figure5_recent_vs_established_gene_summary.tsv"),
  sep = "\t"
)

# ------------------------------------------------------------
# Plot theme
# ------------------------------------------------------------

theme_fig <- theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(colour = "black", size = 16),
    axis.text = element_text(colour = "black", size = 13),
    axis.line = element_line(colour = "black", linewidth = 0.7),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = grid::unit(0.16, "cm"),
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.key.size = grid::unit(0.7, "cm"),
    plot.caption = element_text(size = 10.5, hjust = 0, colour = "black"),
    plot.margin = margin(10, 14, 10, 12)
  )

# ------------------------------------------------------------
# Panel A
# ------------------------------------------------------------

stat_label <- sprintf(
  "Pearson r = %.2f\nP = %.2g",
  unname(PEARSON$estimate),
  PEARSON$p.value
)

lims_A <- range(c(
  A_DATA$AK_median_deltaAF,
  A_DATA$BC_median_deltaAF
), finite = TRUE)

pad_A <- max(diff(lims_A) * 0.12, 0.08)
lims_A <- c(lims_A[1] - pad_A, lims_A[2] + pad_A)

PANEL_A <- ggplot(
  A_DATA,
  aes(AK_median_deltaAF, BC_median_deltaAF)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey70"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey70"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted",
    linewidth = 0.7,
    colour = "grey45"
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "black",
    fill = "grey78",
    linewidth = 0.9
  ) +
  geom_point(
    data = A_DATA[candidate == FALSE],
    shape = 21,
    fill = "grey80",
    colour = "grey35",
    size = 3.7,
    stroke = 0.55
  ) +
  geom_point(
    data = A_DATA[candidate == TRUE],
    shape = 21,
    fill = "#F4A261",
    colour = "black",
    size = 4.5,
    stroke = 0.9
  ) +
  geom_text_repel(
    data = A_DATA[candidate == TRUE],
    aes(label = gene),
    fontface = "italic",
    size = 4.2,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.4,
    segment.colour = "grey35",
    box.padding = 0.7,
    point.padding = 0.45,
    force = 2,
    max.time = 5,
    max.iter = 20000,
    seed = 100
  ) +
  annotate(
    "text",
    x = lims_A[1] + 0.03 * diff(lims_A),
    y = lims_A[2] - 0.03 * diff(lims_A),
    label = stat_label,
    hjust = 0,
    vjust = 1,
    size = 4.5
  ) +
  coord_equal(
    xlim = lims_A,
    ylim = lims_A,
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    title = "A",
    x = "Alaska median ΔAF",
    y = "British Columbia median ΔAF"
  ) +
  theme_fig

# ------------------------------------------------------------
# Panel B data and ordering
# ------------------------------------------------------------

B_SUM <- SUMMARY[gene %chin% CANDIDATE_GENES]
B_LONG <- LONG[gene %chin% CANDIDATE_GENES]

B_WIDE <- dcast(
  B_SUM,
  gene + selection_class + shared_score ~ region,
  value.var = "median_deltaAF"
)

B_WIDE[, same_direction := sign(AK) == sign(BC) & sign(AK) != 0]

B_WIDE[, pattern_order := fifelse(
  same_direction & AK > 0 & BC > 0,
  1L,
  fifelse(
    same_direction & AK < 0 & BC < 0,
    2L,
    3L
  )
)]

setorder(B_WIDE, pattern_order, -shared_score, gene)
GENE_ORDER <- B_WIDE$gene

B_SUM[, gene_order := match(gene, GENE_ORDER)]
B_LONG[, gene_order := match(gene, GENE_ORDER)]

GENE_SPACING <- 2.50
REGION_OFFSET <- c(AK = -0.40, BC = 0.40)

B_SUM[, gene_center := (gene_order - 1) * GENE_SPACING + 1]
B_LONG[, gene_center := (gene_order - 1) * GENE_SPACING + 1]

B_SUM[, x_position := (
  gene_center + REGION_OFFSET[as.character(region)]
)]

B_LONG[, x_position := (
  gene_center + REGION_OFFSET[as.character(region)]
)]

setorder(B_LONG, gene_order, region, population_history, pop)

B_LONG[, jitter_offset := if (.N == 1L) {
  0
} else {
  seq(-0.25, 0.25, length.out = .N)
}, by = .(gene, region)]

B_LONG[, x_jitter := x_position + jitter_offset]

B_SUM[, delta_label := sprintf("%+.2f", median_deltaAF)]
B_SUM[, label_y := fifelse(median_deltaAF >= 0, 1.035, -0.035)]
B_SUM[, label_vjust := fifelse(median_deltaAF >= 0, 0, 1)]

GENE_LABELS <- unique(B_SUM[, .(
  gene,
  gene_order,
  gene_center,
  selection_class
)])

GENE_LABELS[, gene_label := fifelse(
  selection_class == "opposite_or_zero_fallback",
  paste0(gene, "†"),
  gene
)]

setorder(GENE_LABELS, gene_order)

SEPARATORS <- if (nrow(GENE_LABELS) > 1L) {
  (
    GENE_LABELS$gene_center[-nrow(GENE_LABELS)] +
      GENE_LABELS$gene_center[-1L]
  ) / 2
} else {
  numeric(0)
}

REGION_LABEL_DATA <- B_SUM[, .(
  x_position,
  region_label = paste0(
    as.character(region),
    fifelse(abs_median_deltaAF >= THRESHOLD, "*", "")
  )
)]

HAS_FALLBACK <- any(
  B_SUM$selection_class == "opposite_or_zero_fallback"
)

PANEL_B_CAPTION <- paste0(
  "Panel B shows the ten high-effect genes identified from all 25 freshwater populations. ",
  "Blue circles and green triangles denote Alaska and British Columbia; red outlines identify ",
  "recently colonized populations. Purple and black solid segments denote the regional marine ",
  "frequency and the median across all freshwater populations, respectively. Signed numbers ",
  "report regional median ΔAF; * indicates |median ΔAF| ≥ 0.5.",
  if (HAS_FALLBACK) paste0(
    " † No same-direction SNP with complete coverage was available; ",
    "the highest shared-score fallback SNP is shown."
  ) else ""
)

# ------------------------------------------------------------
# Panel B
# ------------------------------------------------------------

PANEL_B <- ggplot() +
  geom_hline(
    yintercept = seq(0, 1, 0.25),
    linewidth = 0.4,
    colour = "grey89"
  ) +
  geom_vline(
    xintercept = SEPARATORS,
    linewidth = 0.4,
    colour = "grey85"
  ) +
  geom_segment(
    data = B_SUM,
    aes(
      x = x_position - 0.31,
      xend = x_position + 0.31,
      y = marine_af,
      yend = marine_af
    ),
    linetype = "solid",
    linewidth = 1.05,
    colour = MARINE_COLOR
  ) +
  geom_point(
    data = B_LONG[population_history == "established"],
    aes(
      x = x_jitter,
      y = af,
      shape = region,
      colour = region
    ),
    size = 3.1,
    stroke = 0.35,
    alpha = 0.62
  ) +
  geom_point(
    data = B_LONG[
      population_history == "recent" & region == "AK"
    ],
    aes(x = x_jitter, y = af),
    shape = 21,
    size = 3.8,
    stroke = 1.0,
    fill = REGION_COLORS[["AK"]],
    colour = RECENT_OUTLINE
  ) +
  geom_point(
    data = B_LONG[
      population_history == "recent" & region == "BC"
    ],
    aes(x = x_jitter, y = af),
    shape = 24,
    size = 4.1,
    stroke = 1.0,
    fill = REGION_COLORS[["BC"]],
    colour = RECENT_OUTLINE
  ) +
  geom_segment(
    data = B_SUM,
    aes(
      x = x_position - 0.33,
      xend = x_position + 0.33,
      y = median_freshwater_af,
      yend = median_freshwater_af
    ),
    linewidth = 1.6,
    colour = "black",
    lineend = "round"
  ) +
  geom_text(
    data = B_SUM,
    aes(
      x = x_position,
      y = label_y,
      label = delta_label,
      vjust = label_vjust
    ),
    size = 3.8,
    fontface = "bold"
  ) +
  geom_text(
    data = REGION_LABEL_DATA,
    aes(
      x = x_position,
      y = -0.125,
      label = region_label
    ),
    size = 3.8,
    fontface = "bold",
    colour = "grey20"
  ) +
  geom_text(
    data = GENE_LABELS,
    aes(
      x = gene_center,
      y = -0.225,
      label = gene_label
    ),
    fontface = "italic",
    size = 4.0,
    colour = "black"
  ) +
  scale_shape_manual(
    values = REGION_SHAPES,
    labels = REGION_LABELS
  ) +
  scale_colour_manual(
    values = REGION_COLORS,
    labels = REGION_LABELS
  ) +
  scale_x_continuous(
    breaks = NULL,
    expand = expansion(add = c(0.75, 0.75))
  ) +
  scale_y_continuous(
    limits = c(-0.27, 1.10),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "B",
    x = NULL,
    y = "Focal-allele frequency",
    caption = PANEL_B_CAPTION
  ) +
  theme_fig +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    plot.margin = margin(12, 16, 64, 12)
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4, alpha = 0.9)
    ),
    shape = guide_legend(
      override.aes = list(size = 4, alpha = 0.9)
    )
  )

# ------------------------------------------------------------
# Export figures and plot data
# ------------------------------------------------------------

fwrite(
  B_SUM,
  file.path(OUTDIR, "Figure5B_plot_summary_all_populations.tsv"),
  sep = "\t"
)

fwrite(
  B_LONG,
  file.path(OUTDIR, "Figure5B_plot_population_data_all_populations.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

ggsave(
  file.path(OUTDIR, "Figure5A_all_populations.pdf"),
  PANEL_A,
  width = 7.0,
  height = 6.3,
  units = "in",
  device = grDevices::cairo_pdf
)

ggsave(
  file.path(OUTDIR, "Figure5A_all_populations.png"),
  PANEL_A,
  width = 7.0,
  height = 6.3,
  units = "in",
  dpi = DPI,
  bg = "white"
)

ggsave(
  file.path(OUTDIR, "Figure5B_all_populations.pdf"),
  PANEL_B,
  width = 15.5,
  height = 7.1,
  units = "in",
  device = grDevices::cairo_pdf,
  limitsize = FALSE
)

ggsave(
  file.path(OUTDIR, "Figure5B_all_populations.png"),
  PANEL_B,
  width = 15.5,
  height = 7.1,
  units = "in",
  dpi = DPI,
  bg = "white",
  limitsize = FALSE
)

FIGURE5 <- PANEL_A / PANEL_B +
  plot_layout(heights = c(0.92, 1.08))

ggsave(
  file.path(OUTDIR, "Figure5_all_populations.pdf"),
  FIGURE5,
  width = 15.8,
  height = 12.6,
  units = "in",
  device = grDevices::cairo_pdf,
  limitsize = FALSE
)

ggsave(
  file.path(OUTDIR, "Figure5_all_populations.png"),
  FIGURE5,
  width = 15.8,
  height = 12.6,
  units = "in",
  dpi = DPI,
  bg = "white",
  limitsize = FALSE
)

# ------------------------------------------------------------
# Final audit
# ------------------------------------------------------------

cat("\n===== Figure 5 all-population audit =====\n")
cat("Freshwater populations in primary analysis:",
    uniqueN(LONG$pop), "\n")
cat("  Alaska:", uniqueN(LONG[region == "AK"]$pop), "\n")
cat("  British Columbia:", uniqueN(LONG[region == "BC"]$pop), "\n")
cat("Established freshwater populations:",
    uniqueN(LONG[population_history == "established"]$pop), "\n")
cat("Recent freshwater populations:",
    uniqueN(LONG[population_history == "recent"]$pop), "\n")
cat("Genes in Panel A:", nrow(A_DATA), "\n")
cat("All-population high-effect genes:", length(CANDIDATE_GENES), "\n")
cat("High-effect genes displayed in Panel B:", length(CANDIDATE_GENES), "\n")
cat("Pearson r:", unname(PEARSON$estimate), "\n")
cat("Pearson P:", PEARSON$p.value, "\n")
cat("Same-direction genes:", sum(A_DATA$same_direction),
    "of", nrow(A_DATA), "\n")
cat("Directional binomial P:", DIRECTION_TEST$p.value, "\n")
cat("All-population high-effect genes:\n")
print(CANDIDATE_GENES)
cat("\nRecent-versus-established follow-up:\n")
print(RECENT_POP_STATS[order(region, pop)])
cat("\nOutput directory:", OUTDIR, "\n")

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "Figure5_R_sessionInfo.txt")
)
