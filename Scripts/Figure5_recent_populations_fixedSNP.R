#!/usr/bin/env Rscript

# Evaluate the five recently colonized populations at the representative SNPs
# selected independently for the established-population Figure 5 analysis.
# This script does NOT reselect SNPs or redefine Figure 5 candidate genes.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Paths and settings
# -----------------------------------------------------------------------------

BASE <- "/mnt/spareHD_2/nu_287/q2_parallelism"
FIG5_DIR <- file.path(BASE, "q2_deltaAF_withAMO_noLB_Jesse")

IN_TOP <- file.path(
  FIG5_DIR,
  "Fig3B_sharedRepresentativeSNP.withAMO_noLB_Jesse.tsv"
)
IN_SUM <- file.path(
  FIG5_DIR,
  "Fig3B_sharedSNP_regionSummary.withAMO_noLB_Jesse.tsv"
)

# Full 72-gene AF table containing all populations, including the recent ones.
IN_AF <- file.path(BASE, "af_long_final_72genes_subunit_with_si.tsv.gz")

OUTDIR <- file.path(FIG5_DIR, "Figure5_recent_populations_fixedSNP")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

THRESHOLD <- 0.50
DPI <- 500

RECENT_META <- data.table(
  pop = c("SC", "CH", "LB", "PACH", "FRED"),
  region = c("AK", "AK", "AK", "BC", "BC")
)

REGION_COLORS <- c(AK = "#1796C4", BC = "#08A77B")
POP_SHAPES <- c(SC = 21, CH = 22, LB = 24, PACH = 23, FRED = 25)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

normalize_pop <- function(x) {
  x <- basename(as.character(x))
  x <- sub("_subset\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("^(?:[0-9]+_)?([A-Za-z]+)(?:_S[0-9]+)?$", "\\1", x, perl = TRUE)
  toupper(x)
}

normalize_gene <- function(x) tolower(trimws(as.character(x)))

required_files <- c(IN_TOP, IN_SUM, IN_AF)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

# -----------------------------------------------------------------------------
# Read the original Figure 5 SNP choices and established-population summaries
# -----------------------------------------------------------------------------

TOP <- fread(IN_TOP)
SUM_RAW <- fread(IN_SUM)

required_top <- c("gene", "snp")
required_sum <- c(
  "gene", "region", "median_deltaAF", "abs_median_deltaAF",
  "median_freshwater_af", "marine_af", "marine_pop", "focal_allele"
)

if (!all(required_top %in% names(TOP))) {
  stop("IN_TOP lacks columns: ", paste(setdiff(required_top, names(TOP)), collapse = ", "))
}
if (!all(required_sum %in% names(SUM_RAW))) {
  stop("IN_SUM lacks columns: ", paste(setdiff(required_sum, names(SUM_RAW)), collapse = ", "))
}

TOP[, gene := normalize_gene(gene)]
TOP[, snp := as.character(snp)]
SUM_RAW[, gene := normalize_gene(gene)]
SUM_RAW[, region := toupper(trimws(as.character(region)))]

# Candidate definition is identical to the established-population Figure 5.
CANDIDATE_GENES <- sort(unique(
  SUM_RAW[abs_median_deltaAF >= THRESHOLD, gene]
))
if (!length(CANDIDATE_GENES)) stop("No Figure 5 candidate genes were identified.")

SELECTED <- unique(TOP[gene %chin% CANDIDATE_GENES, .(
  gene,
  selected_snp = snp
)])

if (SELECTED[, anyDuplicated(gene)] != 0L) {
  stop("More than one representative SNP was found for at least one candidate gene.")
}

# Attach the focal allele. It should be the same in AK and BC for a shared SNP.
FOCAL <- unique(SUM_RAW[gene %chin% CANDIDATE_GENES, .(
  gene,
  focal_allele = toupper(as.character(focal_allele))
)])
if (FOCAL[, uniqueN(focal_allele), by = gene][V1 > 1L, .N] > 0L) {
  stop("AK and BC do not use the same focal allele for at least one candidate gene.")
}
SELECTED <- merge(SELECTED, FOCAL, by = "gene", all.x = TRUE)

ESTABLISHED <- unique(SUM_RAW[gene %chin% CANDIDATE_GENES, .(
  gene,
  region,
  established_median_deltaAF = median_deltaAF,
  established_median_af = median_freshwater_af,
  marine_af,
  marine_pop
)])

# -----------------------------------------------------------------------------
# Read the full AF table and extract only the fixed SNPs in recent populations
# -----------------------------------------------------------------------------

AF <- fread(IN_AF)

required_af <- c("gene", "pop", "af")
if (!all(required_af %in% names(AF))) {
  stop("IN_AF lacks columns: ", paste(setdiff(required_af, names(AF)), collapse = ", "))
}

AF[, gene := normalize_gene(gene)]
AF[, pop := normalize_pop(pop)]
AF[, af := as.numeric(af)]
if ("depth" %in% names(AF)) AF[, depth := as.numeric(depth)]

if (!"snp" %in% names(AF)) {
  chr_col <- intersect(c("chr", "chromosome"), names(AF))[1]
  pos_col <- intersect(c("pos", "position"), names(AF))[1]
  if (is.na(chr_col) || is.na(pos_col)) {
    stop("IN_AF must contain either snp or chromosome/position columns.")
  }
  AF[, snp := paste(
    as.character(get(chr_col)),
    as.character(get(pos_col)),
    gene,
    sep = ":"
  )]
} else {
  AF[, snp := as.character(snp)]
}

# Match primarily by gene plus genomic coordinates. Rebuild a coordinate key so
# differences in gene-symbol capitalization inside SNP strings cannot break joins.
coord_key <- function(x) {
  parts <- tstrsplit(as.character(x), ":", fixed = TRUE)
  if (length(parts) < 2L) stop("Unexpected SNP identifier format; expected chr:pos[:gene].")
  paste(parts[[1]], parts[[2]], sep = ":")
}

SELECTED[, coord := coord_key(selected_snp)]
AF[, coord := coord_key(snp)]

RECENT <- merge(
  AF[pop %chin% RECENT_META$pop & gene %chin% CANDIDATE_GENES],
  SELECTED[, .(gene, coord, selected_snp, expected_focal_allele = focal_allele)],
  by = c("gene", "coord"),
  all = FALSE
)
RECENT <- merge(RECENT, RECENT_META, by = "pop", all.x = TRUE)

if (!nrow(RECENT)) {
  stop("No recent-population observations matched the fixed Figure 5 SNPs.")
}

# Confirm that AF is measured for the same focal allele as in Figure 5.
if ("focal_allele" %in% names(RECENT)) {
  RECENT[, observed_focal_allele := toupper(as.character(focal_allele))]
  mismatch <- RECENT[
    !is.na(observed_focal_allele) & !is.na(expected_focal_allele) &
      observed_focal_allele != expected_focal_allele
  ]
  if (nrow(mismatch)) {
    fwrite(
      mismatch,
      file.path(OUTDIR, "ERROR_focal_allele_mismatches.tsv"),
      sep = "\t"
    )
    stop(
      "The full AF table uses a different focal allele for at least one fixed SNP. ",
      "See ERROR_focal_allele_mismatches.tsv; do not interpret the plot."
    )
  }
}

# One observation per population and selected SNP is expected.
duplicates <- RECENT[, .N, by = .(gene, selected_snp, pop)][N > 1L]
if (nrow(duplicates)) {
  fwrite(duplicates, file.path(OUTDIR, "ERROR_duplicate_rows.tsv"), sep = "\t")
  stop("Duplicate gene-SNP-population rows found; see ERROR_duplicate_rows.tsv.")
}

RECENT <- merge(
  RECENT,
  ESTABLISHED,
  by = c("gene", "region"),
  all.x = TRUE
)

RECENT[, recent_deltaAF := af - marine_af]
RECENT[, direction_concordant := fifelse(
  !is.finite(recent_deltaAF) | !is.finite(established_median_deltaAF),
  NA,
  sign(recent_deltaAF) == sign(established_median_deltaAF)
)]
RECENT[, distance_to_established_median := abs(af - established_median_af)]

# -----------------------------------------------------------------------------
# Audit completeness and summarize correspondence with established populations
# -----------------------------------------------------------------------------

EXPECTED <- CJ(gene = CANDIDATE_GENES, pop = RECENT_META$pop, unique = TRUE)
EXPECTED <- merge(EXPECTED, RECENT_META, by = "pop", all.x = TRUE)
OBSERVED <- unique(RECENT[, .(gene, pop)])
COMPLETENESS <- merge(
  EXPECTED,
  OBSERVED[, observed := TRUE],
  by = c("gene", "pop"),
  all.x = TRUE
)
COMPLETENESS[is.na(observed), observed := FALSE]

POP_SUMMARY <- RECENT[, .(
  genes_observed = uniqueN(gene),
  concordant_genes = sum(direction_concordant %in% TRUE, na.rm = TRUE),
  discordant_genes = sum(direction_concordant %in% FALSE, na.rm = TRUE),
  proportion_concordant = mean(direction_concordant, na.rm = TRUE),
  median_recent_deltaAF = median(recent_deltaAF, na.rm = TRUE)
), by = .(pop, region)]

CORRELATIONS <- RECENT[
  is.finite(recent_deltaAF) & is.finite(established_median_deltaAF),
  {
    if (.N >= 3L && sd(recent_deltaAF) > 0 && sd(established_median_deltaAF) > 0) {
      test <- cor.test(
        recent_deltaAF,
        established_median_deltaAF,
        method = "spearman",
        exact = FALSE
      )
      .(n_genes = .N, spearman_rho = unname(test$estimate), p_value = test$p.value)
    } else {
      .(n_genes = .N, spearman_rho = NA_real_, p_value = NA_real_)
    }
  },
  by = .(pop, region)
]

# -----------------------------------------------------------------------------
# Plot: recent populations at the independently selected Figure 5 SNPs
# -----------------------------------------------------------------------------

# Preserve the candidate ordering used by the final Figure 5 script.
WIDE <- dcast(
  SUM_RAW[gene %chin% CANDIDATE_GENES],
  gene ~ region,
  value.var = "median_deltaAF"
)
WIDE[, same_direction := sign(AK) == sign(BC) & sign(AK) != 0]
WIDE[, pattern_order := fifelse(
  same_direction & AK > 0 & BC > 0, 1L,
  fifelse(same_direction & AK < 0 & BC < 0, 2L, 3L)
)]
WIDE[, shared_score := pmin(abs(AK), abs(BC))]
setorder(WIDE, pattern_order, -shared_score, gene)
GENE_ORDER <- WIDE$gene

RECENT[, gene := factor(gene, levels = GENE_ORDER)]
ESTABLISHED[, gene := factor(gene, levels = GENE_ORDER)]
RECENT[, gene_x := as.numeric(gene)]
ESTABLISHED[, gene_x := as.numeric(gene)]

PANEL_RECENT <- ggplot() +
  geom_hline(
    yintercept = seq(0, 1, 0.25),
    colour = "grey89",
    linewidth = 0.4
  ) +
  # Dotted: regional marine reference frequency.
  geom_segment(
    data = ESTABLISHED,
    aes(
      x = gene_x - 0.34,
      xend = gene_x + 0.34,
      y = marine_af,
      yend = marine_af,
      colour = region
    ),
    linetype = "dotted",
    linewidth = 0.9,
    alpha = 0.85
  ) +
  # Solid: median frequency among established freshwater populations.
  geom_segment(
    data = ESTABLISHED,
    aes(
      x = gene_x - 0.34,
      xend = gene_x + 0.34,
      y = established_median_af,
      yend = established_median_af,
      colour = region
    ),
    linewidth = 1.35,
    lineend = "round"
  ) +
  # Open symbols: the five recent populations.
  geom_point(
    data = RECENT,
    aes(x = gene_x, y = af, colour = region, shape = pop),
    size = 3.4,
    stroke = 1.05,
    fill = "white",
    position = position_dodge(width = 0.58)
  ) +
  scale_colour_manual(
    values = REGION_COLORS,
    breaks = c("AK", "BC"),
    labels = c(AK = "Alaska", BC = "British Columbia")
  ) +
  scale_shape_manual(values = POP_SHAPES) +
  scale_x_continuous(
    breaks = seq_along(GENE_ORDER),
    labels = GENE_ORDER,
    expand = expansion(mult = c(0.025, 0.025))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Focal-allele frequency",
    colour = "Region",
    shape = "Recent population",
    caption = paste0(
      "Open symbols show recently colonized populations at SNPs selected ",
      "independently from established populations. Dotted segments show the ",
      "regional marine frequency; solid segments show the established freshwater median."
    )
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(
      size = 12,
      colour = "black",
      face = "italic",
      angle = 45,
      hjust = 1
    ),
    axis.line = element_line(colour = "black", linewidth = 0.7),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.caption = element_text(size = 10.5, hjust = 0),
    plot.margin = margin(12, 16, 12, 12)
  )

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------

fwrite(
  RECENT[order(gene, region, pop)],
  file.path(OUTDIR, "Figure5_recent_fixedSNP_population_frequencies.tsv"),
  sep = "\t"
)
fwrite(
  COMPLETENESS[order(gene, region, pop)],
  file.path(OUTDIR, "Figure5_recent_fixedSNP_completeness.tsv"),
  sep = "\t"
)
fwrite(
  POP_SUMMARY[order(region, pop)],
  file.path(OUTDIR, "Figure5_recent_direction_summary.tsv"),
  sep = "\t"
)
fwrite(
  CORRELATIONS[order(region, pop)],
  file.path(OUTDIR, "Figure5_recent_established_correlations.tsv"),
  sep = "\t"
)

ggsave(
  file.path(OUTDIR, "Figure5_recent_populations_fixedSNP.pdf"),
  PANEL_RECENT,
  width = 13.5,
  height = 6.4,
  device = cairo_pdf
)
ggsave(
  file.path(OUTDIR, "Figure5_recent_populations_fixedSNP.png"),
  PANEL_RECENT,
  width = 13.5,
  height = 6.4,
  dpi = DPI,
  bg = "white"
)

cat("===== Figure 5 recent-population validation =====\n")
cat("Candidate genes fixed from established populations:", length(CANDIDATE_GENES), "\n")
cat("Candidate genes:", paste(CANDIDATE_GENES, collapse = ", "), "\n")
cat("Recent populations expected:", paste(RECENT_META$pop, collapse = ", "), "\n")
cat("Recent populations observed:", paste(sort(unique(RECENT$pop)), collapse = ", "), "\n")
cat("Missing gene-population observations:", COMPLETENESS[observed == FALSE, .N], "\n\n")
print(POP_SUMMARY[order(region, pop)])
cat("\n===== Correlation with established regional median delta AF =====\n")
print(CORRELATIONS[order(region, pop)])
cat("\n[output]", OUTDIR, "\n")
