#!/usr/bin/env Rscript

# Supplementary Figure: original Figure 5A plus Figure 5B with the five
# recently colonized populations added at the same preselected SNPs.
#
# Established populations retain the original symbols and colours.
# Recent AK populations are blue-filled circles with a red outline;
# recent BC populations are green-filled triangles with a red outline.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# Load the exact final Figure 5 analysis. This preserves Panel A, the 11 genes,
# representative SNPs, focal alleles, gene order, offsets, labels and statistics.
# -----------------------------------------------------------------------------

MAIN_SCRIPT <- "/work/cyu/Kuster2026_genomewide_reanalysis/Figure5_deltaAF_sharedSNP_final.R"

if (!file.exists(MAIN_SCRIPT)) {
  stop("Missing final Figure 5 script: ", MAIN_SCRIPT)
}

source(MAIN_SCRIPT, local = FALSE)

# -----------------------------------------------------------------------------
# Read recent-population frequencies previously validated against the fixed SNPs
# -----------------------------------------------------------------------------

RECENT_DIR <- file.path(INDIR, "Figure5_recent_populations_fixedSNP")
RECENT_FILE <- file.path(
  RECENT_DIR,
  "Figure5_recent_fixedSNP_population_frequencies.tsv"
)

if (!file.exists(RECENT_FILE)) {
  stop(
    "Missing recent-population table: ", RECENT_FILE, "\n",
    "Run Figure5_recent_populations_fixedSNP.R first."
  )
}

RECENT <- fread(RECENT_FILE)
RECENT[, gene := tolower(trimws(gene))]
RECENT[, region := toupper(trimws(region))]
RECENT[, pop := toupper(trimws(pop))]
RECENT[, af := as.numeric(af)]

EXPECTED_RECENT <- c("SC", "CH", "LB", "PACH", "FRED")
RECENT <- RECENT[
  pop %chin% EXPECTED_RECENT &
    gene %chin% CANDIDATE_GENES &
    is.finite(af)
]

if (!nrow(RECENT)) stop("No valid recent-population observations were found.")

# Enforce one row per fixed gene-SNP-population combination.
DUPLICATES <- RECENT[, .N, by = .(gene, selected_snp, pop)][N > 1L]
if (nrow(DUPLICATES)) {
  stop("Duplicate recent gene-SNP-population rows were found.")
}

# Confirm complete coverage: 11 genes x 5 recent populations.
EXPECTED_GRID <- CJ(
  gene = CANDIDATE_GENES,
  pop = EXPECTED_RECENT,
  unique = TRUE
)
OBSERVED_GRID <- unique(RECENT[, .(gene, pop)])
MISSING_GRID <- EXPECTED_GRID[!OBSERVED_GRID, on = .(gene, pop)]
if (nrow(MISSING_GRID)) {
  fwrite(
    MISSING_GRID,
    file.path(RECENT_DIR, "FigureS5_missing_recent_observations.tsv"),
    sep = "\t"
  )
  stop("Recent observations are incomplete; see FigureS5_missing_recent_observations.tsv")
}

# Use the exact horizontal coordinates from the final Figure 5B.
RECENT[, gene_order := match(gene, GENE_ORDER)]
RECENT[, gene_center := (gene_order - 1) * GENE_SPACING + 1]
RECENT[, x_position := gene_center + REGION_OFFSET[region]]

# Separate recent points slightly within each AK/BC position. Their offsets are
# deterministic so that the plot is reproducible.
setorder(RECENT, gene_order, region, pop)
RECENT[, recent_offset := if (.N == 1L) 0 else
         seq(-0.19, 0.19, length.out = .N),
       by = .(gene, region)]
RECENT[, x_recent := x_position + recent_offset]

AK_RECENT <- RECENT[region == "AK"]
BC_RECENT <- RECENT[region == "BC"]

# -----------------------------------------------------------------------------
# Supplementary Panel B
# -----------------------------------------------------------------------------

RECENT_RED <- "#D62728"

PANEL_B_SUPP <- PANEL_B +
  # AK recent: blue fill, red circle outline.
  geom_point(
    data = AK_RECENT,
    aes(x = x_recent, y = af),
    inherit.aes = FALSE,
    shape = 21,
    fill = REGION_COLORS[["AK"]],
    colour = RECENT_RED,
    size = 4.0,
    stroke = 1.15,
    alpha = 0.95
  ) +
  # BC recent: green fill, red triangle outline.
  geom_point(
    data = BC_RECENT,
    aes(x = x_recent, y = af),
    inherit.aes = FALSE,
    shape = 24,
    fill = REGION_COLORS[["BC"]],
    colour = RECENT_RED,
    size = 4.2,
    stroke = 1.15,
    alpha = 0.95
  ) +
  labs(
    caption = paste0(
      PANEL_B_CAPTION,
      " Blue circles and green triangles show Alaska and British Columbia, ",
      "respectively; red outlines identify recently colonized populations. ",
      "Marine reference frequencies are dotted and established-freshwater ",
      "medians are solid."
    )
  )

# -----------------------------------------------------------------------------
# Assemble and save the complete supplementary figure
# -----------------------------------------------------------------------------

FIGURE_SUPP <- PANEL_A / PANEL_B_SUPP +
  plot_layout(heights = c(0.92, 1.12))

SUPP_PDF <- file.path(
  RECENT_DIR,
  "FigureS5_AB_established_plus_recent.pdf"
)
SUPP_PNG <- file.path(
  RECENT_DIR,
  "FigureS5_AB_established_plus_recent.png"
)

ggsave(
  SUPP_PDF,
  FIGURE_SUPP,
  width = 16.2,
  height = 12.7,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  SUPP_PNG,
  FIGURE_SUPP,
  width = 16.2,
  height = 12.7,
  units = "in",
  dpi = DPI,
  bg = "white",
  limitsize = FALSE
)

fwrite(
  RECENT[order(gene_order, region, pop)],
  file.path(RECENT_DIR, "FigureS5_recent_plot_data.tsv"),
  sep = "\t"
)

cat("\n===== Supplementary Figure 5 audit =====\n")
cat("Panel A genes:", nrow(A_DATA), "\n")
cat("Panel B fixed candidate genes:", length(CANDIDATE_GENES), "\n")
cat("Recent populations:", paste(sort(unique(RECENT$pop)), collapse = ", "), "\n")
cat("Recent observations:", nrow(RECENT), "\n")
cat("Expected recent observations:", length(CANDIDATE_GENES) * length(EXPECTED_RECENT), "\n")
cat("PDF:", SUPP_PDF, "\n")
cat("PNG:", SUPP_PNG, "\n")
