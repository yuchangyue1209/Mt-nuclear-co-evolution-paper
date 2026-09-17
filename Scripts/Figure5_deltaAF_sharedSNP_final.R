#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ============================================================
# Figure 3A-B / manuscript Figure 5A-B
#
# IMPORTANT:
# These inputs already contain the original single representative SNP
# selected for each gene and used in both Alaska and British Columbia.
# This script does not select, replace, or reorient SNPs or focal alleles.
# ============================================================

BASE <- "/mnt/spareHD_2/nu_287/q2_parallelism"
INDIR <- file.path(BASE, "q2_deltaAF_withAMO_noLB_Jesse")
OUTDIR <- file.path(INDIR, "Figure5_final_colored_panelB")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

IN_TOP <- file.path(
  INDIR,
  "Fig3B_sharedRepresentativeSNP.withAMO_noLB_Jesse.tsv"
)
IN_LONG <- file.path(
  INDIR,
  "Fig3B_sharedSNP_populationData.withAMO_noLB_Jesse.tsv.gz"
)
IN_SUM <- file.path(
  INDIR,
  "Fig3B_sharedSNP_regionSummary.withAMO_noLB_Jesse.tsv"
)

THRESHOLD <- 0.50
DPI <- 500

REGION_LEVELS <- c("AK", "BC")
REGION_LABELS <- c(AK = "Alaska", BC = "British Columbia")
REGION_SHAPES <- c(AK = 16, BC = 17)
REGION_COLORS <- c(AK = "#1796C4", BC = "#08A77B")

required_files <- c(IN_TOP, IN_LONG, IN_SUM)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

# ============================================================
# Read the original shared-SNP results
# ============================================================

TOP <- fread(IN_TOP)
LONG <- fread(IN_LONG)
SUM_RAW <- fread(IN_SUM)

TOP[, gene := tolower(trimws(gene))]
LONG[, gene := tolower(trimws(gene))]
SUM_RAW[, gene := tolower(trimws(gene))]

LONG[, region := factor(toupper(trimws(region)), levels = REGION_LEVELS)]
SUM_RAW[, region := factor(toupper(trimws(region)), levels = REGION_LEVELS)]

TOP_META <- TOP[, .(
  gene,
  selected_snp = snp,
  shared_score,
  same_nonzero_direction,
  selection_class
)]

SUM <- SUM_RAW[, .(
  gene,
  region,
  n_pop,
  median_deltaAF,
  abs_median_deltaAF,
  median_freshwater_af,
  marine_af,
  marine_pop,
  focal_allele,
  shared_snp
)]

SUM <- merge(SUM, TOP_META, by = "gene", all.x = TRUE)

# At least one region must have |median delta AF| >= 0.5.
CANDIDATE_GENES <- sort(unique(
  SUM[abs_median_deltaAF >= THRESHOLD, gene]
))
if (!length(CANDIDATE_GENES)) stop("No candidate genes passed the threshold.")

# ============================================================
# Theme
# ============================================================

theme_fig <- theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(colour = "black", size = 16),
    axis.text = element_text(colour = "black", size = 13),
    axis.line = element_line(colour = "black", linewidth = 0.7),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.16, "cm"),
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.key.size = unit(0.7, "cm"),
    plot.caption = element_text(size = 10.5, hjust = 0, colour = "black"),
    plot.margin = margin(10, 14, 10, 12)
  )

# ============================================================
# Panel A: original median-delta-AF correlation
# ============================================================

A_DATA <- dcast(
  SUM,
  gene + selected_snp + selection_class + shared_score ~ region,
  value.var = "median_deltaAF"
)
setnames(A_DATA, c("AK", "BC"),
         c("AK_median_deltaAF", "BC_median_deltaAF"))

A_DATA[, candidate := gene %chin% CANDIDATE_GENES]
A_DATA <- A_DATA[
  is.finite(AK_median_deltaAF) & is.finite(BC_median_deltaAF)
]

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

COR_STATS <- data.table(
  test = c("Pearson", "Spearman"),
  estimate = c(unname(PEARSON$estimate), unname(SPEARMAN$estimate)),
  p_value = c(PEARSON$p.value, SPEARMAN$p.value),
  n_gene = nrow(A_DATA),
  SNP_strategy = paste0(
    "same SNP in AK and BC; same-direction SNPs preferred; ",
    "maximize min(|AK median deltaAF|, |BC median deltaAF|)"
  )
)

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
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.45, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed",
             linewidth = 0.45, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              linewidth = 0.7, colour = "grey45") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey78", linewidth = 0.9) +
  geom_point(
    data = A_DATA[candidate == FALSE],
    shape = 21, fill = "grey80", colour = "grey35",
    size = 3.7, stroke = 0.55
  ) +
  geom_point(
    data = A_DATA[candidate == TRUE],
    shape = 21, fill = "#F4A261", colour = "black",
    size = 4.5, stroke = 0.9
  ) +
  geom_text_repel(
    data = A_DATA[candidate == TRUE],
    aes(label = gene),
    fontface = "italic", size = 4.2,
    max.overlaps = Inf, min.segment.length = 0,
    segment.size = 0.4, segment.colour = "grey35",
    box.padding = 0.7, point.padding = 0.45,
    force = 2, max.time = 5, max.iter = 20000, seed = 100
  ) +
  annotate(
    "text",
    x = lims_A[1] + 0.03 * diff(lims_A),
    y = lims_A[2] - 0.03 * diff(lims_A),
    label = stat_label,
    hjust = 0, vjust = 1, size = 4.5
  ) +
  coord_equal(xlim = lims_A, ylim = lims_A, expand = FALSE, clip = "off") +
  labs(
    title = "A",
    x = "Alaska median \u0394AF",
    y = "British Columbia median \u0394AF"
  ) +
  theme_fig

# ============================================================
# Panel B data and ordering
# ============================================================

B_SUM <- SUM[gene %chin% CANDIDATE_GENES]
B_LONG <- LONG[gene %chin% CANDIDATE_GENES]

B_WIDE <- dcast(
  B_SUM,
  gene + selection_class + shared_score + selected_snp ~ region,
  value.var = "median_deltaAF"
)
B_WIDE[, same_direction := sign(AK) == sign(BC) & sign(AK) != 0]
B_WIDE[, pattern_order := fifelse(
  same_direction & AK > 0 & BC > 0, 1L,
  fifelse(same_direction & AK < 0 & BC < 0, 2L, 3L)
)]
setorder(B_WIDE, pattern_order, -shared_score, gene)
GENE_ORDER <- B_WIDE$gene

B_SUM[, gene_order := match(gene, GENE_ORDER)]
B_LONG[, gene_order := match(gene, GENE_ORDER)]

GENE_SPACING <- 2.50
REGION_OFFSET <- c(AK = -0.40, BC = 0.40)

B_SUM[, gene_center := (gene_order - 1) * GENE_SPACING + 1]
B_LONG[, gene_center := (gene_order - 1) * GENE_SPACING + 1]
B_SUM[, x_position := gene_center + REGION_OFFSET[as.character(region)]]
B_LONG[, x_position := gene_center + REGION_OFFSET[as.character(region)]]

# Deterministic jitter prevents points from moving between runs.
setorder(B_LONG, gene_order, region, pop)
B_LONG[, jitter_offset := if (.N == 1L) 0 else
         seq(-0.21, 0.21, length.out = .N), by = .(gene, region)]
B_LONG[, x_jitter := x_position + jitter_offset]

# Signed median delta AF values are retained.
B_SUM[, delta_label := sprintf("%+.2f", median_deltaAF)]
B_SUM[, label_y := fifelse(median_deltaAF >= 0, 1.035, -0.035)]
B_SUM[, label_vjust := fifelse(median_deltaAF >= 0, 0, 1)]

GENE_LABELS <- unique(B_SUM[, .(
  gene, gene_order, gene_center, selection_class
)])

# Reserve * for the region-specific large-effect threshold.
# Use a dagger for a fallback SNP.
GENE_LABELS[, gene_label := fifelse(
  selection_class == "opposite_or_zero_fallback",
  paste0(gene, "\u2020"),
  gene
)]
setorder(GENE_LABELS, gene_order)

SEPARATORS <- if (nrow(GENE_LABELS) > 1L) {
  (GENE_LABELS$gene_center[-nrow(GENE_LABELS)] +
     GENE_LABELS$gene_center[-1L]) / 2
} else numeric(0)

# Add * to AK or BC when that region has |median delta AF| >= 0.5.
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
  "Signed numbers report regional median \u0394AF; ",
  "* indicates |median \u0394AF| \u2265 0.5 in that region.",
  if (HAS_FALLBACK) paste0(
    " \u2020 No same-direction SNP with complete AK and BC coverage was ",
    "available; the highest shared-score fallback SNP is shown."
  ) else ""
)

# ============================================================
# Panel B
# ============================================================

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
  # Dotted line: corresponding regional marine reference.
  geom_segment(
    data = B_SUM,
    aes(
      x = x_position - 0.29,
      xend = x_position + 0.29,
      y = marine_af,
      yend = marine_af
    ),
    linetype = "dotted",
    linewidth = 1.05,
    colour = "grey25"
  ) +
  # AK = blue circles; BC = green triangles.
  geom_point(
    data = B_LONG,
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
  # Solid line: median freshwater focal-allele frequency.
  geom_segment(
    data = B_SUM,
    aes(
      x = x_position - 0.31,
      xend = x_position + 0.31,
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
    aes(x = x_position, y = -0.125, label = region_label),
    size = 3.8,
    fontface = "bold",
    colour = "grey20"
  ) +
  geom_text(
    data = GENE_LABELS,
    aes(x = gene_center, y = -0.225, label = gene_label),
    fontface = "italic",
    size = 4.0,
    colour = "black"
  ) +
  scale_shape_manual(values = REGION_SHAPES, labels = REGION_LABELS) +
  scale_colour_manual(values = REGION_COLORS, labels = REGION_LABELS) +
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
    plot.margin = margin(12, 16, 58, 12)
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4, alpha = 0.9)
    ),
    shape = guide_legend(
      override.aes = list(size = 4, alpha = 0.9)
    )
  )

# ============================================================
# Save plot data, individual panels, and combined figure
# ============================================================

fwrite(COR_STATS,
       file.path(OUTDIR, "Figure5A_correlation_statistics.tsv"), sep = "\t")
fwrite(A_DATA,
       file.path(OUTDIR, "Figure5A_plot_data.tsv"), sep = "\t")
fwrite(B_SUM,
       file.path(OUTDIR, "Figure5B_plot_summary.tsv"), sep = "\t")
fwrite(B_LONG,
       file.path(OUTDIR, "Figure5B_population_data.tsv.gz"),
       sep = "\t", compress = "gzip")

ggsave(
  file.path(OUTDIR, "Figure5A_sharedSNP.pdf"),
  PANEL_A,
  width = 7.0, height = 6.3, units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(OUTDIR, "Figure5A_sharedSNP.png"),
  PANEL_A,
  width = 7.0, height = 6.3, units = "in",
  dpi = DPI, bg = "white"
)

ggsave(
  file.path(OUTDIR, "Figure5B_sharedSNP_colored.pdf"),
  PANEL_B,
  width = 15.5, height = 6.8, units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)
ggsave(
  file.path(OUTDIR, "Figure5B_sharedSNP_colored.png"),
  PANEL_B,
  width = 15.5, height = 6.8, units = "in",
  dpi = DPI, bg = "white",
  limitsize = FALSE
)

FIGURE5_AB <- PANEL_A / PANEL_B +
  plot_layout(heights = c(0.95, 1.05))

ggsave(
  file.path(OUTDIR, "Figure5_AB_sharedSNP_final.pdf"),
  FIGURE5_AB,
  width = 15.8, height = 12.2, units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)
ggsave(
  file.path(OUTDIR, "Figure5_AB_sharedSNP_final.png"),
  FIGURE5_AB,
  width = 15.8, height = 12.2, units = "in",
  dpi = DPI, bg = "white",
  limitsize = FALSE
)

cat("\n===== Figure 5 audit =====\n")
cat("Genes in Panel A:", nrow(A_DATA), "\n")
cat("Candidate genes in Panel B:", length(CANDIDATE_GENES), "\n")
cat("Pearson r:", unname(PEARSON$estimate), "\n")
cat("Pearson P:", PEARSON$p.value, "\n")
cat("Candidate genes:\n")
print(CANDIDATE_GENES)
cat("Output directory:", OUTDIR, "\n")
