#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript Figure5_replot_colors_largefont.R INPUT_RESULTS_DIR OUTPUT_DIR")
INPUT <- normalizePath(args[1], mustWork = TRUE)
OUTDIR <- args[2]
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
read_result <- function(name) {
  p <- file.path(INPUT, name)
  if (!file.exists(p)) stop("Missing input: ", p)
  fread(p)
}
A_DATA <- read_result("Figure5A_plot_data_all_populations.tsv")
SUMMARY <- read_result("Figure5B_plot_summary_all_populations.tsv")
LONG <- read_result("Figure5B_plot_population_data_all_populations.tsv.gz")
COR <- read_result("Figure5A_correlation_statistics_all_populations.tsv")
PEARSON <- list(estimate = COR[test == "Pearson", estimate], p.value = COR[test == "Pearson", p_value])
CANDIDATE_GENES <- A_DATA[candidate == TRUE, gene]
if (!length(CANDIDATE_GENES)) stop("No candidates in saved plot data")
THRESHOLD <- 0.5
parameter_path <- file.path(INPUT, "Figure5_run_parameters.tsv")
if (file.exists(parameter_path)) {
  pars <- fread(parameter_path)
  v <- pars[parameter == "candidate_threshold", value]
  if (length(v) == 1) THRESHOLD <- as.numeric(v)
} else {
  message("No saved run parameters; using the original default candidate threshold 0.5 for region stars.")
}
DPI <- 300
MARINE_COLOR <- "#7B3294"
theme_fig <- theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(colour = "black", size = 22),
    axis.text = element_text(colour = "black", size = 18),
    axis.line = element_line(colour = "black", linewidth = 0.7),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = grid::unit(0.16, "cm"),
    plot.title = element_text(face = "bold", size = 26, hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.key.size = grid::unit(0.7, "cm"),
    plot.caption = element_text(size = 13, hjust = 0, colour = "black"),
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
    fill = "#D989AB",
    colour = "black",
    size = 4.5,
    stroke = 0.9
  ) +
  geom_text_repel(
    data = A_DATA[candidate == TRUE],
    aes(label = gene),
    fontface = "italic",
    size = 5.5,
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
    size = 5.5
  ) +
  coord_equal(
    xlim = lims_A,
    ylim = lims_A,
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    title = "A",
    x = expression(Alaska~median~Delta*AF),
    y = expression(British~Columbia~median~Delta*AF)
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
  "Blue circles denote established freshwater populations; orange circles denote recently colonized populations. ",
  "Purple and black segments denote regional marine frequencies and freshwater medians, respectively. ",
  "AK and BC identify the regional groups; * indicates |median deltaAF| >= 0.5.",
  if (HAS_FALLBACK) " † The highest shared-score fallback SNP is shown where no same-direction SNP with complete coverage was available." else ""
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
    data = B_LONG,
    aes(x = x_jitter, y = af, colour = population_history),
    shape = 16, size = 3.5, alpha = 0.85
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
    data = REGION_LABEL_DATA,
    aes(
      x = x_position,
      y = -0.125,
      label = region_label
    ),
    size = 5.0,
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
    size = 5.0,
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(established = "#1796C4", recent = "#F4A261"),
    breaks = c("established", "recent"),
    labels = c("Established freshwater", "Recently colonized"),
    drop = FALSE
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
    caption = paste(strwrap(PANEL_B_CAPTION, width = 115), collapse = "\n")
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


cat("[OK] Results saved to:", OUTDIR, "\n")
