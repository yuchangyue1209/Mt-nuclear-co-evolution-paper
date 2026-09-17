#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ============================================================
# 1. Paths
# ============================================================

BASE_DIR <- paste0(
  "/path/to/workspace/gene_fst_work_withNorway/",
  "pbs_mt_vs_nuclear_NorwayOutgroup"
)

RAW_PBS_FILE <- file.path(
  BASE_DIR,
  "OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_noAMO.tsv"
)

CONTROLLED_FILE <- file.path(
  BASE_DIR,
  "OXPHOS72_chr21_controlled",
  "OXPHOS72_chr21_controlled_mitonuclear_PBS.tsv"
)

CHR21_FILE <- file.path(
  BASE_DIR,
  "chr21_neutral_PBS_by_population.tsv"
)

ANNOT_FILE <- "/path/to/workspace/codeml_sites_summary.merged.tsv"

OUTDIR <- file.path(
  BASE_DIR,
  "OXPHOS72_chr21_controlled",
  "Figure5_originalAB_controlledC"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Check files
# ============================================================

input_files <- c(
  RAW_PBS_FILE,
  CONTROLLED_FILE,
  CHR21_FILE,
  ANNOT_FILE
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# ============================================================
# 3. Plot settings
# ============================================================

BASE_SIZE <- 16
TITLE_SIZE <- 20
AXIS_SIZE <- 18
TEXT_SIZE <- 14
LABEL_SIZE <- 4.3

# Original Panel A colors
status_colors <- c(
  "AK-enriched" = "#3B6FB6",
  "BC-enriched" = "#2CA25F",
  "Shared" = "#8E44AD",
  "Background" = "grey80"
)

# Complex colors
complex_colors <- c(
  "CI" = "#F8766D",
  "CII" = "#C49A00",
  "CIII" = "#00BA38",
  "CIV" = "#00BFC4",
  "CV" = "#619CFF",
  "Cytochrome C" = "#F564E3",
  "Other" = "grey60"
)

# Panel C region colors matching Panel A
region_colors <- c(
  "AK" = "#3B6FB6",
  "BC" = "#2CA25F"
)

theme_xy <- theme_classic(
  base_size = BASE_SIZE
) +
  theme(
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    ),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "grey95",
      color = "black",
      linewidth = 0.45
    ),
    strip.text = element_text(
      size = 14
    ),
    plot.title = element_text(
      face = "bold",
      size = TITLE_SIZE,
      hjust = 0
    ),
    axis.title = element_text(
      size = AXIS_SIZE
    ),
    axis.text = element_text(
      size = TEXT_SIZE,
      color = "black"
    ),
    legend.title = element_text(
      size = 14
    ),
    legend.text = element_text(
      size = 13
    )
  )

italic_label <- function(x) {
  paste0(
    "italic('",
    x,
    "')"
  )
}

normalize_complex <- function(x) {

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  y <- sub(
    "^COMPLEX[ _-]*",
    "",
    y
  )

  y[y == "I"] <- "CI"
  y[y == "II"] <- "CII"
  y[y == "III"] <- "CIII"
  y[y == "IV"] <- "CIV"
  y[y == "V"] <- "CV"

  y[
    y %in% c(
      "CYTOCHROME C",
      "CYTOCHROME_C",
      "CYTC"
    )
  ] <- "Cytochrome C"

  y[
    !y %in% c(
      "CI",
      "CII",
      "CIII",
      "CIV",
      "CV",
      "Cytochrome C"
    )
  ] <- "Other"

  y
}

first_nonmissing <- function(x) {

  valid <- x[
    !is.na(x) &
    x != ""
  ]

  if (length(valid) == 0) {
    return(NA_character_)
  }

  as.character(valid[1])
}

# ============================================================
# 4. Read data
# ============================================================

raw_dt <- fread(
  RAW_PBS_FILE
)

controlled_dt <- fread(
  CONTROLLED_FILE
)

chr21_dt <- fread(
  CHR21_FILE
)

annot_dt <- fread(
  ANNOT_FILE
)

# ============================================================
# 5. Check columns
# ============================================================

required_raw <- c(
  "gene",
  "focal",
  "region",
  "nu_PBS",
  "mt_PBS"
)

required_controlled <- c(
  "gene",
  "beta_mt",
  "partial_R2_mt",
  "cor_raw",
  "cor_residual",
  "p_perm",
  "candidate"
)

required_chr21 <- c(
  "focal",
  "region",
  "chr21_PBS"
)

if (!all(required_raw %in% names(raw_dt))) {
  stop(
    "Raw PBS file is missing required columns"
  )
}

if (!all(required_controlled %in% names(controlled_dt))) {
  stop(
    "Controlled result file is missing required columns"
  )
}

if (!all(required_chr21 %in% names(chr21_dt))) {
  stop(
    "chr21 PBS file is missing required columns"
  )
}

if (!all(c("gene", "role", "complex") %in% names(annot_dt))) {
  stop(
    "Annotation file must contain gene, role, and complex"
  )
}

# ============================================================
# 6. Prepare annotation
# ============================================================

annot_gene <- annot_dt[
  ,
  .(
    role = first_nonmissing(role),
    complex = first_nonmissing(complex)
  ),
  by = gene
]

annot_gene[
  ,
  complex_plot := normalize_complex(
    complex
  )
]

# ============================================================
# 7. Prepare raw OXPHOS data
# ============================================================

# Use the same populations as the controlled analysis
raw_dt <- raw_dt[
  !focal %in% c(
    "AMO",
    "LB"
  )
]

raw_dt <- merge(
  raw_dt,
  annot_gene,
  by = "gene",
  all.x = TRUE
)

# Keep structural OXPHOS subunits
raw_sub <- raw_dt[
  is.na(role) |
  role == "subunit"
]

raw_sub[
  is.na(complex_plot),
  complex_plot := "Other"
]

raw_sub[
  ,
  complex_plot := factor(
    complex_plot,
    levels = c(
      "CI",
      "CII",
      "CIII",
      "CIV",
      "CV",
      "Cytochrome C",
      "Other"
    )
  )
]

# Average duplicated gene-population records
raw_sub <- raw_sub[
  ,
  .(
    nu_PBS = mean(
      nu_PBS,
      na.rm = TRUE
    ),
    mt_PBS = mean(
      mt_PBS,
      na.rm = TRUE
    ),
    role = first_nonmissing(role),
    complex = first_nonmissing(complex),
    complex_plot = first_nonmissing(
      as.character(complex_plot)
    )
  ),
  by = .(
    gene,
    focal,
    region
  )
]

raw_sub[
  ,
  complex_plot := factor(
    complex_plot,
    levels = c(
      "CI",
      "CII",
      "CIII",
      "CIV",
      "CV",
      "Cytochrome C",
      "Other"
    )
  )
]

# ============================================================
# 8. Panel A
# Original AK versus BC mean nuclear PBS
# ============================================================

gene_summary <- raw_sub[
  is.finite(nu_PBS),
  .(
    mean_PBS = mean(
      nu_PBS,
      na.rm = TRUE
    ),
    n_pop = .N
  ),
  by = .(
    gene,
    region
  )
]

ak_dt <- gene_summary[
  region == "AK",
  .(
    gene,
    mean_PBS_AK = mean_PBS,
    n_pop_AK = n_pop
  )
]

bc_dt <- gene_summary[
  region == "BC",
  .(
    gene,
    mean_PBS_BC = mean_PBS,
    n_pop_BC = n_pop
  )
]

compare_dt <- merge(
  ak_dt,
  bc_dt,
  by = "gene"
)

compare_dt <- merge(
  compare_dt,
  annot_gene,
  by = "gene",
  all.x = TRUE
)

compare_dt <- compare_dt[
  is.na(role) |
  role == "subunit"
]

# Rank genes within each region.
# Genes among the five highest PBS values in both regions are "Shared".
compare_dt[
  ,
  rank_AK := frank(
    -mean_PBS_AK,
    ties.method = "min"
  )
]

compare_dt[
  ,
  rank_BC := frank(
    -mean_PBS_BC,
    ties.method = "min"
  )
]

compare_dt[
  ,
  status := fcase(
    rank_AK <= 5 & rank_BC <= 5, "Shared",
    rank_AK <= 5, "AK-enriched",
    rank_BC <= 5, "BC-enriched",
    default = "Background"
  )
]

compare_dt[
  ,
  status := factor(
    status,
    levels = c(
      "AK-enriched",
      "BC-enriched",
      "Shared",
      "Background"
    )
  )
]

label_A <- compare_dt[
  status != "Background"
]

pA <- ggplot(
  compare_dt,
  aes(
    x = mean_PBS_BC,
    y = mean_PBS_AK
  )
) +
  geom_point(
    aes(
      color = status
    ),
    size = 3.0,
    alpha = 0.88
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = label_A,
    aes(
      label = italic_label(
        gene
      )
    ),
    parse = TRUE,
    size = LABEL_SIZE,
    color = "black",
    box.padding = 0.50,
    point.padding = 0.35,
    force = 1.5,
    max.time = 3,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = status_colors,
    drop = FALSE,
    name = NULL
  ) +
  labs(
    x = "Mean nuclear PBS in British Columbia",
    y = "Mean nuclear PBS in Alaska",
    title = "A"
  ) +
  theme_xy +
  theme(
    legend.position = "right"
  )

# ============================================================
# 9. Panel B
# Original raw nuclear PBS versus mtDNA PBS correlation
# ============================================================

cor_dt <- raw_sub[
  is.finite(nu_PBS) &
  is.finite(mt_PBS),
  {
    if (
      .N >= 5 &&
      var(nu_PBS) > 0 &&
      var(mt_PBS) > 0
    ) {

      correlation_test <- suppressWarnings(
        cor.test(
          nu_PBS,
          mt_PBS,
          method = "pearson"
        )
      )

      .(
        n_pop = .N,
        cor_r = unname(
          correlation_test$estimate
        ),
        p = correlation_test$p.value,
        mean_nu_PBS = mean(
          nu_PBS,
          na.rm = TRUE
        ),
        mean_mt_PBS = mean(
          mt_PBS,
          na.rm = TRUE
        )
      )

    } else {

      .(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(
          nu_PBS,
          na.rm = TRUE
        ),
        mean_mt_PBS = mean(
          mt_PBS,
          na.rm = TRUE
        )
      )
    }
  },
  by = gene
]

cor_dt[
  ,
  p_bh := p.adjust(
    p,
    method = "BH"
  )
]

cor_dt <- merge(
  cor_dt,
  annot_gene,
  by = "gene",
  all.x = TRUE
)

cor_dt[
  is.na(complex_plot),
  complex_plot := "Other"
]

cor_dt <- cor_dt[
  is.finite(cor_r) &
  is.finite(p) &
  p > 0
]

cor_dt[
  ,
  logp := -log10(
    p
  )
]

cor_dt[
  ,
  complex_plot := factor(
    complex_plot,
    levels = c(
      "CI",
      "CII",
      "CIII",
      "CIV",
      "CV",
      "Cytochrome C",
      "Other"
    )
  )
]

# Label significant genes and strongest correlations
label_B <- cor_dt[
  p < 0.05 |
  abs(cor_r) > 0.55
]

label_B <- label_B[order(p, -abs(cor_r))]

label_B <- label_B[
  1:min(
    .N,
    12
  )
]

pB <- ggplot(
  cor_dt,
  aes(
    x = cor_r,
    y = logp
  )
) +
  geom_point(
    aes(color = complex_plot),
    size = 3.0,
    alpha = 0.88
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "grey45",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = label_B,
    aes(
      label = italic_label(
        gene
      )
    ),
    parse = TRUE,
    size = LABEL_SIZE,
    color = "black",
    box.padding = 0.50,
    point.padding = 0.35,
    force = 1.4,
    max.time = 3,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = complex_colors,
    drop = TRUE,
    name = "Complex"
  ) +
  labs(
    x = "Correlation between nuclear PBS and mtDNA PBS",
    y = expression(
      -log[10](
        P
      )
    ),
    title = "B"
  ) +
  theme_xy +
  theme(
    legend.position = "right"
  )

# ============================================================
# 10. Prepare controlled results for Panel C
# ============================================================

controlled_dt <- merge(
  controlled_dt,
  annot_gene,
  by = "gene",
  all.x = TRUE
)

controlled_dt[
  is.na(complex_plot),
  complex_plot := "Other"
]

controlled_dt[
  ,
  candidate := as.logical(
    candidate
  )
]

mixed_complexes <- c(
  "CI",
  "CIII",
  "CIV",
  "CV"
)

# Select the positive gene with the smallest permutation P
# from each mixed-genome complex
top_by_complex <- controlled_dt[
  complex_plot %in% mixed_complexes &
  is.finite(beta_mt) &
  is.finite(p_perm) &
  beta_mt > 0,
  .SD[
    which.min(
      p_perm
    )
  ],
  by = complex_plot
]

top_by_complex[
  ,
  complex_plot := factor(
    complex_plot,
    levels = mixed_complexes
  )
]

setorder(
  top_by_complex,
  complex_plot
)

cat(
  "\nTop controlled gene from each complex:\n"
)

print(
  top_by_complex[
    ,
    .(
      complex_plot,
      gene,
      beta_mt,
      partial_R2_mt,
      cor_residual,
      p_perm
    )
  ]
)

# ============================================================
# 11. Merge chr21 PBS for Panel C
# ============================================================

chr21_use <- unique(
  chr21_dt[
    !focal %in% c(
      "AMO",
      "LB"
    ),
    .(
      focal,
      region,
      chr21_PBS
    )
  ]
)

plotC_input <- raw_sub[
  gene %in% top_by_complex$gene
]

plotC_input <- merge(
  plotC_input,
  chr21_use,
  by = c(
    "focal",
    "region"
  ),
  all.x = TRUE
)

plotC_input <- plotC_input[
  is.finite(nu_PBS) &
  is.finite(mt_PBS) &
  is.finite(chr21_PBS)
]

plotC_input[
  ,
  region := factor(
    region,
    levels = c(
      "AK",
      "BC"
    )
  )
]

# ============================================================
# 12. Calculate controlled residuals
# ============================================================

calculate_residuals <- function(d) {

  d <- copy(
    d
  )

  d[
    ,
    region := droplevels(
      factor(
        region
      )
    )
  ]

  nu_model <- lm(
    nu_PBS ~ chr21_PBS + region,
    data = d
  )

  mt_model <- lm(
    mt_PBS ~ chr21_PBS + region,
    data = d
  )

  d[
    ,
    nu_residual := residuals(
      nu_model
    )
  ]

  d[
    ,
    mt_residual := residuals(
      mt_model
    )
  ]

  d
}

plotC_dt <- plotC_input[
  ,
  calculate_residuals(
    .SD
  ),
  by = gene
]


if ("complex_plot" %in% names(plotC_dt)) {
  plotC_dt[, complex_plot := NULL]
}

plotC_dt <- merge(
  plotC_dt,
  top_by_complex[
    ,
    .(
      gene,
      complex_plot,

      beta_mt,
      partial_R2_mt,
      cor_residual,
      p_perm
    )
  ],
  by = "gene",
  all.x = TRUE
)

plotC_dt[
  ,
  complex_plot := factor(
    complex_plot,
    levels = mixed_complexes
  )
]

# ============================================================
# 13. Parsed facet labels with italic gene names
# ============================================================

plotC_dt[
  ,
  facet_label := paste0(
    "bold('",
    as.character(
      complex_plot
    ),
    "')*':'~italic('",
    gene,
    "')"
  )
]

facet_order <- top_by_complex[
  ,
  paste0(
    "bold('",
    as.character(
      complex_plot
    ),
    "')*':'~italic('",
    gene,
    "')"
  )
]

plotC_dt[
  ,
  facet_label := factor(
    facet_label,
    levels = facet_order
  )
]

# ============================================================
# 14. Panel C statistics and label positions
# ============================================================

panelC_stats <- unique(
  plotC_dt[
    ,
    .(
      facet_label,
      gene,
      complex_plot,
      cor_residual,
      partial_R2_mt,
      p_perm
    )
  ]
)

# Use common positions because facets have fixed scales
global_x_range <- range(
  plotC_dt$mt_residual,
  finite = TRUE
)

global_y_range <- range(
  plotC_dt$nu_residual,
  finite = TRUE
)

panelC_stats[
  ,
  x_position := global_x_range[1] +
    0.04 * diff(
      global_x_range
    )
]

panelC_stats[
  ,
  y_position := global_y_range[2] -
    0.05 * diff(
      global_y_range
    )
]

panelC_stats[
  ,
  stat_label := paste0(
    "Residual r = ",
    sprintf(
      "%.2f",
      cor_residual
    ),
    "\nPermutation P = ",
    formatC(
      p_perm,
      format = "f",
      digits = 3
    )
  )
]

# ============================================================
# 15. Panel C
# Controlled residual relationships
# ============================================================

pC <- ggplot(
  plotC_dt,
  aes(
    x = mt_residual,
    y = nu_residual
  )
) +
  geom_hline(
    yintercept = 0,
    color = "grey82",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey82",
    linewidth = 0.4
  ) +
  geom_point(
    aes(
      color = region
    ),
    size = 2.8,
    alpha = 0.90
  ) +
  geom_smooth(
    data = plotC_dt,
    mapping = aes(
      x = mt_residual,
      y = nu_residual,
      group = facet_label
    ),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    level = 0.95,
    color = "black",
    fill = "grey75",
    alpha = 0.45,
    linewidth = 0.8,
    inherit.aes = FALSE,
    na.rm = TRUE
  ) +
  geom_text(
    data = panelC_stats,
    aes(
      x = x_position,
      y = y_position,
      label = stat_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  facet_wrap(
    ~ facet_label,
    scales = "fixed",
    nrow = 1,
    labeller = label_parsed
  ) +
  scale_color_manual(
    values = region_colors,
    name = "Region"
  ) +
  labs(
    x = paste0(
      "Residual mitochondrial PBS after controlling ",
      "for chr21 PBS and region"
    ),
    y = paste0(
      "Residual nuclear gene PBS after controlling ",
      "for chr21 PBS and region"
    ),
    title = "C"
  ) +
  theme_xy +
  theme(
    legend.position = "top",
    strip.text = element_text(
      size = 14
    )
  )

# ============================================================
# 16. Combine panels
# ============================================================

figure5 <- (
  pA +
  pB
) /
  pC +
  plot_layout(
    heights = c(
      1,
      0.92
    )
  ) &
  theme(
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 25
    )
  )

# ============================================================
# 17. Save individual panels
# ============================================================

ggsave(
  file.path(
    OUTDIR,
    "Figure5A_AK_vs_BC_mean_nuclear_PBS.png"
  ),
  pA,
  width = 7.5,
  height = 6.2,
  dpi = 600
)

ggsave(
  file.path(
    OUTDIR,
    "Figure5A_AK_vs_BC_mean_nuclear_PBS.pdf"
  ),
  pA,
  width = 7.5,
  height = 6.2,
  device = cairo_pdf
)

ggsave(
  file.path(
    OUTDIR,
    "Figure5B_raw_nuclear_mtPBS_correlation.png"
  ),
  pB,
  width = 7.5,
  height = 6.2,
  dpi = 600
)

ggsave(
  file.path(
    OUTDIR,
    "Figure5B_raw_nuclear_mtPBS_correlation.pdf"
  ),
  pB,
  width = 7.5,
  height = 6.2,
  device = cairo_pdf
)

ggsave(
  file.path(
    OUTDIR,
    "Figure5C_chr21_controlled_top_genes.png"
  ),
  pC,
  width = 14,
  height = 5.0,
  dpi = 600
)

ggsave(
  file.path(
    OUTDIR,
    "Figure5C_chr21_controlled_top_genes.pdf"
  ),
  pC,
  width = 14,
  height = 5.0,
  device = cairo_pdf
)

# ============================================================
# 18. Save combined Figure 5
# ============================================================

ggsave(
  file.path(
    OUTDIR,
    "Figure5_originalAB_controlledC.png"
  ),
  figure5,
  width = 15,
  height = 11.5,
  dpi = 600
)

ggsave(
  file.path(
    OUTDIR,
    "Figure5_originalAB_controlledC.pdf"
  ),
  figure5,
  width = 15,
  height = 11.5,
  device = cairo_pdf
)

# ============================================================
# 19. Save source data
# ============================================================

fwrite(
  compare_dt,
  file.path(
    OUTDIR,
    "Figure5A_source_data.tsv"
  ),
  sep = "\t"
)

fwrite(
  cor_dt,
  file.path(
    OUTDIR,
    "Figure5B_source_data.tsv"
  ),
  sep = "\t"
)

fwrite(
  top_by_complex,
  file.path(
    OUTDIR,
    "Figure5C_top_gene_by_complex.tsv"
  ),
  sep = "\t"
)

fwrite(
  plotC_dt,
  file.path(
    OUTDIR,
    "Figure5C_residual_source_data.tsv"
  ),
  sep = "\t"
)

fwrite(
  panelC_stats,
  file.path(
    OUTDIR,
    "Figure5C_statistics.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 20. Print summary
# ============================================================

cat(
  "\nFigure 5 completed successfully.\n"
)

cat(
  "\nOutput directory:\n",
  OUTDIR,
  "\n"
)

cat(
  "\nPanel A labeled genes:\n"
)

print(
  label_A[
    ,
    .(
      gene,
      mean_PBS_AK,
      mean_PBS_BC,
      status
    )
  ]
)

cat(
  "\nPanel B raw significant genes:\n"
)

print(
  cor_dt[
    p < 0.05,
    .(
      gene,
      complex_plot,
      cor_r,
      p
    )
  ][
    order(
      p
    )
  ]
)

cat(
  "\nPanel C top controlled genes:\n"
)

print(
  top_by_complex[
    ,
    .(
      complex_plot,
      gene,
      beta_mt,
      partial_R2_mt,
      cor_residual,
      p_perm
    )
  ]
)

cat(
  "\nWritten files:\n"
)

print(
  list.files(
    OUTDIR,
    full.names = TRUE
  )
)
