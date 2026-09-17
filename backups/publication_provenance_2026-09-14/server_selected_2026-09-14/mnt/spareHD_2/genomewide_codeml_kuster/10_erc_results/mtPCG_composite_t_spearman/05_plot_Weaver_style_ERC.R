#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

erc_dir <- paste0(
  "/mnt/spareHD_2/genomewide_codeml_kuster/",
  "10_erc_results/mtPCG_composite_t_spearman"
)

figure_dir <- paste0(
  "/mnt/spareHD_2/genomewide_codeml_kuster/",
  "11_erc_figures/Figure_ERC_Weaver_style"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

summary_file <- file.path(
  erc_dir,
  "weaver_gene_set_erc_summary.tsv"
)

composite_file <- file.path(
  erc_dir,
  "weaver_gene_set_branch_composites.tsv"
)

loo_mt_file <- file.path(
  erc_dir,
  "CIV_leave_one_mt_gene_out.tsv"
)

loo_nuclear_file <- file.path(
  erc_dir,
  "CIV_leave_one_nuclear_gene_out.tsv"
)

loo_branch_file <- file.path(
  erc_dir,
  "CIV_leave_one_branch_out.tsv"
)

summary_data <- fread(summary_file)
composite_data <- fread(composite_file)
loo_mt <- fread(loo_mt_file)
loo_nuclear <- fread(loo_nuclear_file)
loo_branch <- fread(loo_branch_file)

# ----------------------------------------------------------
# Plot settings
# ----------------------------------------------------------

color_mito <- "#B2182B"
color_control <- "#2166AC"
color_civ <- "#D6604D"
color_terminal <- "#2166AC"
color_internal <- "#B2182B"
color_zero <- "grey45"

base_theme <- theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(
      size = 15,
      colour = "black"
    ),
    axis.text = element_text(
      size = 13,
      colour = "black"
    ),
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 15, 10, 10)
  )

# ----------------------------------------------------------
# Panel A: Functional gene sets
# ----------------------------------------------------------

functional_order <- c(
  "Cyto_ARS_control",
  "Cyto_ribosomal_control",
  "OXPHOS_assembly_factors",
  "Nmt_ARS",
  "Nmt_ribosomal",
  "Nmt_OXPHOS_structural",
  "Indirect_nmt",
  "Direct_nmt"
)

functional_labels <- c(
  Cyto_ARS_control = "Cyto-ARS",
  Cyto_ribosomal_control = "Cyto-ribo",
  OXPHOS_assembly_factors = "OXPHOS assembly",
  Nmt_ARS = "Nmt-ARS",
  Nmt_ribosomal = "Nmt-ribo",
  Nmt_OXPHOS_structural = "N-mt OXPHOS",
  Indirect_nmt = "Indirect n-mt",
  Direct_nmt = "Direct n-mt"
)

panel_a_data <- summary_data[
  nuclear_set %in% functional_order
]

panel_a_data[
  ,
  display_label := functional_labels[nuclear_set]
]

panel_a_data[
  ,
  display_label := factor(
    display_label,
    levels = functional_labels[functional_order]
  )
]

panel_a_data[
  ,
  significance := fifelse(
    is.finite(ci_lower) & ci_lower > 0,
    "Significant",
    "Not significant"
  )
]

panel_a <- ggplot(
  panel_a_data,
  aes(
    x = bootstrap_mean_rs,
    y = display_label,
    colour = significance
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = color_zero
  ) +
  geom_errorbar(
    aes(
      xmin = ci_lower,
      xmax = ci_upper
    ),
    orientation = "y",
    width = 0,
    linewidth = 1
  ) +
  geom_point(size = 3.8) +
  scale_colour_manual(
    values = c(
      "Significant" = color_mito,
      "Not significant" = color_control
    ),
    breaks = c(
      "Significant",
      "Not significant"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(-0.45, 0.60),
    breaks = seq(-0.4, 0.6, 0.2)
  ) +
  labs(
    title = "A",
    x = expression(
      "mt–nuclear ERC (" * r[s] * ")"
    ),
    y = NULL,
    colour = NULL
  ) +
  base_theme +
  theme(
    legend.position = "bottom"
  )

# ----------------------------------------------------------
# Panel B: Within-complex ERC
# ----------------------------------------------------------

complex_order <- c(
  "OXPHOS_CI",
  "OXPHOS_CII",
  "OXPHOS_CIII",
  "OXPHOS_CIV",
  "OXPHOS_CV"
)

complex_labels <- c(
  OXPHOS_CI = "Complex I",
  OXPHOS_CII = "Complex II†",
  OXPHOS_CIII = "Complex III",
  OXPHOS_CIV = "Complex IV",
  OXPHOS_CV = "Complex V"
)

panel_b_data <- summary_data[
  nuclear_set %in% complex_order
]

panel_b_data[
  ,
  display_label := complex_labels[nuclear_set]
]

panel_b_data[
  ,
  display_label := factor(
    display_label,
    levels = rev(complex_labels[complex_order])
  )
]

panel_b_data[
  ,
  significant_bootstrap :=
    ci_lower > 0
]

panel_b_data[
  ,
  significance := fifelse(
    is.finite(ci_lower) & ci_lower > 0,
    "Significant",
    "Not significant"
  )
]

panel_b_data[
  ,
  significance_label := fifelse(
    significant_bootstrap,
    "*",
    ""
  )
]

panel_b <- ggplot(
  panel_b_data,
  aes(
    x = bootstrap_mean_rs,
    y = display_label,
    colour = significance
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = color_zero
  ) +
  geom_errorbar(
    aes(
      xmin = ci_lower,
      xmax = ci_upper
    ),
    orientation = "y",
    width = 0,
    linewidth = 1
  ) +
  geom_point(size = 3.8) +
  geom_text(
    aes(
      x = ci_upper + 0.035,
      label = significance_label
    ),
    colour = "black",
    size = 6,
    vjust = 0.35
  ) +
  scale_colour_manual(
    values = c(
      "Significant" = color_mito,
      "Not significant" = color_control
    ),
    breaks = c(
      "Significant",
      "Not significant"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(-0.40, 0.62),
    breaks = seq(-0.4, 0.6, 0.2)
  ) +
  labs(
    title = "B",
    x = expression(
      "Within-complex ERC (" * r[s] * ")"
    ),
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "bottom"
  )

# ----------------------------------------------------------
# Panel C: CIV branch-level relationship
# ----------------------------------------------------------

panel_c_data <- composite_data[
  nuclear_set == "OXPHOS_CIV"
]

branch_is_terminal <- function(branch_label) {

  numbers <- as.integer(
    regmatches(
      branch_label,
      gregexpr("[0-9]+", branch_label)
    )[[1]]
  )

  if (length(numbers) < 2) {
    return(NA)
  }

  any(numbers <= 27)
}

panel_c_data[
  ,
  branch_type := fifelse(
    vapply(
      branch,
      branch_is_terminal,
      logical(1)
    ),
    "Terminal",
    "Internal"
  )
]

civ_rs <- cor(
  panel_c_data$nuclear_composite_rer,
  panel_c_data$mt_composite_rer,
  method = "spearman",
  use = "complete.obs"
)

panel_c <- ggplot(
  panel_c_data,
  aes(
    x = mt_composite_rer,
    y = nuclear_composite_rer,
    colour = branch_type
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    colour = "grey80"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    colour = "grey80"
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8,
    colour = "grey35",
    fill = "grey75",
    alpha = 0.35
  ) +
  geom_point(
    size = 2.8,
    alpha = 0.9
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.12,
    vjust = 1.4,
    size = 4.6,
    label = sprintf(
      "Spearman rₛ = %.3f",
      civ_rs
    ),
    parse = FALSE
  ) +
  scale_colour_manual(
    values = c(
      Terminal = color_terminal,
      Internal = color_internal
    )
  ) +
  labs(
    title = "C",
    x = "Mitochondrial CIV composite RER",
    y = "Nuclear CIV composite RER",
    colour = "Branch"
  ) +
  base_theme +
  theme(
    legend.position = c(0.22, 0.84),
    legend.background = element_blank()
  )

# ----------------------------------------------------------
# Panel D: CIV leave-one-out robustness
# ----------------------------------------------------------

full_rs <- civ_rs

loo_mt_plot <- loo_mt[
  ,
  .(
    analysis = "Leave one\nmt gene out",
    rs
  )
]

loo_nuclear_plot <- loo_nuclear[
  ,
  .(
    analysis = "Leave one\nnuclear gene out",
    rs
  )
]

loo_branch_plot <- loo_branch[
  ,
  .(
    analysis = "Leave one\nbranch out",
    rs
  )
]

panel_d_data <- rbindlist(
  list(
    loo_mt_plot,
    loo_nuclear_plot,
    loo_branch_plot
  )
)

panel_d_data[
  ,
  analysis := factor(
    analysis,
    levels = c(
      "Leave one\nmt gene out",
      "Leave one\nnuclear gene out",
      "Leave one\nbranch out"
    )
  )
]

panel_d <- ggplot(
  panel_d_data,
  aes(
    x = analysis,
    y = rs
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = color_zero,
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = full_rs,
    colour = color_civ,
    linewidth = 0.9
  ) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    fill = "grey92",
    colour = "grey35",
    linewidth = 0.7
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 2.2,
    alpha = 0.75,
    colour = "grey20"
  ) +
  annotate(
    "text",
    x = 3.45,
    y = full_rs,
    label = sprintf(
      "Full CIV rₛ = %.3f",
      full_rs
    ),
    parse = FALSE,
    hjust = 1,
    vjust = -0.7,
    size = 4.1,
    colour = color_civ
  ) +
  coord_cartesian(
    ylim = c(0, 0.46),
    clip = "off"
  ) +
  labs(
    title = "D",
    x = NULL,
    y = expression(
      "CIV ERC (" * r[s] * ")"
    )
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(
      size = 11.5
    )
  )

# ----------------------------------------------------------
# Combine and export
# ----------------------------------------------------------

combined_figure <- (
  panel_a | panel_b
) +
  plot_layout(
    guides = "collect",
    widths = c(1.15, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

pdf_file <- file.path(
  figure_dir,
  "Figure_ERC_Weaver_style_A_B_direct_indirect.pdf"
)

png_file <- file.path(
  figure_dir,
  "Figure_ERC_Weaver_style_A_B_direct_indirect.png"
)

ggsave(
  pdf_file,
  combined_figure,
  width = 12.5,
  height = 6.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  png_file,
  combined_figure,
  width = 12.5,
  height = 6.2,
  units = "in",
  dpi = 400,
  bg = "white"
)

# Save individual panels too
ggsave(
  file.path(figure_dir, "Panel_A_functional_sets.pdf"),
  panel_a,
  width = 7,
  height = 5.3,
  device = cairo_pdf
)

ggsave(
  file.path(figure_dir, "Panel_B_OXPHOS_complexes.pdf"),
  panel_b,
  width = 7,
  height = 5.3,
  device = cairo_pdf
)

ggsave(
  file.path(figure_dir, "Panel_C_CIV_scatter.pdf"),
  panel_c,
  width = 6.5,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(figure_dir, "Panel_D_CIV_robustness.pdf"),
  panel_d,
  width = 6.5,
  height = 5.5,
  device = cairo_pdf
)

cat("[OK] PDF:", pdf_file, "\n")
cat("[OK] PNG:", png_file, "\n")
