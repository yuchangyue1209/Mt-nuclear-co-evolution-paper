#fig5


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

INPUT_FILE <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_noAMO.tsv"
ANNOT_FILE <- "/work/cyu/codeml_sites_summary.merged.tsv"

OUTDIR <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/validation_results_OXPHOS72_noAMO/Figure5_three_panels"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. Read data
# ==============================================================================

dt <- fread(INPUT_FILE)
annot <- fread(ANNOT_FILE)

annot_gene <- unique(annot[, .(gene, role, complex)])
dt <- merge(dt, annot_gene, by = "gene", all.x = TRUE)

# keep OXPHOS subunits
dt_sub <- dt[role == "subunit"]

# ==============================================================================
# 3. Shared theme: only x/y axes, no full border
# ==============================================================================

theme_xy <- theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right"
  )

# ==============================================================================
# 4. Panel A: AK PBS vs BC PBS scatter
# ==============================================================================

gene_summary <- dt_sub[is.finite(nu_PBS), .(
  mean_PBS = mean(nu_PBS, na.rm = TRUE),
  n_pop = .N
), by = .(gene, region)]

ak_dt <- gene_summary[region == "AK", .(
  gene,
  mean_PBS_AK = mean_PBS,
  n_pop_AK = n_pop
)]

bc_dt <- gene_summary[region == "BC", .(
  gene,
  mean_PBS_BC = mean_PBS,
  n_pop_BC = n_pop
)]

compare_dt <- merge(ak_dt, bc_dt, by = "gene")
compare_dt <- merge(compare_dt, annot_gene, by = "gene", all.x = TRUE)
compare_dt <- compare_dt[role == "subunit"]

cutoff_ak <- quantile(compare_dt$mean_PBS_AK, 0.95, na.rm = TRUE)
cutoff_bc <- quantile(compare_dt$mean_PBS_BC, 0.95, na.rm = TRUE)

compare_dt[, status := "Background"]

compare_dt[
  mean_PBS_AK > cutoff_ak & mean_PBS_BC > cutoff_bc,
  status := "Shared"
]

compare_dt[
  mean_PBS_AK > cutoff_ak & mean_PBS_BC <= cutoff_bc,
  status := "AK-specific"
]

compare_dt[
  mean_PBS_AK <= cutoff_ak & mean_PBS_BC > cutoff_bc,
  status := "BC-specific"
]

label_A <- compare_dt[status != "Background"]

pA <- ggplot(compare_dt, aes(x = mean_PBS_BC, y = mean_PBS_AK)) +
  geom_point(aes(color = status), size = 2.6, alpha = 0.85) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = label_A,
    aes(label = gene),
    fontface = "italic",
    size = 3.6,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "AK-specific" = "#3B6FB6",
    "BC-specific" = "#2CA25F",
    "Shared" = "#8E44AD",
    "Background" = "grey80"
  )) +
  labs(
    x = "Mean nuclear PBS in British Columbia",
    y = "Mean nuclear PBS in Alaska",
    color = NULL,
    title = "A"
  ) +
  theme_xy

# ==============================================================================
# 5. Panel B: gene-level nuclear PBS vs mt PBS correlation volcano
# ==============================================================================

cor_dt <- dt_sub[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))
      .(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    } else {
      .(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    }
  },
  by = gene
]

cor_dt[, p_bh := p.adjust(p, method = "BH")]
cor_dt <- merge(cor_dt, annot_gene, by = "gene", all.x = TRUE)

cor_dt <- cor_dt[
  is.finite(cor_r) &
    is.finite(p) &
    p > 0
]

cor_dt[, logp := -log10(p)]

cor_dt[, sig := fifelse(
  p_bh < 0.10, "FDR < 0.10",
  fifelse(p < 0.05, "p < 0.05", "NS")
)]

label_B <- cor_dt[p < 0.05 | abs(cor_r) > 0.55]
label_B <- label_B[order(p)]
label_B <- label_B[1:min(.N, 12)]

pB <- ggplot(cor_dt, aes(x = cor_r, y = logp)) +
  geom_point(
    aes(color = complex, shape = sig),
    size = 2.7,
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = label_B,
    aes(label = gene),
    fontface = "italic",
    size = 3.4,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  labs(
    x = "Correlation between nuclear PBS and mtDNA PBS",
    y = expression(-log[10](p)),
    color = "Complex",
    shape = NULL,
    title = "B"
  ) +
  theme_xy

# ==============================================================================
# 6. Panel C: population-level mt PBS vs nuclear PBS regression
# ==============================================================================

plotC_dt <- dt_sub[
  is.finite(nu_PBS) &
    is.finite(mt_PBS)
]

# region-level regression stats
statC <- plotC_dt[
  ,
  {
    fit <- lm(nu_PBS ~ mt_PBS)
    sm <- summary(fit)
    .(
      r2 = sm$r.squared,
      p = coef(sm)[2, 4],
      x = min(mt_PBS, na.rm = TRUE) + 0.05 * diff(range(mt_PBS, na.rm = TRUE)),
      y = max(nu_PBS, na.rm = TRUE) - 0.08 * diff(range(nu_PBS, na.rm = TRUE))
    )
  },
  by = region
]

statC[, label := paste0(
  "R² = ", sprintf("%.2f", r2),
  "\nP = ", formatC(p, format = "e", digits = 2)
)]

pC <- ggplot(plotC_dt, aes(x = mt_PBS, y = nu_PBS)) +
  geom_point(
    aes(color = complex),
    size = 1.8,
    alpha = 0.55
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.8
  ) +
  geom_text(
    data = statC,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  facet_wrap(~ region, scales = "free") +
  labs(
    x = "mtDNA PBS",
    y = "Nuclear OXPHOS PBS",
    color = "Complex",
    title = "C"
  ) +
  theme_xy

# ==============================================================================
# 7. Save individual panels
# ==============================================================================

ggsave(
  file.path(OUTDIR, "Fig5A_AK_vs_BC_nuclear_PBS.png"),
  pA,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig5B_gene_correlation_nuclearPBS_mtPBS.png"),
  pB,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig5C_population_level_mtPBS_nuPBS_regression.png"),
  pC,
  width = 8,
  height = 5.5,
  dpi = 300
)

# ==============================================================================
# 8. Combined figure
# ==============================================================================

fig5 <- (pA + pB) / pC +
  plot_layout(heights = c(1, 0.9))

ggsave(
  file.path(OUTDIR, "Figure5_three_panels_combined.png"),
  fig5,
  width = 14,
  height = 11,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Figure5_three_panels_combined.pdf"),
  fig5,
  width = 14,
  height = 11
)

# ==============================================================================
# 9. Save tables
# ==============================================================================

fwrite(
  compare_dt,
  file.path(OUTDIR, "Fig5A_AK_BC_PBS_gene_summary.tsv"),
  sep = "\t"
)

fwrite(
  cor_dt,
  file.path(OUTDIR, "Fig5B_gene_correlation_nuclearPBS_mtPBS.tsv"),
  sep = "\t"
)

fwrite(
  statC,
  file.path(OUTDIR, "Fig5C_region_regression_stats.tsv"),
  sep = "\t"
)

# ==============================================================================
# 10. Print plots
# ==============================================================================

print(pA)
print(pB)
print(pC)
print(fig5)

cat("\nDone. Output files saved in:\n")
cat(OUTDIR, "\n")

cat("\nCheck with:\n")
cat("ls -lh ", OUTDIR, "\n", sep = "")



#

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

INPUT_FILE <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_noAMO.tsv"
ANNOT_FILE <- "/work/cyu/codeml_sites_summary.merged.tsv"

OUTDIR <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/validation_results_OXPHOS72_noAMO/Figure5_three_panels"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. Global plotting settings
# ==============================================================================

BASE_SIZE  <- 16
TITLE_SIZE <- 20
AXIS_SIZE  <- 18
TEXT_SIZE  <- 15
LABEL_SIZE <- 4.5

italic_label <- function(x) {
  paste0("italic('", x, "')")
}

theme_xy <- theme_classic(base_size = BASE_SIZE) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = TEXT_SIZE),
    plot.title = element_text(face = "bold", size = TITLE_SIZE),
    axis.title = element_text(size = AXIS_SIZE),
    axis.text = element_text(size = TEXT_SIZE),
    legend.text = element_text(size = TEXT_SIZE),
    legend.title = element_text(size = TEXT_SIZE),
    legend.position = "right"
  )

# ==============================================================================
# 3. Read data
# ==============================================================================

dt <- fread(INPUT_FILE)
annot <- fread(ANNOT_FILE)

annot_gene <- unique(annot[, .(gene, role, complex)])
dt <- merge(dt, annot_gene, by = "gene", all.x = TRUE)

dt_sub <- dt[role == "subunit"]

# ==============================================================================
# 4. Panel A: AK PBS vs BC PBS scatter
# ==============================================================================

gene_summary <- dt_sub[is.finite(nu_PBS), .(
  mean_PBS = mean(nu_PBS, na.rm = TRUE),
  n_pop = .N
), by = .(gene, region)]

ak_dt <- gene_summary[region == "AK", .(
  gene,
  mean_PBS_AK = mean_PBS,
  n_pop_AK = n_pop
)]

bc_dt <- gene_summary[region == "BC", .(
  gene,
  mean_PBS_BC = mean_PBS,
  n_pop_BC = n_pop
)]

compare_dt <- merge(ak_dt, bc_dt, by = "gene")
compare_dt <- merge(compare_dt, annot_gene, by = "gene", all.x = TRUE)
compare_dt <- compare_dt[role == "subunit"]

cutoff_ak <- quantile(compare_dt$mean_PBS_AK, 0.95, na.rm = TRUE)
cutoff_bc <- quantile(compare_dt$mean_PBS_BC, 0.95, na.rm = TRUE)

compare_dt[, status := "Background"]
compare_dt[mean_PBS_AK > cutoff_ak & mean_PBS_BC > cutoff_bc, status := "Shared"]
compare_dt[mean_PBS_AK > cutoff_ak & mean_PBS_BC <= cutoff_bc, status := "AK-specific"]
compare_dt[mean_PBS_AK <= cutoff_ak & mean_PBS_BC > cutoff_bc, status := "BC-specific"]

label_A <- compare_dt[status != "Background"]

pA <- ggplot(compare_dt, aes(x = mean_PBS_BC, y = mean_PBS_AK)) +
  geom_point(aes(color = status), size = 3.0, alpha = 0.85) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = label_A,
    aes(label = italic_label(gene)),
    parse = TRUE,
    size = LABEL_SIZE,
    color = "black",
    box.padding = 0.45,
    point.padding = 0.35,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "AK-specific" = "#3B6FB6",
    "BC-specific" = "#2CA25F",
    "Shared" = "#8E44AD",
    "Background" = "grey80"
  )) +
  labs(
    x = "Mean nuclear PBS in British Columbia",
    y = "Mean nuclear PBS in Alaska",
    color = NULL,
    title = "A"
  ) +
  theme_xy

# ==============================================================================
# 5. Panel B: gene-level nuclear PBS vs mt PBS correlation volcano
# ==============================================================================

cor_dt <- dt_sub[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))
      .(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    } else {
      .(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    }
  },
  by = gene
]

cor_dt[, p_bh := p.adjust(p, method = "BH")]
cor_dt <- merge(cor_dt, annot_gene, by = "gene", all.x = TRUE)
cor_dt <- cor_dt[is.finite(cor_r) & is.finite(p) & p > 0]
cor_dt[, logp := -log10(p)]

cor_dt[, sig := fifelse(
  p_bh < 0.10, "FDR < 0.10",
  fifelse(p < 0.05, "p < 0.05", "NS")
)]

label_B <- cor_dt[p < 0.05 | abs(cor_r) > 0.55]
label_B <- label_B[order(p)][1:min(.N, 12)]

pB <- ggplot(cor_dt, aes(x = cor_r, y = logp)) +
  geom_point(
    aes(color = complex, shape = sig),
    size = 3.0,
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "grey40",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = label_B,
    aes(label = italic_label(gene)),
    parse = TRUE,
    size = LABEL_SIZE,
    color = "black",
    box.padding = 0.45,
    point.padding = 0.35,
    max.overlaps = Inf
  ) +
  labs(
    x = "Correlation between nuclear PBS and mtDNA PBS",
    y = expression(-log[10](p)),
    color = "Complex",
    shape = NULL,
    title = "B"
  ) +
  theme_xy

# ==============================================================================
# 6. Panel C: population-level mt PBS vs nuclear PBS regression
# ==============================================================================

plotC_dt <- dt_sub[is.finite(nu_PBS) & is.finite(mt_PBS)]

pC <- ggplot(plotC_dt, aes(x = mt_PBS, y = nu_PBS)) +
  geom_point(
    aes(color = complex),
    size = 2.0,
    alpha = 0.55
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.9
  ) +
  facet_wrap(~ region, scales = "free") +
  labs(
    x = "mtDNA PBS",
    y = "Nuclear OXPHOS PBS",
    color = "Complex",
    title = "C"
  ) +
  theme_xy +
  theme(
    legend.position = "right"
  )

# ==============================================================================
# 7. Combined figure
# ==============================================================================

fig5 <- (pA + pB) / pC +
  plot_layout(heights = c(1, 0.9))

# ==============================================================================
# 8. Save outputs
# ==============================================================================

ggsave(
  file.path(OUTDIR, "Figure5_three_panels_combined.png"),
  fig5,
  width = 15,
  height = 12,
  dpi = 600
)

ggsave(
  file.path(OUTDIR, "Figure5_three_panels_combined.pdf"),
  fig5,
  width = 15,
  height = 12,
  device = cairo_pdf
)

ggsave(file.path(OUTDIR, "Fig5A_AK_vs_BC_nuclear_PBS.pdf"), pA, width = 7.5, height = 6, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig5B_gene_correlation_nuclearPBS_mtPBS.pdf"), pB, width = 7.5, height = 6, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig5C_population_level_mtPBS_nuPBS_regression.pdf"), pC, width = 12, height = 5.8, device = cairo_pdf)

fwrite(compare_dt, file.path(OUTDIR, "Fig5A_AK_BC_PBS_gene_summary.tsv"), sep = "\t")
fwrite(cor_dt, file.path(OUTDIR, "Fig5B_gene_correlation_nuclearPBS_mtPBS.tsv"), sep = "\t")



cat("\nDone. Output files saved in:\n")
cat(OUTDIR, "\n")











#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ============================================================
# Input paths
# ============================================================

PBS_FILE <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/merged_mt_nuclear_PBS_NorwayOutgroup_by_gene_population.tsv"
ANNOT_FILE <- "/work/cyu/codeml_sites_summary.merged.tsv"

OUTDIR <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/non_subunit_gene_set_figures"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Global font settings
# ============================================================

BASE_SIZE  <- 16
TITLE_SIZE <- 18
AXIS_SIZE  <- 17
TEXT_SIZE  <- 14
LABEL_SIZE <- 5

# ============================================================
# Read data
# ============================================================

dt <- fread(PBS_FILE)
annot <- fread(ANNOT_FILE)

annot_gene <- unique(annot[, .(gene, role, complex)])

dt <- merge(dt, annot_gene, by = "gene", all.x = TRUE)
dt <- dt[focal != "AMO"]

role_levels <- c(
  "assembly_factor",
  "Nmt-ARS",
  "Nmt-ribo",
  "cyto-ARS",
  "cyto-ribo"
)

role_labels <- c(
  "assembly_factor" = "Assembly factors",
  "Nmt-ARS" = "Mitochondrial ARS",
  "Nmt-ribo" = "Mitochondrial ribosomal proteins",
  "cyto-ARS" = "Cytosolic ARS",
  "cyto-ribo" = "Cytosolic ribosomal proteins"
)

dt <- dt[role %in% role_levels]

# ============================================================
# Helper functions
# ============================================================

make_region_table <- function(x) {
  gene_summary <- x[
    is.finite(nu_PBS),
    .(
      mean_PBS = mean(nu_PBS, na.rm = TRUE),
      n_pop = .N
    ),
    by = .(gene, region, role)
  ]

  ak <- gene_summary[region == "AK", .(
    gene,
    role,
    mean_PBS_AK = mean_PBS
  )]

  bc <- gene_summary[region == "BC", .(
    gene,
    role,
    mean_PBS_BC = mean_PBS
  )]

  out <- merge(ak, bc, by = c("gene", "role"))

  cutoff_ak <- quantile(out$mean_PBS_AK, 0.95, na.rm = TRUE)
  cutoff_bc <- quantile(out$mean_PBS_BC, 0.95, na.rm = TRUE)

  out[, status := "Background"]

  out[
    mean_PBS_AK > cutoff_ak & mean_PBS_BC > cutoff_bc,
    status := "Shared"
  ]

  out[
    mean_PBS_AK > cutoff_ak & mean_PBS_BC <= cutoff_bc,
    status := "AK-specific"
  ]

  out[
    mean_PBS_AK <= cutoff_ak & mean_PBS_BC > cutoff_bc,
    status := "BC-specific"
  ]

  out
}

make_corr_table <- function(x) {
  res <- x[
    is.finite(nu_PBS) & is.finite(mt_PBS),
    {
      if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
        ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))
        .(
          n_pop = .N,
          cor_r = unname(ct$estimate),
          p = ct$p.value
        )
      } else {
        .(
          n_pop = .N,
          cor_r = NA_real_,
          p = NA_real_
        )
      }
    },
    by = .(gene, role)
  ]

  res[, p_bh := p.adjust(p, method = "BH")]
  res[, sig := ifelse(is.finite(p) & p < 0.05, "p < 0.05", "NS")]

  res
}

italic_label <- function(x) {
  paste0("italic('", x, "')")
}

common_theme <- function() {
  theme_classic(base_size = BASE_SIZE) +
    theme(
      axis.title = element_text(size = AXIS_SIZE),
      axis.text = element_text(size = TEXT_SIZE),
      legend.text = element_text(size = TEXT_SIZE),
      legend.title = element_text(size = TEXT_SIZE),
      plot.title = element_text(size = TITLE_SIZE, face = "bold"),
      legend.position = "bottom"
    )
}

# ============================================================
# Plot one role
# ============================================================

plot_one_role <- function(role_i, panel_letter) {
  x <- dt[role == role_i]
  label_i <- role_labels[[role_i]]

  region_dt <- make_region_table(x)
  corr_dt <- make_corr_table(x)

  fwrite(region_dt, file.path(OUTDIR, paste0(role_i, "_AK_vs_BC.tsv")), sep = "\t")
  fwrite(corr_dt, file.path(OUTDIR, paste0(role_i, "_correlation.tsv")), sep = "\t")

  label_dt <- region_dt[status != "Background"]

  pA <- ggplot(
    region_dt,
    aes(x = mean_PBS_BC, y = mean_PBS_AK, color = status)
  ) +
    geom_point(size = 3.2, alpha = 0.8) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_text_repel(
      data = label_dt,
      aes(label = italic_label(gene)),
      parse = TRUE,
      size = LABEL_SIZE,
      color = "black",
      box.padding = 0.45,
      point.padding = 0.35,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = c(
      "AK-specific" = "#377EB8",
      "BC-specific" = "#4DAF4A",
      "Shared" = "#984EA3",
      "Background" = "grey80"
    )) +
    common_theme() +
    labs(
      title = paste0(panel_letter, ". ", label_i),
      x = "Mean nuclear PBS in British Columbia",
      y = "Mean nuclear PBS in Alaska",
      color = NULL
    )

  corr_label <- corr_dt[
    is.finite(cor_r) & is.finite(p)
  ][order(p)][1:min(5, .N)]

  pB <- ggplot(
    corr_dt[is.finite(cor_r) & is.finite(p)],
    aes(x = cor_r, y = -log10(p))
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "grey60"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_point(
      aes(shape = sig),
      size = 3.2,
      alpha = 0.9
    ) +
    geom_text_repel(
      data = corr_label,
      aes(label = italic_label(gene)),
      parse = TRUE,
      size = LABEL_SIZE,
      color = "black",
      box.padding = 0.45,
      point.padding = 0.35,
      max.overlaps = Inf
    ) +
    scale_shape_manual(values = c(
      "NS" = 16,
      "p < 0.05" = 17
    )) +
    common_theme() +
    labs(
      x = "Correlation between nuclear PBS and mtDNA PBS",
      y = expression(-log[10](p)),
      shape = NULL
    )

  combined <- pA + pB + plot_layout(widths = c(1, 1))

  ggsave(
    file.path(OUTDIR, paste0("Fig_", role_i, ".png")),
    combined,
    width = 16,
    height = 7,
    dpi = 600
  )

  ggsave(
    file.path(OUTDIR, paste0("Fig_", role_i, ".pdf")),
    combined,
    width = 16,
    height = 7,
    device = cairo_pdf
  )

  combined
}

# ============================================================
# Generate panels A-E
# ============================================================

plots <- list()
panel_letters <- c("A", "B", "C", "D", "E")

for (i in seq_along(role_levels)) {
  plots[[i]] <- plot_one_role(role_levels[i], panel_letters[i])
}

final_plot <- wrap_plots(plots, ncol = 1)

# ============================================================
# Save final figure
# ============================================================

ggsave(
  file.path(OUTDIR, "Fig_non_subunit_gene_sets_PBS_A_to_E.png"),
  final_plot,
  width = 16,
  height = 32,
  dpi = 600
)

ggsave(
  file.path(OUTDIR, "Fig_non_subunit_gene_sets_PBS_A_to_E.pdf"),
  final_plot,
  width = 16,
  height = 32,
  device = cairo_pdf
)

# For RStudio Server: print individual panels if full plot is too large
print(plots[[1]])
print(plots[[2]])
print(plots[[3]])
print(plots[[4]])
print(plots[[5]])

cat("\nDONE\n")
cat("Output directory:\n")
cat(OUTDIR, "\n")











#drop lb





#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ============================================================
# Figure S5: non-subunit mitochondrial-associated gene sets
# Norway outgroup; LB and AMO excluded
# ============================================================

# -------------------------
# Input paths
# -------------------------
PBS_FILE <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup_dropLB/merged_mt_nuclear_PBS_NorwayOutgroup_dropLB_by_gene_population.tsv"
ANNOT_FILE <- "/work/cyu/codeml_sites_summary.merged.tsv"

OUTDIR <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup_dropLB/FigS5_non_subunit_gene_sets_dropLB_noAMO"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Plot settings
# -------------------------
BASE_SIZE  <- 16
TITLE_SIZE <- 18
AXIS_SIZE  <- 17
TEXT_SIZE  <- 14
LABEL_SIZE <- 4.4

# Set TRUE if you want ggrepel labels.
# If RStudio Server gives viewport errors, keep FALSE.
ADD_REPEL_LABELS <- FALSE

# Use simple geom_text labels instead of ggrepel for stability
ADD_SIMPLE_LABELS <- TRUE

theme_s5 <- theme_classic(base_size = BASE_SIZE) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = AXIS_SIZE),
    axis.text = element_text(size = TEXT_SIZE),
    legend.text = element_text(size = TEXT_SIZE),
    legend.title = element_text(size = TEXT_SIZE),
    plot.title = element_text(size = TITLE_SIZE, face = "bold"),
    legend.position = "bottom"
  )

# ============================================================
# 1. Read data
# ============================================================

dt <- fread(PBS_FILE)
annot <- fread(ANNOT_FILE)

annot_gene <- unique(annot[, .(gene, role, complex)])
dt <- merge(dt, annot_gene, by = "gene", all.x = TRUE)

# Drop AMO and LB
dt <- dt[!focal %in% c("AMO", "LB")]

role_levels <- c(
  "assembly_factor",
  "Nmt-ARS",
  "Nmt-ribo",
  "cyto-ARS",
  "cyto-ribo"
)

role_labels <- c(
  "assembly_factor" = "Assembly factors",
  "Nmt-ARS" = "Mitochondrial ARS",
  "Nmt-ribo" = "Mitochondrial ribosomal proteins",
  "cyto-ARS" = "Cytosolic ARS",
  "cyto-ribo" = "Cytosolic ribosomal proteins"
)

dt <- dt[role %in% role_levels]

cat("\n[Check] populations used:\n")
print(dt[, .(
  n_pop = uniqueN(focal),
  populations = paste(sort(unique(focal)), collapse = ", ")
), by = region])

cat("\n[Check] genes by role:\n")
print(dt[, .(n_gene = uniqueN(gene)), by = role][order(role)])

# ============================================================
# 2. Helper functions
# ============================================================

make_region_table <- function(x) {
  gene_summary <- x[
    is.finite(nu_PBS),
    .(
      mean_PBS = mean(nu_PBS, na.rm = TRUE),
      median_PBS = median(nu_PBS, na.rm = TRUE),
      max_PBS = max(nu_PBS, na.rm = TRUE),
      n_pop = .N
    ),
    by = .(gene, region, role)
  ]

  ak <- gene_summary[region == "AK", .(
    gene,
    role,
    mean_PBS_AK = mean_PBS,
    median_PBS_AK = median_PBS,
    max_PBS_AK = max_PBS,
    n_pop_AK = n_pop
  )]

  bc <- gene_summary[region == "BC", .(
    gene,
    role,
    mean_PBS_BC = mean_PBS,
    median_PBS_BC = median_PBS,
    max_PBS_BC = max_PBS,
    n_pop_BC = n_pop
  )]

  out <- merge(ak, bc, by = c("gene", "role"))

  if (nrow(out) == 0) {
    return(out)
  }

  cutoff_ak <- quantile(out$mean_PBS_AK, 0.95, na.rm = TRUE)
  cutoff_bc <- quantile(out$mean_PBS_BC, 0.95, na.rm = TRUE)

  out[, cutoff_AK_95 := cutoff_ak]
  out[, cutoff_BC_95 := cutoff_bc]

  out[, status := "Background"]

  out[
    mean_PBS_AK > cutoff_ak & mean_PBS_BC > cutoff_bc,
    status := "Shared"
  ]

  out[
    mean_PBS_AK > cutoff_ak & mean_PBS_BC <= cutoff_bc,
    status := "AK-specific"
  ]

  out[
    mean_PBS_AK <= cutoff_ak & mean_PBS_BC > cutoff_bc,
    status := "BC-specific"
  ]

  out[, PBS_diff_AK_minus_BC := mean_PBS_AK - mean_PBS_BC]
  out[, abs_PBS_diff := abs(PBS_diff_AK_minus_BC)]

  out[]
}

make_corr_table <- function(x) {
  res <- x[
    is.finite(nu_PBS) & is.finite(mt_PBS),
    {
      if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
        ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))

        .(
          n_pop = .N,
          cor_r = unname(ct$estimate),
          p = ct$p.value,
          mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
          mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
        )
      } else {
        .(
          n_pop = .N,
          cor_r = NA_real_,
          p = NA_real_,
          mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
          mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
        )
      }
    },
    by = .(gene, role)
  ]

  res[, p_bh := p.adjust(p, method = "BH")]
  res[, sig := fifelse(
    is.finite(p_bh) & p_bh < 0.10, "FDR < 0.10",
    fifelse(is.finite(p) & p < 0.05, "p < 0.05", "NS")
  )]

  res[]
}

# ============================================================
# 3. Plot one role
# ============================================================

plot_one_role <- function(role_i, panel_letter) {
  x <- dt[role == role_i]
  label_i <- role_labels[[role_i]]

  cat("\n[Running] ", role_i, "\n", sep = "")

  region_dt <- make_region_table(x)
  corr_dt <- make_corr_table(x)

  # Save tables
  fwrite(
    region_dt,
    file.path(OUTDIR, paste0(role_i, "_AK_vs_BC_meanPBS_dropLB_noAMO.tsv")),
    sep = "\t"
  )

  fwrite(
    corr_dt,
    file.path(OUTDIR, paste0(role_i, "_nuPBS_mtPBS_correlation_dropLB_noAMO.tsv")),
    sep = "\t"
  )

  if (nrow(region_dt) == 0 || nrow(corr_dt) == 0) {
    warning("No data for role: ", role_i)
    return(NULL)
  }

  label_region <- region_dt[status != "Background"]

  # If no status outliers due to small category, label top AK/BC diff genes
  if (nrow(label_region) == 0) {
    label_region <- region_dt[order(-abs_PBS_diff)][1:min(.N, 4)]
  }

  # -------------------------
  # Left panel: AK vs BC mean PBS
  # -------------------------
  p_left <- ggplot(
    region_dt,
    aes(x = mean_PBS_BC, y = mean_PBS_AK)
  ) +
    geom_point(aes(color = status), size = 3.2, alpha = 0.85) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.6
    ) +
    scale_color_manual(values = c(
      "AK-specific" = "#377EB8",
      "BC-specific" = "#4DAF4A",
      "Shared" = "#984EA3",
      "Background" = "grey80"
    )) +
    labs(
      title = paste0(panel_letter, ". ", label_i),
      x = "Mean nuclear PBS in British Columbia",
      y = "Mean nuclear PBS in Alaska",
      color = NULL
    ) +
    theme_s5

  if (ADD_REPEL_LABELS && nrow(label_region) > 0) {
    p_left <- p_left +
      geom_text_repel(
        data = label_region,
        aes(label = gene),
        fontface = "italic",
        size = LABEL_SIZE,
        color = "black",
        box.padding = 0.45,
        point.padding = 0.35,
        max.overlaps = Inf
      )
  }

  if (!ADD_REPEL_LABELS && ADD_SIMPLE_LABELS && nrow(label_region) > 0) {
    p_left <- p_left +
      geom_text(
        data = label_region,
        aes(label = gene),
        fontface = "italic",
        size = 4,
        color = "black",
        vjust = -0.7,
        check_overlap = TRUE
      )
  }

  # -------------------------
  # Right panel: correlation volcano
  # -------------------------
  corr_plot_dt <- corr_dt[is.finite(cor_r) & is.finite(p) & p > 0]
  corr_plot_dt[, logp := -log10(p)]

  corr_label <- corr_plot_dt[order(p)][1:min(.N, 5)]

  p_right <- ggplot(
    corr_plot_dt,
    aes(x = cor_r, y = logp)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "grey60",
      linewidth = 0.6
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.6
    ) +
    geom_point(
      aes(shape = sig),
      size = 3.2,
      alpha = 0.9
    ) +
    scale_shape_manual(values = c(
      "NS" = 16,
      "p < 0.05" = 17,
      "FDR < 0.10" = 18
    )) +
    labs(
      x = "Correlation between nuclear PBS and mtDNA PBS",
      y = expression(-log[10](p)),
      shape = NULL
    ) +
    theme_s5

  if (ADD_REPEL_LABELS && nrow(corr_label) > 0) {
    p_right <- p_right +
      geom_text_repel(
        data = corr_label,
        aes(label = gene),
        fontface = "italic",
        size = LABEL_SIZE,
        color = "black",
        box.padding = 0.45,
        point.padding = 0.35,
        max.overlaps = Inf
      )
  }

  if (!ADD_REPEL_LABELS && ADD_SIMPLE_LABELS && nrow(corr_label) > 0) {
    p_right <- p_right +
      geom_text(
        data = corr_label,
        aes(label = gene),
        fontface = "italic",
        size = 4,
        color = "black",
        vjust = -0.7,
        check_overlap = TRUE
      )
  }

  combined <- p_left + p_right + plot_layout(widths = c(1, 1))

  # Save each role
  ggsave(
    file.path(OUTDIR, paste0("FigS5_", role_i, "_dropLB_noAMO.png")),
    combined,
    width = 16,
    height = 7,
    dpi = 600
  )

  ggsave(
    file.path(OUTDIR, paste0("FigS5_", role_i, "_dropLB_noAMO.pdf")),
    combined,
    width = 16,
    height = 7,
    device = cairo_pdf
  )

  # PDF print version for RStudio Server
  pdf(
    file.path(OUTDIR, paste0("PRINT_FigS5_", role_i, "_dropLB_noAMO.pdf")),
    width = 16,
    height = 7
  )
  print(combined)
  dev.off()

  combined
}

# ============================================================
# 4. Generate panels A-E
# ============================================================

plots <- list()
panel_letters <- c("A", "B", "C", "D", "E")

for (i in seq_along(role_levels)) {
  plots[[i]] <- plot_one_role(role_levels[i], panel_letters[i])
}

plots <- Filter(Negate(is.null), plots)

final_plot <- wrap_plots(plots, ncol = 1)

# ============================================================
# 5. Save final Figure S5
# ============================================================

ggsave(
  file.path(OUTDIR, "FigS5_non_subunit_gene_sets_PBS_A_to_E_dropLB_noAMO.png"),
  final_plot,
  width = 16,
  height = 32,
  dpi = 600
)

ggsave(
  file.path(OUTDIR, "FigS5_non_subunit_gene_sets_PBS_A_to_E_dropLB_noAMO.pdf"),
  final_plot,
  width = 16,
  height = 32,
  device = cairo_pdf
)

# Safe print-to-PDF version
pdf(
  file.path(OUTDIR, "PRINT_FigS5_non_subunit_gene_sets_PBS_A_to_E_dropLB_noAMO.pdf"),
  width = 16,
  height = 32
)
print(final_plot)
dev.off()

# ============================================================
# 6. Summary tables
# ============================================================

summary_by_role <- dt[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  .(
    n_gene = uniqueN(gene),
    n_rows = .N,
    n_pop = uniqueN(focal),
    n_positive_cor = {
      tmp <- make_corr_table(.SD)
      tmp[is.finite(cor_r) & cor_r > 0, .N]
    },
    n_nominal_p05 = {
      tmp <- make_corr_table(.SD)
      tmp[is.finite(p) & p < 0.05, .N]
    }
  ),
  by = role
]

fwrite(
  summary_by_role,
  file.path(OUTDIR, "FigS5_summary_by_role_dropLB_noAMO.tsv"),
  sep = "\t"
)

cat("\nDONE\n")
cat("Output directory:\n")
cat(OUTDIR, "\n")
cat("\nCheck with:\n")
cat("ls -lh ", OUTDIR, "\n", sep = "")