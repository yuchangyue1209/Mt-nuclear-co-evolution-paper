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