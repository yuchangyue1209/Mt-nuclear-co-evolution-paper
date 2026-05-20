#figs3
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