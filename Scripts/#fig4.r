#fig4
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
})

# =========================
# Paths
# =========================

IN_LM <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm/LM_perSNP_mtCluster_plus_treePC12.LDpruned_r2_0.2.tsv.gz"

IN_ENRICH <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm/gene_enrichment/Gene_enrichment_overall.tsv.gz"

OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm/Figure4_unified_singlePanels"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# =========================
# Unified style
# =========================

COL_CLUSTER <- c(
  "C1" = "#F8766D",
  "C2" = "#7CAE00",
  "C3" = "#00BFC4",
  "C4" = "#C77CFF"
)

COL_WITHIN <- c(
  "Within" = "#F8766D",
  "Between" = "#00BFC4"
)

theme_fig <- theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm")
  )

# =========================
# Panel A: mtCluster volcano
# =========================

LM <- fread(cmd = paste("zcat", shQuote(IN_LM)))

LM <- LM[
  is.finite(p_cluster) &
    p_cluster > 0 &
    p_cluster <= 1
]

LM[, logp := -log10(p_cluster)]
LM[, driver_mu := as.numeric(driver_mu)]

VOL <- LM[
  is.finite(driver_mu) &
    is.finite(logp)
]

VOL[, sig := fifelse(p_cluster < 0.05, "p < 0.05", "NS")]

TOP_VOL <- VOL[logp > 5]
TOP_VOL <- TOP_VOL[order(region, -logp)]
TOP_VOL <- TOP_VOL[, .SD[!duplicated(gene)], by = region]
TOP_VOL <- TOP_VOL[, head(.SD, 7), by = region]

# cap extreme values only for plotting
VOL[, logp_plot := pmin(logp, 60)]
TOP_VOL[, logp_plot := pmin(logp, 60)]

pA <- ggplot(VOL, aes(x = driver_mu, y = logp_plot)) +
  geom_point(
    aes(color = driver_cluster, shape = sig),
    alpha = 0.75,
    size = 1.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = TOP_VOL,
    aes(y = logp_plot, label = gene),
    size = 3.5,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.3,
    max.overlaps = Inf,
    force = 2
  ) +
  facet_wrap(~ region, nrow = 1) +
  coord_cartesian(ylim = c(0, 65)) +
  scale_color_manual(values = COL_CLUSTER, name = "Driver cluster") +
  scale_shape_manual(
    values = c("NS" = 16, "p < 0.05" = 17),
    name = NULL
  ) +
  labs(
    title = "A",
    x = expression("Driver-cluster mean " * Delta * "AF"),
    y = expression(-log[10](p))
  ) +
  theme_fig +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
  )

ggsave(
  file.path(OUTDIR, "Fig4A_mtCluster_volcano_unified.png"),
  pA,
  width = 8.2,
  height = 5.8,
  dpi = 300
)

print(pA)

# =========================
# Panel B: gene enrichment heatmap
# =========================

G <- fread(IN_ENRICH)

G[, log2OR := log2(or_overall)]
G[!is.finite(log2OR), log2OR := NA_real_]

G[, sig_label := fifelse(
  q_overall < 0.10,
  "+",
  fifelse(p_overall < 0.05, "*", "")
)]

top_genes <- unique(
  G[order(p_overall), head(gene, 20), by = region]$V1
)

P <- G[gene %in% top_genes]

gene_order <- P[
  ,
  .(best_p = min(p_overall, na.rm = TRUE)),
  by = gene
][order(best_p)]$gene

P[, gene := factor(gene, levels = rev(gene_order))]
P[, region := factor(region, levels = c("AK", "BC"))]

pB <- ggplot(P, aes(x = region, y = gene, fill = log2OR)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig_label), size = 3.4, color = "black") +
  scale_fill_gradient2(
    low = "#4C78A8",
    mid = "white",
    high = "#E45756",
    midpoint = 0,
    name = "log2(OR)"
  ) +
  labs(
    title = "B",
    x = NULL,
    y = NULL
  ) +
  theme_fig +
  theme(
    axis.text.y = element_text(
      size = 8.5,
      face = "italic"
    ),
    legend.position = "right"
  )

ggsave(
  file.path(OUTDIR, "Fig4B_gene_enrichment_heatmap_unified.png"),
  pB,
  width = 4.8,
  height = 6.4,
  dpi = 300
)

print(pB)

# =========================
# Panel C: profile distance barplot
# =========================

BAR <- data.table(
  region = c("AK", "AK", "BC", "BC"),
  group = c("Within", "Between", "Within", "Between"),
  mean = c(0.395, 0.355, 0.510, 0.625),
  se = c(0.055, 0.030, 0.070, 0.025),
  p_label = c("p = 0.82", "p = 0.82", "p = 0.029 *", "p = 0.029 *")
)

BAR[, group := factor(group, levels = c("Within", "Between"))]
BAR[, region := factor(region, levels = c("AK", "BC"))]

LAB <- BAR[
  ,
  .(
    p_label = unique(p_label)[1],
    x = 1.5,
    y = max(mean + se) + 0.045
  ),
  by = region
]

pC <- ggplot(BAR, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.62, color = NA, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.18,
    linewidth = 0.5
  ) +
  geom_text(
    data = LAB,
    aes(x = x, y = y, label = p_label),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~ region, nrow = 1) +
  scale_fill_manual(values = COL_WITHIN) +
  coord_cartesian(ylim = c(0, 0.72)) +
  labs(
    title = "C",
    x = NULL,
    y = expression("1 - correlation of " * Delta * "AF profiles")
  ) +
  theme_fig +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

ggsave(
  file.path(OUTDIR, "Fig4C_deltaAF_profile_bar_unified.png"),
  pC,
  width = 6.0,
  height = 4.8,
  dpi = 300
)

print(pC)

# =========================
# Panel D: Venn-style overlap
# =========================

VENN <- data.table(
  label = c(
    "610\nmtCluster",
    "830",
    "18217\nGeography",
    "23079\nOXPHOS SNPs"
  ),
  x = c(0.30, 0.50, 0.68, 0.50),
  y = c(0.50, 0.50, 0.50, 0.72)
)

pD <- ggplot() +
  geom_ellipse(
    aes(x0 = 0.50, y0 = 0.50, a = 0.42, b = 0.30, angle = 0),
    fill = "grey92",
    color = "grey45",
    linewidth = 0.7
  ) +
  
  geom_ellipse(
    aes(x0 = 0.36, y0 = 0.50, a = 0.20, b = 0.16, angle = 0),
    fill = "#F8766D",
    alpha = 0.28,
    color = "grey35",
    linewidth = 0.7
  ) +
  
  geom_ellipse(
    aes(x0 = 0.64, y0 = 0.50, a = 0.23, b = 0.17, angle = 0),
    fill = "#00BFC4",
    alpha = 0.28,
    color = "grey35",
    linewidth = 0.7
  ) +
  geom_text(
    data = VENN,
    aes(x = x, y = y, label = label),
    size = 4.6,
    fontface = "bold",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 0.04,
    y = 0.94,
    label = "D",
    fontface = "bold",
    size = 7,
    hjust = 0
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0.15, 0.96)) +
  theme_void()

ggsave(
  file.path(OUTDIR, "Fig4D_overlap_venn_unified.png"),
  pD,
  width = 5.8,
  height = 4.6,
  dpi = 300
)

print(pD)

cat("\nDone. Files saved in:\n")
cat(OUTDIR, "\n")
cat("\nCheck with:\n")
cat("ls -lh ", OUTDIR, "\n", sep = "")










suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})

# ============================================================
# Paths
# ============================================================

GLM_DIR <- "/work/cyu/ldx_all_subunits/ld/mtlineage_quasibinomial_GLM_OXPHOS72/separate_geo_mtlineage_GLM"

IN_GEO <- file.path(
  GLM_DIR,
  "OXPHOS72_geography_quasibinomial_GLM_results.tsv"
)

IN_MT <- file.path(
  GLM_DIR,
  "OXPHOS72_mtlineage_quasibinomial_GLM_results.tsv"
)

PAIRWISE_FILE <- "/work/cyu/ldx_all_subunits/ld/pairwise_deltaAF_similarity_ldpruned_effectFiltered_10000perm/pairwise_deltaAF_similarity_within_between_effect0.2_10000perm.tsv"

SUMMARY_FILE <- "/work/cyu/ldx_all_subunits/ld/pairwise_deltaAF_similarity_ldpruned_effectFiltered_10000perm/pairwise_deltaAF_similarity_summary_effect0.2_10000perm.tsv"

OUTDIR <- "/work/cyu/ldx_all_subunits/ld/mtlineage_quasibinomial_GLM_OXPHOS72/Figure4_threePanels_originalStyle_realSE_print"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parameters
# ============================================================

Q_THRESHOLD <- 0.05
EFFECT_THRESHOLD <- 0.20

# ============================================================
# Unified style
# ============================================================

COL_CAND <- c(
  "Not candidate" = "grey75",
  "FDR only" = "#4C78A8",
  "Large effect only" = "#E45756",
  "Candidate" = "#C77CFF"
)

COL_WITHIN <- c(
  "Within" = "#F8766D",
  "Between" = "#00BFC4"
)

theme_fig <- theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    plot.subtitle = element_text(size = 13, hjust = 0),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm")
  )

# ============================================================
# Panel A: Geography GLM volcano
# ============================================================

GEO <- fread(IN_GEO)

GEO <- GEO[
  is.finite(q_geo) &
    q_geo > 0 &
    q_geo <= 1 &
    is.finite(abs_geo_AF_diff)
]

GEO[, logq := -log10(q_geo)]

GEO[, plot_group := "Not candidate"]
GEO[significant_geo == TRUE & large_effect_geo == FALSE, plot_group := "FDR only"]
GEO[significant_geo == FALSE & large_effect_geo == TRUE, plot_group := "Large effect only"]
GEO[candidate_geo == TRUE, plot_group := "Candidate"]

GEO[, plot_group := factor(
  plot_group,
  levels = c("Not candidate", "FDR only", "Large effect only", "Candidate")
)]

# label each candidate gene once
TOP_GEO <- GEO[candidate_geo == TRUE][order(q_geo, -abs_geo_AF_diff)]
TOP_GEO <- TOP_GEO[, .SD[1], by = gene]

GEO[, logq_plot := pmin(logq, 2)]
TOP_GEO[, logq_plot := pmin(logq, 2)]

pA <- ggplot(GEO, aes(x = abs_geo_AF_diff, y = logq_plot)) +
  geom_point(
    aes(color = plot_group),
    alpha = 0.75,
    size = 1.6
  ) +
  geom_hline(
    yintercept = -log10(Q_THRESHOLD),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = EFFECT_THRESHOLD,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = TOP_GEO,
    aes(label = gene),
    size = 3.8,
    color = "black",
    fontface = "italic",
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.3,
    max.overlaps = Inf,
    force = 2
  ) +
  coord_cartesian(ylim = c(0, 1.6)) +
  scale_color_manual(values = COL_CAND, name = NULL) +
  labs(
    title = "A",
    subtitle = "Geography-associated SNPs",
    x = "Absolute AK-BC allele-frequency difference",
    y = expression(-log[10]("FDR q-value"))
  ) +
  theme_fig +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
  )

# ============================================================
# Panel B: mt-lineage GLM volcano
# ============================================================

MT <- fread(IN_MT)

MT <- MT[
  is.finite(q_mt) &
    q_mt > 0 &
    q_mt <= 1 &
    is.finite(max_lineage_AF_diff)
]

MT[, logq := -log10(q_mt)]

MT[, plot_group := "Not candidate"]
MT[significant_mt == TRUE & large_effect_mt == FALSE, plot_group := "FDR only"]
MT[significant_mt == FALSE & large_effect_mt == TRUE, plot_group := "Large effect only"]
MT[candidate_mt == TRUE, plot_group := "Candidate"]

MT[, plot_group := factor(
  plot_group,
  levels = c("Not candidate", "FDR only", "Large effect only", "Candidate")
)]

TOP_MT <- MT[candidate_mt == TRUE][order(q_mt, -max_lineage_AF_diff)]
TOP_MT <- TOP_MT[, .SD[1], by = gene]

# cap extreme small-effect SNPs only for plotting
MT[, logq_plot := pmin(logq, 6)]
TOP_MT[, logq_plot := pmin(logq, 6)]

pB <- ggplot(MT, aes(x = max_lineage_AF_diff, y = logq_plot)) +
  geom_point(
    aes(color = plot_group),
    alpha = 0.75,
    size = 1.6
  ) +
  geom_hline(
    yintercept = -log10(Q_THRESHOLD),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = EFFECT_THRESHOLD,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = TOP_MT,
    aes(label = gene),
    size = 3.8,
    color = "black",
    fontface = "italic",
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.3,
    max.overlaps = Inf,
    force = 2
  ) +
  coord_cartesian(ylim = c(0, 6)) +
  scale_color_manual(values = COL_CAND, name = NULL) +
  labs(
    title = "B",
    subtitle = "Mitochondrial-lineage-associated SNPs",
    x = "Maximum allele-frequency difference among mt lineages",
    y = expression(-log[10]("FDR q-value"))
  ) +
  theme_fig +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
  )

# ============================================================
# Panel C: signed deltaAF profile distance with real SE
# ============================================================

PAIRWISE <- fread(PAIRWISE_FILE)
SUMMARY <- fread(SUMMARY_FILE)

P <- PAIRWISE[metric == "Signed_deltaAF"]

P[, type := as.character(type)]
P[type %in% c("within", "WITHIN"), type := "Within"]
P[type %in% c("between", "BETWEEN"), type := "Between"]

P[, type := factor(type, levels = c("Within", "Between"))]
P[, region := factor(region, levels = c("AK", "BC"))]

BAR <- P[, .(
  mean = mean(dist, na.rm = TRUE),
  se = sd(dist, na.rm = TRUE) / sqrt(.N),
  sd = sd(dist, na.rm = TRUE),
  n_pairs = .N
), by = .(region, type)]

setorder(BAR, region, type)

cat("\n[INFO] Panel C BAR table with real SE:\n")
print(BAR)

# ============================================================
# Manual p labels for Panel C
# ============================================================

LAB <- data.table(
  region = factor(c("AK", "BC"), levels = c("AK", "BC")),
  x = 1.5,
  p_label = c("p = 0.107", "p = 0.0491 *")
)

YLAB <- BAR[, .(
  y = max(mean + se, na.rm = TRUE) + 0.045
), by = region]

YLAB[, region := factor(region, levels = c("AK", "BC"))]

LAB <- merge(LAB, YLAB, by = "region", all.x = TRUE)

pC <- ggplot(BAR, aes(x = type, y = mean, fill = type)) +
  geom_col(width = 0.62, color = NA, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.18,
    linewidth = 0.5
  ) +
  geom_text(
    data = LAB,
    aes(x = x, y = y, label = p_label),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~ region, nrow = 1) +
  scale_fill_manual(values = COL_WITHIN) +
  coord_cartesian(
    ylim = c(0, max(BAR$mean + BAR$se, na.rm = TRUE) + 0.12)
  ) +
  labs(
    title = "C",
    x = NULL,
    y = expression("1 - correlation of signed " * Delta * "AF")
  ) +
  theme_fig +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
  )

# ============================================================
# Save single panels
# ============================================================

ggsave(
  file.path(OUTDIR, "Fig4A_geography_GLM_volcano_unified.png"),
  pA,
  width = 7.2,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig4A_geography_GLM_volcano_unified.pdf"),
  pA,
  width = 7.2,
  height = 5.5
)

ggsave(
  file.path(OUTDIR, "Fig4B_mtlineage_GLM_volcano_unified.png"),
  pB,
  width = 7.2,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig4B_mtlineage_GLM_volcano_unified.pdf"),
  pB,
  width = 7.2,
  height = 5.5
)

ggsave(
  file.path(OUTDIR, "Fig4C_signed_deltaAF_profile_bar_realSE_unified.png"),
  pC,
  width = 6.0,
  height = 4.8,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig4C_signed_deltaAF_profile_bar_realSE_unified.pdf"),
  pC,
  width = 6.0,
  height = 4.8
)

fwrite(
  BAR,
  file.path(OUTDIR, "Fig4C_signed_deltaAF_profile_bar_realSE_values.tsv"),
  sep = "\t"
)

# ============================================================
# Combined three-panel figure
# ============================================================

OUT_PANEL_PNG <- file.path(
  OUTDIR,
  "Fig4_threePanels_GLM_profile_realSE_unified.png"
)

OUT_PANEL_PDF <- file.path(
  OUTDIR,
  "Fig4_threePanels_GLM_profile_realSE_unified.pdf"
)

draw_three_panel <- function(use_label = TRUE) {
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        widths = unit(c(1, 1), "null"),
        heights = unit(c(1, 0.85), "null")
      )
    )
  )
  
  print(
    pA + theme(legend.position = "none"),
    vp = viewport(layout.pos.row = 1, layout.pos.col = 1)
  )
  
  print(
    pB + theme(legend.position = "none"),
    vp = viewport(layout.pos.row = 1, layout.pos.col = 2)
  )
  
  print(
    pC,
    vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2)
  )
}

png(
  filename = OUT_PANEL_PNG,
  width = 3000,
  height = 2400,
  res = 300
)

draw_three_panel()

dev.off()

pdf(
  file = OUT_PANEL_PDF,
  width = 10,
  height = 8
)

draw_three_panel()

dev.off()

# ============================================================
# Print preview in R window
# ============================================================
# RStudio/R web plot pane sometimes gives:
# "Viewport has zero dimension(s)" when ggrepel is printed.
# So preview uses the same plots but removes ggrepel labels.
# Official saved PNG/PDF above keep the gene labels.

remove_ggrepel <- function(p) {
  keep <- !sapply(p$layers, function(x) {
    inherits(x$geom, "GeomTextRepel") || inherits(x$geom, "GeomLabelRepel")
  })
  p$layers <- p$layers[keep]
  p
}

pA_win <- remove_ggrepel(pA) + theme(legend.position = "none")
pB_win <- remove_ggrepel(pB) + theme(legend.position = "none")
pC_win <- pC

draw_three_panel_window <- function() {
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        widths = unit(c(1, 1), "null"),
        heights = unit(c(1, 0.85), "null")
      )
    )
  )
  
  print(
    pA_win,
    vp = viewport(layout.pos.row = 1, layout.pos.col = 1)
  )
  
  print(
    pB_win,
    vp = viewport(layout.pos.row = 1, layout.pos.col = 2)
  )
  
  print(
    pC_win,
    vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2)
  )
}

graphics.off()
draw_three_panel_window()

cat("\nDone. Files saved in:\n")
cat(OUTDIR, "\n")
cat("\nCheck with:\n")
cat("ls -lh ", OUTDIR, "\n", sep = "")