#figure3  panela-d

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# =========================
# Input files
# =========================
global_file <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_perSNP.tsv.gz"
ak_file <- "/mnt/spareHD_2/nu_287/q2_parallelism/within_region_strict_rawGene_min10/within_AK_perSNP.strict.tsv.gz"
bc_file <- "/mnt/spareHD_2/nu_287/q2_parallelism/within_region_strict_rawGene_min10/within_BC_perSNP.strict.tsv.gz"
pair_file <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_tables_gene_snp_pop/INTERSECT_sharedSNP_AK_BC_twoLines.tsv"

out_png <- "/mnt/spareHD_2/nu_287/q2_parallelism/Fig3_ABCD_final_pastel.png"
out_pdf <- "/mnt/spareHD_2/nu_287/q2_parallelism/Fig3_ABCD_final_pastel.pdf"

# =========================
# Parameters
# =========================
cut_dAF <- 0.5

# =========================
# Soft pastel colors
# =========================
col_shared <- "#5B7DB1"      # 柔和蓝
col_region <- "#BFC9D9"      # 浅灰蓝
col_high   <- "#E49A8D"      # 柔和珊瑚粉
col_bar    <- "#7F7F7F"
col_line   <- "#6A8DFF"

# =========================
# Theme
# =========================
theme_pub <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10, colour = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      plot.tag = element_text(face = "bold", size = 15),
      plot.margin = margin(8, 10, 8, 10)
    )
}

# =========================
# Panel A: all retained SNPs
# =========================
SNP <- fread(cmd = paste("zcat", shQuote(global_file)))
SNP <- SNP[is.finite(dAF_AK) & is.finite(dAF_BC)]

pearson <- cor.test(SNP$dAF_AK, SNP$dAF_BC, method = "pearson")

xA_min <- min(SNP$dAF_AK, na.rm = TRUE)
xA_max <- max(SNP$dAF_AK, na.rm = TRUE)
yA_min <- min(SNP$dAF_BC, na.rm = TRUE)
yA_max <- max(SNP$dAF_BC, na.rm = TRUE)

panelA <- ggplot(SNP, aes(dAF_AK, dAF_BC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey20") +
  geom_point(alpha = 0.18, size = 0.8, colour = "grey45") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, colour = col_line) +
  annotate(
    "text",
    x = xA_min + 0.05 * (xA_max - xA_min),
    y = yA_max - 0.05 * (yA_max - yA_min),
    hjust = 0, vjust = 1,
    label = sprintf(
      "n = %s\nPearson's r = %.2f\np < 2.2 x 10^-16",
      format(nrow(SNP), big.mark = ","),
      unname(pearson$estimate)
    ),
    size = 3.1
  ) +
  labs(
    x = expression(Delta*AF[AK]~"(AK freshwater mean - RS)"),
    y = expression(Delta*AF[BC]~"(BC freshwater mean - SAY)")
  ) +
  theme_pub()

# =========================
# Load within-region signed dAF
# =========================
AK <- fread(cmd = paste("zcat", shQuote(ak_file)))
BC <- fread(cmd = paste("zcat", shQuote(bc_file)))

setnames(AK, c("dAF_med", "concord_prop"), c("dAF_med_AK", "concord_AK"))
setnames(BC, c("dAF_med", "concord_prop"), c("dAF_med_BC", "concord_BC"))

AB <- merge(
  AK[, .(snp, dAF_med_AK, concord_AK)],
  BC[, .(snp, dAF_med_BC, concord_BC)],
  by = "snp"
)

AB[, c("chr2", "pos2") := tstrsplit(snp, ":", keep = c(1, 2))]
AB[, pos2 := as.numeric(pos2)]

# =========================
# Load 163 comparable shared SNP table
# =========================
PAIR <- fread(pair_file, header = FALSE)
setnames(PAIR, c(
  "region", "gene", "snp", "chr", "pos",
  "maj_sign", "concordance", "mean_abs_delta", "pops_nonzero"
))

PAIR_AK <- PAIR[region == "AK", .(
  snp,
  gene_pair = gene,
  maj_sign_AK = maj_sign
)]

PAIR_BC <- PAIR[region == "BC", .(
  snp,
  gene_pair = gene,
  maj_sign_BC = maj_sign
)]

COMP <- merge(PAIR_AK, PAIR_BC, by = "snp")
stopifnot(all(COMP$gene_pair.x == COMP$gene_pair.y))
COMP[, gene := gene_pair.x]

COMP[, c("chr2", "pos2") := tstrsplit(snp, ":", keep = c(1, 2))]
COMP[, pos2 := as.numeric(pos2)]

subset_dt <- merge(
  COMP,
  AB[, .(chr2, pos2, dAF_med_AK, dAF_med_BC, concord_AK, concord_BC)],
  by = c("chr2", "pos2"),
  all.x = TRUE
)

subset_dt[, concordant_dir := maj_sign_AK == maj_sign_BC]
subset_dt[, class := ifelse(concordant_dir, "Shared parallel", "Region-specific / discordant")]

subset_dt[, high_effect :=
  concordant_dir &
  !is.na(dAF_med_AK) &
  !is.na(dAF_med_BC) &
  abs(dAF_med_AK) >= cut_dAF &
  abs(dAF_med_BC) >= cut_dAF
]

n_subset   <- nrow(subset_dt)
n_parallel <- sum(subset_dt$concordant_dir, na.rm = TRUE)
n_other    <- n_subset - n_parallel
prop_parallel <- n_parallel / n_subset
binom_p <- binom.test(n_parallel, n_subset, p = 0.5)$p.value
n_high <- sum(subset_dt$high_effect, na.rm = TRUE)

# ranges for annotations
xB_min <- min(subset_dt$dAF_med_AK, na.rm = TRUE)
xB_max <- max(subset_dt$dAF_med_AK, na.rm = TRUE)
yB_min <- min(subset_dt$dAF_med_BC, na.rm = TRUE)
yB_max <- max(subset_dt$dAF_med_BC, na.rm = TRUE)

# =========================
# Panel B
# =========================
panelB <- ggplot(subset_dt, aes(dAF_med_AK, dAF_med_BC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey20") +
  geom_point(
    data = subset_dt[class == "Region-specific / discordant"],
    shape = 21, stroke = 0.5, size = 2.1,
    fill = "white", colour = col_region, alpha = 0.95
  ) +
  geom_point(
    data = subset_dt[class == "Shared parallel"],
    shape = 16, size = 2.1,
    colour = col_shared, alpha = 0.95
  ) +
  annotate(
    "text",
    x = xB_min + 0.04 * (xB_max - xB_min),
    y = yB_max - 0.03 * (yB_max - yB_min),
    hjust = 0, vjust = 1,
    label = sprintf(
      "Comparable SNPs = %d\nShared parallel = %d (%.1f%%)\nRegion-specific / discordant = %d (%.1f%%)",
      n_subset, n_parallel, 100 * prop_parallel,
      n_other, 100 * (n_other / n_subset)
    ),
    size = 2.7
  ) +
  labs(
    x = expression("AK within-region median "*Delta*AF),
    y = expression("BC within-region median "*Delta*AF)
  ) +
  theme_pub()

# =========================
# Panel C
# =========================
bar_dt <- data.table(
  category = c("Region-specific /\ndiscordant", "Shared parallel"),
  n = c(n_other, n_parallel)
)
bar_dt[, prop := n / sum(n)]
bar_dt[, label := sprintf("%d\n(%.1f%%)", n, 100 * prop)]

panelC <- ggplot(bar_dt, aes(category, n)) +
  geom_col(width = 0.65, fill = col_bar) +
  geom_text(aes(label = label), vjust = -0.2, size = 3.6) +
  expand_limits(y = max(bar_dt$n) * 1.22) +
  annotate(
    "text",
    x = 0.78,
    y = max(bar_dt$n) * 1.12,
    hjust = 0,
    label = "Binomial p = 1.3e-09",
    size = 3.1
  ) +
  labs(
    x = NULL,
    y = "Number of SNPs"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# =========================
# Panel D
# =========================
panelD <- ggplot(subset_dt, aes(dAF_med_AK, dAF_med_BC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = c(-cut_dAF, cut_dAF), linetype = "dotted", linewidth = 0.4, colour = "grey20") +
  geom_hline(yintercept = c(-cut_dAF, cut_dAF), linetype = "dotted", linewidth = 0.4, colour = "grey20") +
  geom_point(
    data = subset_dt[class == "Region-specific / discordant"],
    shape = 16, size = 1.8,
    colour = col_region, alpha = 0.7
  ) +
  geom_point(
    data = subset_dt[class == "Shared parallel"],
    shape = 16, size = 2.0,
    colour = col_shared, alpha = 0.9
  ) +
  geom_point(
    data = subset_dt[high_effect == TRUE],
    shape = 17, size = 3.2,
    colour = col_high, alpha = 1
  ) +
  annotate(
    "text",
    x = xB_min + 0.04 * (xB_max - xB_min),
    y = yB_max - 0.03 * (yB_max - yB_min),
    hjust = 0, vjust = 1,
    label = sprintf(
      "High-effect shared SNPs = %d\nabs(ΔAF) >= %.1f in both regions",
      n_high, cut_dAF
    ),
    size = 2.8
  ) +
  labs(
    x = expression("AK within-region median "*Delta*AF),
    y = expression("BC within-region median "*Delta*AF)
  ) +
  theme_pub()

# =========================
# Combine
# =========================
fig <- (panelA + panelB) / (panelC + panelD) +
  plot_annotation(tag_levels = "A")

print(fig)
ggsave(out_png, fig, width = 12, height = 9, dpi = 300)
ggsave(out_pdf, fig, width = 12, height = 9)