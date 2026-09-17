#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

erc_dir <- "/path/to/data/genomewide_codeml_kuster/10_erc_results/mtPCG_composite_t_spearman"
figure_dir <- "/path/to/data/genomewide_codeml_kuster/11_erc_figures/Figure_ERC_Weaver_style"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
x <- fread(file.path(erc_dir, "weaver_gene_set_erc_summary.tsv"))

red <- "#CA0029"; blue <- "#0877BA"
base_theme <- theme_classic(base_size = 17) + theme(
  axis.title = element_text(size = 17), axis.text = element_text(size = 15, colour = "black"),
  plot.title = element_text(size = 20, face = "bold"), legend.position = "bottom",
  legend.title = element_blank(), plot.margin = margin(8, 16, 8, 8))

functional_order <- c("Cyto_ARS_control", "Cyto_ribosomal_control", "OXPHOS_assembly_factors",
  "Nmt_ARS", "Nmt_ribosomal", "Nmt_OXPHOS_structural", "Indirect_nmt", "Direct_nmt")
functional_labels <- c(Cyto_ARS_control = "Cyto-ARS", Cyto_ribosomal_control = "Cyto-ribo",
  OXPHOS_assembly_factors = "OXPHOS assembly", Nmt_ARS = "Nmt-ARS", Nmt_ribosomal = "Nmt-ribo",
  Nmt_OXPHOS_structural = "N-mt OXPHOS", Indirect_nmt = "Indirect n-mt", Direct_nmt = "Direct n-mt")
a <- x[nuclear_set %in% functional_order]
a[, label := factor(functional_labels[nuclear_set], levels = functional_labels[functional_order])]
a[, significance := ifelse(is.finite(ci_lower) & ci_lower > 0, "Significant", "Not significant")]

panel_a <- ggplot(a, aes(bootstrap_mean_rs, label, colour = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, linewidth = 1) +
  geom_point(size = 4) +
  scale_colour_manual(values = c(Significant = red, `Not significant` = blue), drop = FALSE) +
  scale_x_continuous(limits = c(-.5, .65), breaks = seq(-.4, .6, .2)) +
  labs(title = "A", x = expression("mt–nuclear ERC (" * r[s] * ")"), y = NULL) + base_theme

complex_order <- paste0("OXPHOS_C", c("I", "II", "III", "IV", "V"))
complex_labels <- c(OXPHOS_CI = "Complex I", OXPHOS_CII = "Complex II†",
  OXPHOS_CIII = "Complex III", OXPHOS_CIV = "Complex IV", OXPHOS_CV = "Complex V")
b <- x[nuclear_set %in% complex_order]
b[, label := factor(complex_labels[nuclear_set], levels = rev(complex_labels[complex_order]))]
b[, significance := ifelse(is.finite(ci_lower) & ci_lower > 0, "Significant", "Not significant")]
b[, star := ifelse(significance == "Significant", "*", "")]

panel_b <- ggplot(b, aes(bootstrap_mean_rs, label, colour = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, linewidth = 1) +
  geom_point(size = 4) +
  geom_text(aes(x = ci_upper + .035, label = star), colour = "black", size = 7, vjust = .35) +
  scale_colour_manual(values = c(Significant = red, `Not significant` = blue), drop = FALSE) +
  scale_x_continuous(limits = c(-.45, .67), breaks = seq(-.4, .6, .2)) +
  labs(title = "B", x = expression("Within-complex ERC (" * r[s] * ")"), y = NULL) + base_theme

figure <- (panel_a | panel_b) + plot_layout(guides = "collect", widths = c(1.15, 1)) &
  theme(legend.position = "bottom")
pdf <- file.path(figure_dir, "Figure_ERC_Weaver_style_A_B_direct_indirect.pdf")
png <- file.path(figure_dir, "Figure_ERC_Weaver_style_A_B_direct_indirect.png")
ggsave(pdf, figure, width = 12.5, height = 6.2, device = cairo_pdf)
ggsave(png, figure, width = 12.5, height = 6.2, dpi = 400, bg = "white")
cat("[OK]", pdf, "\n[OK]", png, "\n")
