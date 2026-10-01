#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: Rscript Figure5C_recent_vs_established.R INPUT_RESULTS_DIR OUTPUT_DIR')
input <- normalizePath(args[1], mustWork = TRUE)
outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
f <- file.path(input, 'Figure5B_plot_population_data_all_populations.tsv.gz')
if (!file.exists(f)) stop('Missing saved candidate population data: ', f)
d <- fread(f)
required <- c('gene', 'region', 'population_history', 'pop', 'deltaAF', 'snp', 'marine_pop')
if (!all(required %in% names(d))) stop('Missing required data columns')
if (anyDuplicated(d[, .(gene, region, pop)])) stop('Duplicate gene-region-population observations')
if (any(!is.finite(d$deltaAF))) stop('Non-finite deltaAF values in saved data')
if (!all(d$population_history %in% c('established', 'recent'))) stop('Unexpected population history')
audit <- d[, .(n_snp = uniqueN(snp), n_marine_reference = uniqueN(marine_pop)), by = .(gene, region)]
if (any(audit$n_snp != 1L | audit$n_marine_reference != 1L)) stop('SNP or marine reference differs within gene-region')
summ <- d[, .(median_deltaAF = median(deltaAF), n_populations = uniqueN(pop)),
          by = .(gene, region, population_history)]
wide <- dcast(summ, gene + region ~ population_history, value.var = 'median_deltaAF')
if (!all(c('established', 'recent') %in% names(wide))) stop('Missing established or recent group')
if (any(!is.finite(wide$established) | !is.finite(wide$recent))) stop('Some gene-regions lack one population group')
counts <- dcast(summ, gene + region ~ population_history, value.var = 'n_populations')
setnames(counts, c('established', 'recent'), c('n_established', 'n_recent'))
wide <- merge(wide, counts, by = c('gene', 'region'))
wide[, same_direction := sign(established) == sign(recent) & sign(established) != 0]
wide[, recent_minus_established := recent - established]
fwrite(wide, file.path(outdir, 'Figure5C_plot_data.tsv'), sep = '\t')

# Descriptive correlations across the preselected candidate genes.
# These are region-level gene summaries, not correlations for individual populations.
cor_stats <- wide[, {
  pe <- cor.test(established, recent, method = 'pearson')
  sp <- cor.test(established, recent, method = 'spearman', exact = FALSE)
  list(n_genes = .N, pearson_r = unname(pe$estimate), pearson_P = pe$p.value,
       spearman_rho = unname(sp$estimate), spearman_P = sp$p.value,
       same_direction_genes = sum(same_direction),
       smaller_recent_shift_genes = sum(abs(recent) < abs(established)))
}, by = region]
fwrite(cor_stats, file.path(outdir, 'Figure5C_correlation_statistics.tsv'), sep = '\t')


# Same limits and physical scale on both axes.
lim <- c(-1.05, 1.05)
p <- ggplot(wide, aes(established, recent)) +
  geom_hline(yintercept = 0, colour = 'grey65', linewidth = 0.6, linetype = 'dashed') +
  geom_vline(xintercept = 0, colour = 'grey65', linewidth = 0.6, linetype = 'dashed') +
  geom_abline(intercept = 0, slope = 1, colour = 'grey40', linewidth = 0.8, linetype = 'dotted') +
  geom_point(aes(shape = region), colour = '#009E73', size = 5, stroke = 1.1) +
  geom_text_repel(aes(label = gene), size = 7.5, fontface = 'italic',
    box.padding = 0.7, point.padding = 0.55, force = 3,
    min.segment.length = 0, segment.colour = 'grey50', segment.size = 0.45,
    max.overlaps = Inf, max.time = 10, max.iter = 50000, seed = 105,
    xlim = c(-1, 1), ylim = c(-1, 1)) +
  scale_shape_manual(values = c(AK = 16, BC = 1), breaks = c('AK', 'BC'),
                     labels = c('Alaska', 'British Columbia')) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.5)) +
  coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = FALSE) +
  labs(title = 'C', x = expression(Established~median~Delta*AF),
       y = expression(Recent~median~Delta*AF), shape = NULL) +
  theme_classic(base_size = 28) +
  theme(axis.title = element_text(size = 32, colour = 'black'),
        axis.text = element_text(size = 26, colour = 'black'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.8),
        axis.ticks = element_line(colour = 'black', linewidth = 0.7),
        plot.title = element_text(size = 36, face = 'bold', hjust = 0),
        legend.position = 'top', legend.text = element_text(size = 26),
        legend.key.size = grid::unit(0.9, 'cm'),
        plot.margin = margin(16, 20, 16, 16))

ggsave(file.path(outdir, 'Figure5C_recent_vs_established.pdf'), p,
       width = 10, height = 10, units = 'in', device = grDevices::cairo_pdf, bg = 'white')
ggsave(file.path(outdir, 'Figure5C_recent_vs_established.png'), p,
       width = 10, height = 10, units = 'in', dpi = 300, bg = 'white')
cat('[OK]', nrow(wide), 'gene-region points; outputs saved to', outdir, '\n')
