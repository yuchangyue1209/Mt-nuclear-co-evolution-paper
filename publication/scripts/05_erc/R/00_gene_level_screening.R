#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

erc_dir <- "/path/to/data/genomewide_codeml_kuster/10_erc_results/mtPCG_composite_t_spearman"
x <- fread(file.path(erc_dir, "erc_mtPCG_composite_annotated.tsv"))
x[, r := as.numeric(r)]
x <- x[is.finite(r)]

x[, class := gsub("_", "-", primary_class)]
x[class %in% c("direct-n-mt", "direct n-mt"), class := "direct_n-mt"]
x[class %in% c("indirect-n-mt", "indirect n-mt"), class := "indirect_n-mt"]
x[class %in% c("non-n-mt", "non n-mt"), class := "non-n-mt"]

background <- sort(x[class == "non-n-mt", r])
n_background <- length(background)
thresholds <- quantile(background, c(.90, .95, .975, .99), na.rm = TRUE)
names(thresholds) <- c("top10", "top05", "top025", "top01")

x[, background_percentile := findInterval(r, background) / n_background]
x[, empirical_p_upper := vapply(r, function(value)
  (1 + sum(background >= value)) / (n_background + 1), numeric(1))]
x[, empirical_padj := p.adjust(empirical_p_upper, "BH")]

for (label in names(thresholds))
  x[, (paste0("above_", label)) := r >= thresholds[[label]]]

fwrite(x[order(-r)], file.path(erc_dir, "erc_all_genes_empirical_statistics.tsv.gz"), sep = "\t")
fwrite(x[r >= thresholds[["top05"]]][order(-r)],
       file.path(erc_dir, "erc_candidates_background_top05.tsv"), sep = "\t")
fwrite(x[r >= thresholds[["top01"]]][order(-r)],
       file.path(erc_dir, "erc_candidates_background_top01.tsv"), sep = "\t")

categories <- rbindlist(list(
  x[, .(nuclear_gene, group = "Kuster", category = class)],
  x[!is.na(own_role), .(nuclear_gene, group = "functional_role", category = own_role)],
  x[!is.na(own_complex), .(nuclear_gene, group = "OXPHOS_complex", category = own_complex)],
  x[!is.na(core_status), .(nuclear_gene, group = "core_status", category = core_status)]
))

tests <- unique(categories[, .(group, category)])[, {
  genes <- categories[group == .BY$group & category == .BY$category, nuclear_gene]
  focal <- x[nuclear_gene %in% genes, r]
  wt <- suppressWarnings(wilcox.test(focal, background, alternative = "greater", exact = FALSE))
  auc <- as.numeric(wt$statistic) / (length(focal) * length(background))
  .(n = length(focal), median_r = median(focal), cliffs_delta = 2 * auc - 1, p = wt$p.value)
}, by = .(group, category)]
tests[, padj := p.adjust(p, "BH")]
fwrite(tests[order(padj)], file.path(erc_dir, "erc_category_distribution_tests.tsv"), sep = "\t")

cat("[OK] Gene-level empirical screening completed\n")
