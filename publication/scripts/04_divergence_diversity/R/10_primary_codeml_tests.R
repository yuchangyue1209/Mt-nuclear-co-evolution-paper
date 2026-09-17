#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
root <- "/path/to/data/genomewide_codeml_kuster"
input <- file.path(root, "07_codeml_genomewide/codeml_master_analysis.tsv")
outdir <- file.path(root, "08_statistics")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
x <- read.delim(input, check.names=FALSE)

comparisons <- list(
  direct_vs_non=c("Kuster_class", "direct_n-mt", "non-n-mt"),
  indirect_vs_non=c("Kuster_class", "indirect_n-mt", "non-n-mt"),
  direct_vs_indirect=c("Kuster_class", "direct_n-mt", "indirect_n-mt"),
  Nmt_ribo_vs_cyto_ribo=c("own_role", "Nmt-ribo", "cyto-ribo"),
  Nmt_ARS_vs_cyto_ARS=c("own_role", "Nmt-ARS", "cyto-ARS"),
  subunit_vs_assembly=c("own_role", "subunit", "assembly_factor"),
  core_vs_noncore=c("core_status", "nu_core", "nu_noncore")
)
metrics <- c("tree_length_dN", "tree_length_dS", "omega_ES1", "omega_ES2")

rows <- list(); k <- 1
for (contrast in names(comparisons)) {
  spec <- comparisons[[contrast]]; column <- spec[1]; g1 <- spec[2]; g2 <- spec[3]
  for (metric in metrics) {
    a <- x[[metric]][!is.na(x[[column]]) & x[[column]] == g1]; a <- a[is.finite(a)]
    b <- x[[metric]][!is.na(x[[column]]) & x[[column]] == g2]; b <- b[is.finite(b)]
    test <- suppressWarnings(wilcox.test(a, b, exact=FALSE)); U <- unname(test$statistic)
    rows[[k]] <- data.frame(contrast, group1=g1, group2=g2, metric, n1=length(a), n2=length(b), median1=median(a), median2=median(b), mean1=mean(a), mean2=mean(b), rank_biserial=2*U/(length(a)*length(b))-1, p=test$p.value)
    k <- k+1
  }
}
results <- do.call(rbind, rows)
results$p_bh_global <- p.adjust(results$p, "BH")
results$p_bh_within_metric <- ave(results$p, results$metric, FUN=function(p) p.adjust(p, "BH"))
write.table(results, file.path(outdir, "primary_rate_comparisons.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(results, row.names=FALSE)

