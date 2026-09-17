#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "erc_helpers.R"))

set.seed(20260810)
erc_dir <- "/path/to/data/genomewide_codeml_kuster/10_erc_results/mtPCG_composite_t_spearman"
n_boot <- 10000L

nuclear <- read_rate_matrix(file.path(erc_dir, "nuclear_relative_rates.tsv.gz"))
mt <- read_rate_matrix(file.path(erc_dir, "mt13_relative_rates.tsv"))
rownames(mt) <- toupper(sub("^MT[-_]", "", rownames(mt)))
branches <- intersect(colnames(nuclear), colnames(mt))
nuclear <- nuclear[, branches, drop = FALSE]
mt <- mt[, branches, drop = FALSE]

annotation <- fread(file.path(erc_dir, "erc_mtPCG_composite_annotated.tsv"))
setnames(annotation, identify_gene_column(annotation), "analysis_gene_id")
annotation <- unique(annotation, by = "analysis_gene_id")
for (column in c("primary_class", "own_role", "own_complex", "core_status")) {
  if (!column %in% names(annotation)) annotation[, (column) := NA_character_]
  annotation[, (column) := clean_label(get(column))]
}
annotation[primary_class %in% c("direct-n-mt", "direct_n-mt", "direct n-mt"), primary_class := "direct_n-mt"]
annotation[primary_class %in% c("indirect-n-mt", "indirect_n-mt", "indirect n-mt"), primary_class := "indirect_n-mt"]

genes <- function(role = NULL, complex = NULL, class = NULL, core = NULL) {
  keep <- rep(TRUE, nrow(annotation))
  if (!is.null(role)) keep <- keep & annotation$own_role == role
  if (!is.null(complex)) keep <- keep & annotation$own_complex == complex
  if (!is.null(class)) keep <- keep & annotation$primary_class == class
  if (!is.null(core)) keep <- keep & annotation$core_status == core
  intersect(annotation$analysis_gene_id[which(keep %in% TRUE)], rownames(nuclear))
}

sets <- list(
  Nmt_OXPHOS_structural = genes(role = "subunit"),
  OXPHOS_assembly_factors = genes(role = "assembly_factor"),
  Nmt_ribosomal = genes(role = "Nmt-ribo"), Nmt_ARS = genes(role = "Nmt-ARS"),
  Cyto_ribosomal_control = genes(role = "cyto-ribo"), Cyto_ARS_control = genes(role = "cyto-ARS"),
  Direct_nmt = genes(class = "direct_n-mt"), Indirect_nmt = genes(class = "indirect_n-mt"),
  OXPHOS_core = genes(core = "nu_core"), OXPHOS_noncore = genes(core = "nu_noncore"),
  OXPHOS_CI = genes(role = "subunit", complex = "CI"),
  OXPHOS_CII = genes(role = "subunit", complex = "CII"),
  OXPHOS_CIII = genes(role = "subunit", complex = "CIII"),
  OXPHOS_CIV = genes(role = "subunit", complex = "CIV"),
  OXPHOS_CV = genes(role = "subunit", complex = "CV")
)
mt_sets <- list(
  all = rownames(mt), CI = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"),
  CII = rownames(mt), CIII = "CYTB", CIV = c("COX1", "COX2", "COX3"), CV = c("ATP6", "ATP8")
)

definitions <- data.table(
  analysis = c(rep("functional_set", 10), "within_complex_control", rep("within_complex", 4)),
  nuclear_set = names(sets),
  mt_set = c(rep("all", 10), "CI", "CII", "CIII", "CIV", "CV")
)
# Correct the final five labels to CI, CII, CIII, CIV and CV.
definitions[11:15, mt_set := c("CI", "CII", "CIII", "CIV", "CV")]
definitions[11:15, analysis := c("within_complex", "within_complex_control", rep("within_complex", 3))]

boot_one <- function(nuclear_genes, mt_genes) {
  nuclear_genes <- intersect(nuclear_genes, rownames(nuclear)); mt_genes <- intersect(mt_genes, rownames(mt))
  vapply(seq_len(n_boot), function(i) safe_spearman(
    composite_rate(nuclear, sample(nuclear_genes, length(nuclear_genes), TRUE)),
    composite_rate(mt, sample(mt_genes, length(mt_genes), TRUE))), numeric(1))
}

summaries <- list(); bootstraps <- list()
for (i in seq_len(nrow(definitions))) {
  d <- definitions[i]; ng <- sets[[d$nuclear_set]]; mg <- mt_sets[[d$mt_set]]
  observed <- safe_spearman(composite_rate(nuclear, ng), composite_rate(mt, mg))
  b <- boot_one(ng, mg); finite <- b[is.finite(b)]
  summaries[[i]] <- data.table(analysis = d$analysis, nuclear_set = d$nuclear_set,
    nuclear_n = length(ng), mt_n = length(mg), observed_rs = observed,
    bootstrap_mean_rs = mean(finite), bootstrap_median_rs = median(finite),
    ci_lower = quantile(finite, .025), ci_upper = quantile(finite, .975),
    bootstrap_p_positive = (1 + sum(finite <= 0)) / (length(finite) + 1))
  bootstraps[[i]] <- data.table(iteration = seq_len(n_boot), analysis = d$analysis,
    nuclear_set = d$nuclear_set, rs = b)
}
summary_table <- rbindlist(summaries)
summary_table[, bootstrap_padj_positive := p.adjust(bootstrap_p_positive, "BH")]
fwrite(summary_table, file.path(erc_dir, "weaver_gene_set_erc_summary.tsv"), sep = "\t")
fwrite(rbindlist(bootstraps), file.path(erc_dir, "weaver_gene_set_erc_bootstrap_10000.tsv.gz"), sep = "\t")

composites <- rbindlist(lapply(seq_len(nrow(definitions)), function(i) {
  d <- definitions[i]
  data.table(analysis = d$analysis, nuclear_set = d$nuclear_set, branch = branches,
    nuclear_composite_rer = composite_rate(nuclear, sets[[d$nuclear_set]]),
    mt_composite_rer = composite_rate(mt, mt_sets[[d$mt_set]]))
}))
fwrite(composites, file.path(erc_dir, "weaver_gene_set_branch_composites.tsv"), sep = "\t")
print(summary_table)
cat("[OK] Weaver-style gene-set ERC completed\n")
