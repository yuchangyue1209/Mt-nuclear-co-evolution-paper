#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "erc_helpers.R"))

set.seed(20260810)
erc_dir <- "/path/to/data/genomewide_codeml_kuster/10_erc_results/mtPCG_composite_t_spearman"
n_perm <- 100000L
nuclear <- read_rate_matrix(file.path(erc_dir, "nuclear_relative_rates.tsv.gz"))
mt <- read_rate_matrix(file.path(erc_dir, "mt13_relative_rates.tsv")); rownames(mt) <- toupper(sub("^MT[-_]", "", rownames(mt)))
branches <- intersect(colnames(nuclear), colnames(mt)); nuclear <- nuclear[, branches]; mt <- mt[, branches]
ann <- fread(file.path(erc_dir, "erc_mtPCG_composite_annotated.tsv")); setnames(ann, identify_gene_column(ann), "gene")

complexes <- c("CI", "CII", "CIII", "CIV", "CV")
nuclear_sets <- setNames(lapply(complexes, function(z) intersect(ann[own_role == "subunit" & own_complex == z, gene], rownames(nuclear))), complexes)
mt_sets <- list(CI = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"), CII = rownames(mt),
                CIII = "CYTB", CIV = c("COX1", "COX2", "COX3"), CV = c("ATP6", "ATP8"))
terminal <- vapply(branches, branch_is_terminal, logical(1)); terminal_i <- which(terminal); internal_i <- which(!terminal)

results <- list(); nulls <- list()
for (z in complexes) {
  nr <- composite_rate(nuclear, nuclear_sets[[z]]); mr <- composite_rate(mt, mt_sets[[z]])
  observed <- safe_spearman(nr, mr); unrestricted <- numeric(n_perm); stratified <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    unrestricted[i] <- safe_spearman(nr, sample(mr))
    permuted <- mr; permuted[terminal_i] <- sample(mr[terminal_i]); permuted[internal_i] <- sample(mr[internal_i])
    stratified[i] <- safe_spearman(nr, permuted)
  }
  results[[z]] <- data.table(complex = z, nuclear_n = length(nuclear_sets[[z]]), mt_n = length(mt_sets[[z]]),
    observed_rs = observed, unrestricted_p = (1 + sum(unrestricted >= observed))/(n_perm + 1),
    stratified_p = (1 + sum(stratified >= observed))/(n_perm + 1))
  nulls[[z]] <- data.table(iteration = seq_len(n_perm), complex = z,
                            unrestricted_rs = unrestricted, stratified_rs = stratified)
}
results <- rbindlist(results); results[, unrestricted_padj := p.adjust(unrestricted_p, "BH")]
results[, stratified_padj := p.adjust(stratified_p, "BH")]
fwrite(results, file.path(erc_dir, "OXPHOS_all_complex_permutation_summary.tsv"), sep = "\t")
fwrite(rbindlist(nulls), file.path(erc_dir, "OXPHOS_all_complex_permutation_null_100000.tsv.gz"), sep = "\t")

civ_n <- nuclear_sets$CIV; civ_m <- mt_sets$CIV
nr <- composite_rate(nuclear, civ_n); mr <- composite_rate(mt, civ_m)
loo_mt <- rbindlist(lapply(civ_m, function(g) data.table(excluded_mt_gene = g, rs = safe_spearman(nr, composite_rate(mt, setdiff(civ_m, g))))))
loo_n <- rbindlist(lapply(civ_n, function(g) data.table(excluded_nuclear_gene = g, rs = safe_spearman(composite_rate(nuclear, setdiff(civ_n, g)), mr))))
loo_branch <- rbindlist(lapply(seq_along(branches), function(i) data.table(excluded_branch = branches[i], rs = safe_spearman(nr[-i], mr[-i]))))
fwrite(loo_mt, file.path(erc_dir, "CIV_leave_one_mt_gene_out.tsv"), sep = "\t")
fwrite(loo_n, file.path(erc_dir, "CIV_leave_one_nuclear_gene_out.tsv"), sep = "\t")
fwrite(loo_branch, file.path(erc_dir, "CIV_leave_one_branch_out.tsv"), sep = "\t")
print(results)
cat("[OK] Complex permutations and CIV robustness completed\n")
