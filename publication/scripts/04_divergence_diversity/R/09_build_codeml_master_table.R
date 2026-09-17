#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- "/path/to/data/genomewide_codeml_kuster"
alignment_file <- file.path(root, "05_gene_alignments/qa/gene_alignment_summary.with_class.tsv")
input_qa_file <- file.path(root, "06_codeml/qa/codeml_input_QA.tsv")
result_file <- file.path(root, "07_codeml_genomewide/codeml_results.tsv")
site_file <- file.path(root, "07_codeml_genomewide/codeml_site_information.tsv")
old_label_file <- "/path/to/workspace/Kuster2026_genomewide_reanalysis/old287_vs_genomewide/old287_three_labels_codeml20347.tsv"
output_file <- file.path(root, "07_codeml_genomewide/codeml_master_analysis.tsv")

read_tsv <- function(path) read.delim(path, check.names=FALSE, quote="", comment.char="")
alignment <- read_tsv(alignment_file)
input_qa <- read_tsv(input_qa_file)
result <- read_tsv(result_file)
sites <- read_tsv(site_file)
labels <- read_tsv(old_label_file)

master <- alignment[alignment$codeml_recommended == "yes", ]
stopifnot(nrow(master) == 17965)

master <- merge(master, input_qa[, c("gene_id", "variable_nt_sites", "variable_codon_sites", "codeml_input_status")], by="gene_id", all.x=TRUE, sort=FALSE)
master <- merge(master, result[, c("gene_id", "exit_status", "elapsed_seconds", "lnL", "kappa", "omega", "tree_length_dN", "tree_length_dS")], by="gene_id", all.x=TRUE, sort=FALSE)
master <- merge(master, sites[, c("gene_id", "N_sites", "S_sites", "site_parse_status")], by="gene_id", all.x=TRUE, sort=FALSE)

numeric_columns <- c("variable_nt_sites", "variable_codon_sites", "elapsed_seconds", "lnL", "kappa", "omega", "tree_length_dN", "tree_length_dS", "N_sites", "S_sites")
for (column in numeric_columns) master[[column]] <- suppressWarnings(as.numeric(master[[column]]))

invariant <- master$codeml_input_status == "invariant"
variable <- master$codeml_input_status == "variable"
master$tree_length_dN[invariant] <- 0
master$tree_length_dS[invariant] <- 0

master$rate_status <- NA_character_
master$rate_status[invariant | (variable & master$tree_length_dN <= 0 & master$tree_length_dS <= 0)] <- "dN0_dS0"
master$rate_status[variable & master$tree_length_dN > 0 & master$tree_length_dS <= 0] <- "dNpositive_dS0"
master$rate_status[variable & master$tree_length_dN <= 0 & master$tree_length_dS > 0] <- "dN0_dSpositive"
master$rate_status[variable & master$tree_length_dN > 0 & master$tree_length_dS > 0] <- "dNpositive_dSpositive"

master$expected_N_changes <- master$N_sites * master$tree_length_dN
master$expected_S_changes <- master$S_sites * master$tree_length_dS
master$omega_at_lower_boundary <- !is.na(master$omega) & master$omega <= 0.0001001
master$omega_at_upper_boundary <- !is.na(master$omega) & master$omega >= 998.999

usable1 <- variable & is.finite(master$omega) & master$tree_length_dS > 0 & master$expected_S_changes >= 1 & !master$omega_at_upper_boundary
usable2 <- usable1 & master$expected_S_changes >= 2
master$omega_ES1 <- master$omega_ES2 <- NA_real_
master$omega_ES1[usable1] <- master$omega[usable1]
master$omega_ES2[usable2] <- master$omega[usable2]
master$omega_ES1[usable1 & master$omega_at_lower_boundary] <- 0
master$omega_ES2[usable2 & master$omega_at_lower_boundary] <- 0
master$omega_reliable_ES1 <- !is.na(master$omega_ES1)
master$omega_reliable_ES2 <- !is.na(master$omega_ES2)

label_columns <- c("gene_id", "old_name", "new_symbol", "own_role", "own_complex", "core_status")
labels <- labels[!duplicated(labels$gene_id), label_columns]
master <- merge(master, labels, by="gene_id", all.x=TRUE, sort=FALSE)
master$old_target_gene <- ifelse(is.na(master$old_name), "no", "yes")

write.table(master, output_file, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
cat("genes", nrow(master), "\n")
cat("variable", sum(variable), "\n")
cat("invariant", sum(invariant), "\n")
cat("omega_ES1", sum(master$omega_reliable_ES1), "\n")
cat("output", output_file, "\n")

