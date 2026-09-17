options(stringsAsFactors = FALSE)

root <- "/mnt/spareHD_2/genomewide_codeml_kuster"

alignment_file <- file.path(
    root,
    paste0(
        "05_gene_alignments/qa/",
        "gene_alignment_summary.with_class.tsv"
    )
)

input_qa_file <- file.path(
    root,
    "06_codeml/qa/codeml_input_QA.tsv"
)

result_file <- file.path(
    root,
    "07_codeml_genomewide/codeml_results.tsv"
)

site_file <- file.path(
    root,
    "07_codeml_genomewide/codeml_site_information.tsv"
)

old_label_file <- paste0(
    "/work/cyu/Kuster2026_genomewide_reanalysis/",
    "old287_vs_genomewide/",
    "old287_three_labels_codeml20347.tsv"
)

output_file <- file.path(
    root,
    "07_codeml_genomewide/",
    "codeml_master_analysis.tsv"
)

alignment <- read.delim(
    alignment_file,
    check.names = FALSE
)

input_qa <- read.delim(
    input_qa_file,
    check.names = FALSE
)

result <- read.delim(
    result_file,
    check.names = FALSE
)

sites <- read.delim(
    site_file,
    check.names = FALSE
)

old_labels <- read.delim(
    old_label_file,
    check.names = FALSE
)

# Main analysis gene universe: 17,965 adequately covered genes
master <- alignment[
    alignment$codeml_recommended == "yes",
]

if (nrow(master) != 17965) {
    stop(
        "Expected 17,965 main genes, found ",
        nrow(master)
    )
}

input_keep <- input_qa[, c(
    "gene_id",
    "variable_nt_sites",
    "variable_codon_sites",
    "codeml_input_status"
)]

master <- merge(
    master,
    input_keep,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

result_columns <- c(
    "gene_id",
    "exit_status",
    "elapsed_seconds",
    "lnL",
    "kappa",
    "omega",
    "tree_length_dN",
    "tree_length_dS"
)

master <- merge(
    master,
    result[, result_columns],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

site_columns <- c(
    "gene_id",
    "N_sites",
    "S_sites",
    "site_parse_status"
)

master <- merge(
    master,
    sites[, site_columns],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

numeric_columns <- c(
    "variable_nt_sites",
    "variable_codon_sites",
    "elapsed_seconds",
    "lnL",
    "kappa",
    "omega",
    "tree_length_dN",
    "tree_length_dS",
    "N_sites",
    "S_sites"
)

for (column in numeric_columns) {
    master[[column]] <- suppressWarnings(
        as.numeric(master[[column]])
    )
}

is_invariant <- (
    master$codeml_input_status == "invariant"
)

is_variable <- (
    master$codeml_input_status == "variable"
)

# Invariant genes have observed zero divergence but no identifiable omega.
master$tree_length_dN[is_invariant] <- 0
master$tree_length_dS[is_invariant] <- 0

master$rate_status <- NA_character_

master$rate_status[is_invariant] <- "dN0_dS0"

master$rate_status[
    is_variable &
    master$tree_length_dN <= 0 &
    master$tree_length_dS <= 0
] <- "dN0_dS0"

master$rate_status[
    is_variable &
    master$tree_length_dN > 0 &
    master$tree_length_dS <= 0
] <- "dNpositive_dS0"

master$rate_status[
    is_variable &
    master$tree_length_dN <= 0 &
    master$tree_length_dS > 0
] <- "dN0_dSpositive"

master$rate_status[
    is_variable &
    master$tree_length_dN > 0 &
    master$tree_length_dS > 0
] <- "dNpositive_dSpositive"

# Expected total nonsynonymous and synonymous changes across the tree.
master$expected_N_changes <- (
    master$N_sites *
    master$tree_length_dN
)

master$expected_S_changes <- (
    master$S_sites *
    master$tree_length_dS
)

# PAML boundaries are flags, not literal biological estimates.
master$omega_at_lower_boundary <- (
    !is.na(master$omega) &
    master$omega <= 0.0001001
)

master$omega_at_upper_boundary <- (
    !is.na(master$omega) &
    master$omega >= 998.999
)

master$omega_ES1 <- NA_real_
master$omega_ES2 <- NA_real_

usable_ES1 <- (
    is_variable &
    is.finite(master$omega) &
    master$tree_length_dS > 0 &
    master$expected_S_changes >= 1 &
    !master$omega_at_upper_boundary
)

usable_ES2 <- (
    usable_ES1 &
    master$expected_S_changes >= 2
)

master$omega_ES1[usable_ES1] <- master$omega[usable_ES1]
master$omega_ES2[usable_ES2] <- master$omega[usable_ES2]

# A lower-bound omega with dN rounded to zero represents no
# detectable nonsynonymous divergence.
master$omega_ES1[
    usable_ES1 &
    master$omega_at_lower_boundary
] <- 0

master$omega_ES2[
    usable_ES2 &
    master$omega_at_lower_boundary
] <- 0

master$omega_reliable_ES1 <- !is.na(master$omega_ES1)
master$omega_reliable_ES2 <- !is.na(master$omega_ES2)

# Add old functional labels for the original targeted gene set.
old_keep <- old_labels[, c(
    "gene_id",
    "old_name",
    "new_symbol",
    "own_role",
    "own_complex",
    "core_status"
)]

old_keep <- old_keep[
    !duplicated(old_keep$gene_id),
]

master <- merge(
    master,
    old_keep,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

master$old_target_gene <- ifelse(
    is.na(master$old_name),
    "no",
    "yes"
)

write.table(
    master,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

cat("===== Master dimensions =====\n")
cat("Genes:", nrow(master), "\n")
cat("Variable genes:", sum(is_variable), "\n")
cat("Invariant genes:", sum(is_invariant), "\n")
cat("Old target genes retained:", sum(master$old_target_gene=="yes"), "\n")

cat("\n===== Rate status =====\n")
print(table(master$rate_status))

cat("\n===== PAML omega boundaries =====\n")
cat(
    "Lower boundary:",
    sum(master$omega_at_lower_boundary, na.rm=TRUE),
    "\n"
)
cat(
    "Upper boundary:",
    sum(master$omega_at_upper_boundary, na.rm=TRUE),
    "\n"
)

cat("\n===== Reliable omega counts =====\n")
cat(
    "Expected S changes >=1:",
    sum(master$omega_reliable_ES1),
    "\n"
)
cat(
    "Expected S changes >=2:",
    sum(master$omega_reliable_ES2),
    "\n"
)

cat("\n===== Main genes by Kuster class =====\n")
print(table(master$Kuster_class))

cat("\n===== Kuster class x rate status =====\n")
print(table(
    master$Kuster_class,
    master$rate_status
))

cat("\n===== Median dN and dS, including invariant zeros =====\n")
print(aggregate(
    cbind(
        tree_length_dN,
        tree_length_dS
    ) ~ Kuster_class,
    data = master,
    FUN = median,
    na.rm = TRUE
))

cat("\n===== Median reliable omega: ES >= 1 =====\n")
print(tapply(
    master$omega_ES1,
    master$Kuster_class,
    median,
    na.rm = TRUE
))

cat("\n===== Old functional roles retained =====\n")
print(table(
    master$own_role,
    useNA = "ifany"
))

cat("\nOutput:", output_file, "\n")
