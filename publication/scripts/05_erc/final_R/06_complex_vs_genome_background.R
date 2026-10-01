#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(20260929)

erc_dir <- paste0(
  "/path/to/data/genomewide_codeml_kuster/",
  "10_erc_results/mtPCG_composite_t_spearman"
)

out_dir <- file.path(
  erc_dir,
  "complex_vs_genome_background"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_permutations <- as.integer(Sys.getenv("N_PERM", "100000"))
batch_size <- as.integer(Sys.getenv("BATCH_SIZE", "1000"))

if (!is.finite(n_permutations) || n_permutations < 1000L) {
  stop("N_PERM must be an integer >= 1000")
}

if (!is.finite(batch_size) || batch_size < 1L) {
  stop("BATCH_SIZE must be a positive integer")
}

read_rate_matrix <- function(path) {
  dat <- fread(path)
  id_candidates <- c(
    "nuclear_gene", "gene_id", "gene", "mt_gene", "Gene", names(dat)[1]
  )
  id_column <- id_candidates[id_candidates %in% names(dat)][1]
  if (is.na(id_column)) stop("Cannot identify gene-ID column in ", path)
  ids <- as.character(dat[[id_column]])
  rate_columns <- setdiff(names(dat), id_column)
  for (column_name in rate_columns) {
    set(dat, j = column_name, value = suppressWarnings(as.numeric(dat[[column_name]])))
  }
  result <- as.matrix(dat[, ..rate_columns])
  storage.mode(result) <- "double"
  rownames(result) <- ids
  result[rowSums(is.finite(result)) > 1, , drop = FALSE]
}

clean_label <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a", ".")] <- NA_character_
  x
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L) return(NA_real_)
  if (sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

composite_rate <- function(rate_matrix, genes) {
  genes <- intersect(genes, rownames(rate_matrix))
  if (length(genes) == 0L) return(rep(NA_real_, ncol(rate_matrix)))
  colMeans(rate_matrix[genes, , drop = FALSE], na.rm = TRUE)
}

nuclear <- read_rate_matrix(file.path(erc_dir, "nuclear_relative_rates.tsv.gz"))
mt <- read_rate_matrix(file.path(erc_dir, "mt13_relative_rates.tsv"))
rownames(mt) <- toupper(sub("^MT[-_]", "", rownames(mt)))

shared_branches <- intersect(colnames(nuclear), colnames(mt))
nuclear <- nuclear[, shared_branches, drop = FALSE]
mt <- mt[, shared_branches, drop = FALSE]

annotation <- fread(file.path(erc_dir, "erc_mtPCG_composite_annotated.tsv"))
gene_candidates <- c(
  "analysis_gene_id", "nuclear_gene", "stickleback_gene", "gene_id", "gene"
)
gene_column <- gene_candidates[gene_candidates %in% names(annotation)][1]
if (is.na(gene_column)) stop("Cannot identify gene-ID column in ERC annotation")
setnames(annotation, gene_column, "analysis_gene_id")
annotation <- unique(annotation, by = "analysis_gene_id")

for (column_name in c("primary_class", "own_role", "own_complex")) {
  if (!column_name %in% names(annotation)) annotation[, (column_name) := NA_character_]
  annotation[, (column_name) := clean_label(get(column_name))]
}

annotation[
  primary_class %in% c("non_n_mt", "non_n-mt", "non-n-mt", "non n-mt"),
  primary_class := "non-n-mt"
]

mt_complexes <- list(
  CI = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"),
  CIII = "CYTB",
  CIV = c("COX1", "COX2", "COX3"),
  CV = c("ATP6", "ATP8")
)

complex_labels <- c(
  CI = "Complex I",
  CIII = "Complex III",
  CIV = "Complex IV",
  CV = "Complex V"
)

get_nuclear_complex <- function(complex_name) {
  selected <- annotation[
    own_role == "subunit" & own_complex == complex_name,
    analysis_gene_id
  ]
  intersect(selected, rownames(nuclear))
}

non_nmt_background <- intersect(
  annotation[primary_class == "non-n-mt", analysis_gene_id],
  rownames(nuclear)
)

if (length(non_nmt_background) == 0L) {
  stop("The non-n-mt background is empty after annotation matching")
}

draw_null <- function(
  background_genes,
  focal_genes,
  mt_genes,
  iterations,
  batch_size
) {
  focal_genes <- intersect(focal_genes, rownames(nuclear))
  mt_genes <- intersect(mt_genes, rownames(mt))
  background_genes <- setdiff(
    intersect(background_genes, rownames(nuclear)),
    focal_genes
  )

  group_size <- length(focal_genes)
  if (group_size == 0L) stop("Focal nuclear gene set is empty")
  if (length(mt_genes) == 0L) stop("Mitochondrial gene set is empty")
  if (length(background_genes) < group_size) {
    stop("Background has fewer genes than the focal set")
  }

  mt_composite <- composite_rate(mt, mt_genes)
  result <- rep(NA_real_, iterations)

  starts <- seq.int(1L, iterations, by = batch_size)
  for (start in starts) {
    stop_at <- min(start + batch_size - 1L, iterations)
    for (iteration in start:stop_at) {
      sampled_genes <- sample(
        background_genes,
        size = group_size,
        replace = FALSE
      )
      result[iteration] <- safe_spearman(
        composite_rate(nuclear, sampled_genes),
        mt_composite
      )
    }
  }

  result
}

backgrounds <- list(
  Non_nmt = non_nmt_background,
  Genomewide_nuclear = rownames(nuclear)
)

summary_rows <- list()
null_rows <- list()
counter <- 0L

cat("[load] nuclear genes =", nrow(nuclear), "branches =", ncol(nuclear), "\n")
cat("[load] non-n-mt genes =", length(non_nmt_background), "\n")
cat("[run] permutations per comparison =", n_permutations, "\n")

for (complex_name in names(mt_complexes)) {
  focal_genes <- get_nuclear_complex(complex_name)
  mitochondrial_genes <- intersect(mt_complexes[[complex_name]], rownames(mt))
  observed_rs <- safe_spearman(
    composite_rate(nuclear, focal_genes),
    composite_rate(mt, mitochondrial_genes)
  )

  for (background_name in names(backgrounds)) {
    counter <- counter + 1L
    cat(
      "[permutation]", complex_labels[[complex_name]], "vs", background_name,
      "nuclear =", length(focal_genes), "mt =", length(mitochondrial_genes), "\n"
    )

    null_values <- draw_null(
      backgrounds[[background_name]],
      focal_genes,
      mitochondrial_genes,
      n_permutations,
      batch_size
    )
    finite_null <- null_values[is.finite(null_values)]
    empirical_p <- (1 + sum(finite_null >= observed_rs)) /
      (length(finite_null) + 1)

    summary_rows[[counter]] <- data.table(
      complex = complex_name,
      complex_label = complex_labels[[complex_name]],
      background = background_name,
      nuclear_n = length(focal_genes),
      mt_n = length(mitochondrial_genes),
      nuclear_genes = paste(focal_genes, collapse = ","),
      mt_genes = paste(mitochondrial_genes, collapse = ","),
      observed_rs = observed_rs,
      null_mean_rs = mean(finite_null),
      null_median_rs = median(finite_null),
      null_ci_lower = quantile(finite_null, 0.025, names = FALSE),
      null_ci_upper = quantile(finite_null, 0.975, names = FALSE),
      empirical_p_greater = empirical_p,
      n_permutations = length(finite_null)
    )

    null_rows[[counter]] <- data.table(
      complex = complex_name,
      complex_label = complex_labels[[complex_name]],
      background = background_name,
      iteration = seq_along(null_values),
      null_rs = null_values
    )
  }
}

summary_table <- rbindlist(summary_rows)
summary_table[
  ,
  empirical_padj_greater := p.adjust(empirical_p_greater, method = "BH"),
  by = background
]
summary_table[
  ,
  significance := fifelse(
    empirical_padj_greater < 0.001, "***",
    fifelse(
      empirical_padj_greater < 0.01, "**",
      fifelse(empirical_padj_greater < 0.05, "*", "ns")
    )
  )
]

null_table <- rbindlist(null_rows)

summary_file <- file.path(
  out_dir,
  "complex_vs_genome_background_summary.tsv"
)
null_file <- file.path(
  out_dir,
  "complex_vs_genome_background_null_distributions.tsv.gz"
)

fwrite(summary_table, summary_file, sep = "\t")
fwrite(null_table, null_file, sep = "\t", compress = "gzip")

plot_data <- copy(null_table)
plot_data[, complex_label := factor(
  complex_label,
  levels = unname(complex_labels)
)]
plot_data[, background_label := fifelse(
  background == "Non_nmt",
  "Non-n-mt background",
  "Genome-wide nuclear background"
)]

observed_data <- copy(summary_table)
observed_data[, complex_label := factor(
  complex_label,
  levels = unname(complex_labels)
)]
observed_data[, background_label := fifelse(
  background == "Non_nmt",
  "Non-n-mt background",
  "Genome-wide nuclear background"
)]

p <- ggplot(plot_data, aes(x = complex_label, y = null_rs)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_violin(
    fill = "grey85",
    colour = "grey45",
    linewidth = 0.45,
    width = 0.85,
    trim = TRUE
  ) +
  geom_point(
    data = observed_data,
    aes(y = observed_rs),
    shape = 23,
    size = 4.2,
    stroke = 0.8,
    fill = "#CA0029",
    colour = "black"
  ) +
  geom_text(
    data = observed_data,
    aes(
      y = pmin(0.72, pmax(observed_rs, null_ci_upper) + 0.07),
      label = significance
    ),
    size = 5,
    fontface = "bold"
  ) +
  facet_wrap(~background_label, ncol = 1) +
  coord_cartesian(ylim = c(-0.65, 0.75), clip = "off") +
  labs(
    x = NULL,
    y = expression("Within-complex ERC (" * r[s] * ")")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16),
    strip.background = element_rect(fill = "grey92", colour = "black"),
    strip.text = element_text(size = 15, face = "bold"),
    plot.margin = margin(10, 18, 10, 10)
  )

pdf_file <- file.path(out_dir, "complex_vs_genome_background.pdf")
png_file <- file.path(out_dir, "complex_vs_genome_background.png")
ggsave(pdf_file, p, width = 8.5, height = 8.5, device = cairo_pdf)
ggsave(png_file, p, width = 8.5, height = 8.5, dpi = 400, bg = "white")

cat("\n===== Complex-versus-background ERC tests =====\n")
print(summary_table)
cat("\n[output]", summary_file, "\n")
cat("[output]", null_file, "\n")
cat("[output]", pdf_file, "\n")
cat("[output]", png_file, "\n")
cat("[note] Complex II was not tested because it has no mtDNA-encoded subunits.\n")
