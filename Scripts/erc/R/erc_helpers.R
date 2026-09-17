suppressPackageStartupMessages(library(data.table))

read_rate_matrix <- function(path) {
  x <- fread(path)
  ids <- as.character(x[[1]])
  rate_columns <- names(x)[-1]
  x[, (rate_columns) := lapply(.SD, function(z) suppressWarnings(as.numeric(z))),
    .SDcols = rate_columns]
  result <- as.matrix(x[, ..rate_columns])
  storage.mode(result) <- "double"
  rownames(result) <- ids
  result[rowSums(is.finite(result)) >= 5, , drop = FALSE]
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5 || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

composite_rate <- function(rate_matrix, genes) {
  genes <- intersect(genes, rownames(rate_matrix))
  if (!length(genes)) stop("Empty gene set")
  colMeans(rate_matrix[genes, , drop = FALSE], na.rm = TRUE)
}

clean_label <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  x
}

identify_gene_column <- function(x) {
  candidates <- c("nuclear_gene", "gene_id", "stickleback_gene",
                  "stickleback_gene_id", "Gene stable ID", "gene")
  answer <- candidates[candidates %in% names(x)][1]
  if (is.na(answer)) stop("Could not identify gene-ID column")
  answer
}

branch_is_terminal <- function(label, n_tips = 27L) {
  hits <- regmatches(label, gregexpr("[0-9]+", label))[[1]]
  nodes <- suppressWarnings(as.integer(hits))
  if (length(nodes) < 2) return(NA)
  any(nodes <= n_tips)
}
