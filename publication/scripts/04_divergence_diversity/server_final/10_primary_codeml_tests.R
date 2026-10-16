
# Normalize a legacy classification header by its biological contents.
normalize_mitochondrial_class <- function(d) {
  if (!is.data.frame(d) || "mitochondrial_class" %in% names(d)) return(d)
  choices <- grep("_class$", names(d), value=TRUE, ignore.case=TRUE)
  choices <- setdiff(choices, "primary_class")
  valid <- c("direct_n-mt", "indirect_n-mt", "non-n-mt", "mtDNA",
             "direct n-mt", "indirect n-mt", "non n-mt")
  choices <- choices[vapply(choices, function(k) {
    v <- unique(as.character(d[[k]])); v <- v[!is.na(v) & nzchar(v)]
    length(v)>0L && all(v %in% valid)
  }, logical(1))]
  if (length(choices)>1L) stop("Ambiguous mitochondrial classification columns")
  if (length(choices)==1L) {
    if (inherits(d,"data.table")) data.table::setnames(d,choices,"mitochondrial_class")
    else names(d)[names(d)==choices] <- "mitochondrial_class"
  }
  d
}
read_analysis_delim <- function(...) normalize_mitochondrial_class(utils::read.delim(...))
read_analysis_table <- function(...) normalize_mitochondrial_class(data.table::fread(...))

options(stringsAsFactors = FALSE)

master_file <- paste0(
    "/path/to/data/genomewide_codeml/",
    "07_codeml_genomewide/",
    "codeml_master_analysis.tsv"
)

outdir <- paste0(
    "/path/to/data/genomewide_codeml/",
    "08_statistics"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

x <- read_analysis_delim(
    master_file,
    check.names = FALSE
)

metrics <- c(
    tree_length_dN = "dN",
    tree_length_dS = "dS",
    omega_ES1 = "omega_ES1",
    omega_ES2 = "omega_ES2"
)

comparisons <- list(
    direct_vs_non = list(
        column = "mitochondrial_class",
        group1 = "direct_n-mt",
        group2 = "non-n-mt",
        family = "mitochondrial_association"
    ),
    indirect_vs_non = list(
        column = "mitochondrial_class",
        group1 = "indirect_n-mt",
        group2 = "non-n-mt",
        family = "mitochondrial_association"
    ),
    direct_vs_indirect = list(
        column = "mitochondrial_class",
        group1 = "direct_n-mt",
        group2 = "indirect_n-mt",
        family = "mitochondrial_association"
    ),
    Nmt_ribo_vs_cyto_ribo = list(
        column = "own_role",
        group1 = "Nmt-ribo",
        group2 = "cyto-ribo",
        family = "functional"
    ),
    Nmt_ARS_vs_cyto_ARS = list(
        column = "own_role",
        group1 = "Nmt-ARS",
        group2 = "cyto-ARS",
        family = "functional"
    ),
    subunit_vs_assembly = list(
        column = "own_role",
        group1 = "subunit",
        group2 = "assembly_factor",
        family = "functional"
    ),
    core_vs_noncore = list(
        column = "core_status",
        group1 = "nu_core",
        group2 = "nu_noncore",
        family = "structural"
    )
)

wilcoxon_result <- function(
    data,
    column,
    group1,
    group2,
    metric
) {
    a <- data[
        data[[column]] == group1,
        metric
    ]

    b <- data[
        data[[column]] == group2,
        metric
    ]

    a <- a[is.finite(a)]
    b <- b[is.finite(b)]

    if (length(a) < 2 || length(b) < 2) {
        return(NULL)
    }

    test <- suppressWarnings(
        wilcox.test(
            a,
            b,
            alternative = "two.sided",
            exact = FALSE
        )
    )

    U <- unname(test$statistic)

    rank_biserial <- (
        2 * U / (length(a) * length(b))
    ) - 1

    data.frame(
        n_group1 = length(a),
        n_group2 = length(b),
        median_group1 = median(a),
        median_group2 = median(b),
        mean_group1 = mean(a),
        mean_group2 = mean(b),
        median_difference = median(a) - median(b),
        rank_biserial = rank_biserial,
        p = test$p.value,
        stringsAsFactors = FALSE
    )
}

results <- list()
index <- 1

for (comparison_name in names(comparisons)) {
    comparison <- comparisons[[comparison_name]]

    for (metric_column in names(metrics)) {
        result <- wilcoxon_result(
            data = x,
            column = comparison$column,
            group1 = comparison$group1,
            group2 = comparison$group2,
            metric = metric_column
        )

        if (is.null(result)) {
            next
        }

        result$contrast <- comparison_name
        result$family <- comparison$family
        result$group1 <- comparison$group1
        result$group2 <- comparison$group2
        result$metric_column <- metric_column
        result$metric <- metrics[[metric_column]]

        results[[index]] <- result
        index <- index + 1
    }
}

result_table <- do.call(
    rbind,
    results
)

result_table$p_bh_global <- p.adjust(
    result_table$p,
    method = "BH"
)

result_table$p_bh_within_metric <- ave(
    result_table$p,
    result_table$metric,
    FUN = function(p) p.adjust(p, method = "BH")
)

result_table <- result_table[, c(
    "family",
    "contrast",
    "group1",
    "group2",
    "metric",
    "n_group1",
    "n_group2",
    "median_group1",
    "median_group2",
    "median_difference",
    "mean_group1",
    "mean_group2",
    "rank_biserial",
    "p",
    "p_bh_within_metric",
    "p_bh_global"
)]

write.table(
    result_table,
    file.path(outdir, "primary_rate_comparisons.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Test whether the probability of observing dN=0 differs.
zero_results <- list()
index <- 1

for (comparison_name in names(comparisons)) {
    comparison <- comparisons[[comparison_name]]

    group <- x[[comparison$column]]

    selected <- (
        group %in% c(
            comparison$group1,
            comparison$group2
        )
    )

    table_data <- table(
        group[selected],
        x$tree_length_dN[selected] <= 0
    )

    if (
        nrow(table_data) != 2 ||
        ncol(table_data) != 2
    ) {
        next
    }

    fisher <- fisher.test(table_data)

    zero_results[[index]] <- data.frame(
        family = comparison$family,
        contrast = comparison_name,
        group1 = comparison$group1,
        group2 = comparison$group2,
        group1_zero_fraction = mean(
            x$tree_length_dN[
                group == comparison$group1
            ] <= 0,
            na.rm = TRUE
        ),
        group2_zero_fraction = mean(
            x$tree_length_dN[
                group == comparison$group2
            ] <= 0,
            na.rm = TRUE
        ),
        odds_ratio = unname(fisher$estimate),
        p = fisher$p.value,
        stringsAsFactors = FALSE
    )

    index <- index + 1
}

zero_table <- do.call(
    rbind,
    zero_results
)

zero_table$p_bh <- p.adjust(
    zero_table$p,
    method = "BH"
)

write.table(
    zero_table,
    file.path(outdir, "dN_zero_fraction_tests.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("===== Primary rate comparisons =====\n")

print(
    result_table[
        result_table$metric %in% c("dN", "dS", "omega_ES1"),
    ],
    row.names = FALSE
)

cat("\n===== dN zero-fraction comparisons =====\n")
print(
    zero_table,
    row.names = FALSE
)

cat("\n===== Outputs =====\n")
cat(
    file.path(outdir, "primary_rate_comparisons.tsv"),
    "\n"
)
cat(
    file.path(outdir, "dN_zero_fraction_tests.tsv"),
    "\n"
)
