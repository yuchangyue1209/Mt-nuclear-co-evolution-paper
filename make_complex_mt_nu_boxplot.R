#!/usr/bin/env Rscript

#' Complex MT/NU Boxplot Generator
#' 
#' This script creates boxplots comparing omega (dN/dS) values between mitochondrial (mt)
#' and nuclear (nu) genes across different OXPHOS complexes (I-V).
#' 
#' Usage: Rscript make_complex_mt_nu_boxplot.R <dnds_from_mlc.tsv> <nu_bed> <out_png>
#' 
#' @param dnds_from_mlc.tsv Input file with dN/dS data from MLC analysis
#' @param nu_bed BED file with nuclear OXPHOS gene annotations including complex information
#' @param out_png Output PNG file path for the generated plot
#' 
#' @author Your Name
#' @date Sys.Date()

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

#' Display usage information and exit
show_usage <- function() {
  cat("Usage: Rscript make_complex_mt_nu_boxplot.R <dnds_from_mlc.tsv> <nu_bed> <out_png>\n")
  cat("\nArguments:\n")
  cat("  dnds_from_mlc.tsv  Input file with dN/dS data from MLC analysis\n")
  cat("  nu_bed            BED file with nuclear OXPHOS gene annotations\n")
  cat("  out_png           Output PNG file path\n")
  cat("\nExample:\n")
  cat("  Rscript make_complex_mt_nu_boxplot.R dnds_from_mlc.tsv nuOXPHOS_genes_with_complex_core.bed Fig_complex_mt_nu.png\n")
}

# Validate command line arguments
if (length(args) < 3) {
  cat("Error: Insufficient arguments provided.\n\n")
  show_usage()
  quit(status = 1)
}

fn_mlc <- args[1]
fn_bed <- args[2]
fn_out <- args[3]

# Validate input files exist
if (!file.exists(fn_mlc)) {
  stop("Error: Input file '", fn_mlc, "' does not exist.")
}
if (!file.exists(fn_bed)) {
  stop("Error: BED file '", fn_bed, "' does not exist.")
}

# Check if output directory exists, create if necessary
out_dir <- dirname(fn_out)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("Created output directory:", out_dir, "\n")
}

# Load required packages with better error handling
suppressPackageStartupMessages({
  required_packages <- c("ggplot2", "dplyr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, repos = "https://cran.r-project.org/")
  }
  
  library(ggplot2)
  library(dplyr)
})

#' Map mitochondrial gene names to OXPHOS complex numbers
#' 
#' @param gene_name Character vector of gene names
#' @return Character vector of complex assignments (I-V or NA)
map_mt_gene_to_complex <- function(gene_name) {
  gene_upper <- toupper(gene_name)
  
  complex_mapping <- case_when(
    grepl("^ND[1-6]$|^ND4L$", gene_upper) ~ "I",
    gene_upper == "CYTB" ~ "III",
    gene_upper %in% c("COX1", "COX2", "COX3") ~ "IV",
    gene_upper %in% c("ATP6", "ATP8") ~ "V",
    TRUE ~ NA_character_
  )
  
  return(complex_mapping)
}

#' Validate and clean the input data
#' 
#' @param data Data frame to validate
#' @param required_cols Required column names
#' @return Validated data frame
validate_data <- function(data, required_cols) {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Remove rows with missing omega values
  initial_rows <- nrow(data)
  data <- data[!is.na(data$omega_tree), ]
  removed_rows <- initial_rows - nrow(data)
  
  if (removed_rows > 0) {
    cat("Removed", removed_rows, "rows with missing omega values\n")
  }
  
  return(data)
}

# --- Data Loading and Processing ---
cat("Loading data files...\n")

# Read MLC data
tryCatch({
  mlc <- read.table(fn_mlc, sep = "\t", header = TRUE, quote = "", 
                    comment.char = "", stringsAsFactors = FALSE)
  cat("Loaded MLC data:", nrow(mlc), "rows\n")
}, error = function(e) {
  stop("Error reading MLC file: ", e$message)
})

# Read BED data
tryCatch({
  bed <- read.table(fn_bed, sep = "\t", header = FALSE, 
                    stringsAsFactors = FALSE, quote = "")
  colnames(bed) <- c("chr", "start", "end", "gene", "score", "strand", 
                     "complex", "core", "contact")
  cat("Loaded BED data:", nrow(bed), "rows\n")
}, error = function(e) {
  stop("Error reading BED file: ", e$message)
})

# --- Data Filtering and Processing ---
cat("Processing data...\n")

# Filter for gene-level M0 omega values in mt and nu domains
d0 <- mlc %>%
  filter(level == "gene", 
         model == "M0", 
         domain %in% c("nu", "mt")) %>%
  rename(gene = tag) %>%  # Standardize column name
  validate_data(c("gene", "domain", "omega_tree"))

cat("Filtered data:", nrow(d0), "rows (", sum(d0$domain == "mt"), "mt,", 
    sum(d0$domain == "nu"), "nu)\n")

# --- Process Nuclear Genes ---
nu_df <- d0 %>%
  filter(domain == "nu") %>%
  left_join(bed %>% 
            select(gene, complex) %>% 
            distinct() %>%
            mutate(complex = toupper(complex)), 
            by = "gene")

# Report missing complex annotations
missing_nu_genes <- nu_df %>%
  filter(is.na(complex)) %>%
  pull(gene) %>%
  unique()

if (length(missing_nu_genes) > 0) {
  cat("Warning: Nuclear genes missing complex annotation:\n")
  cat("  ", paste(missing_nu_genes, collapse = ", "), "\n")
}

# --- Process Mitochondrial Genes ---
mt_df <- d0 %>%
  filter(domain == "mt") %>%
  mutate(complex = map_mt_gene_to_complex(gene))

# Report MT gene complex assignments
cat("Mitochondrial gene complex assignments:\n")
mt_assignments <- mt_df %>%
  filter(!is.na(complex)) %>%
  select(gene, complex) %>%
  arrange(complex, gene)
print(mt_assignments)

# --- Combine and Prepare Final Dataset ---
df <- bind_rows(
  nu_df %>% select(domain, gene, complex, omega_tree),
  mt_df %>% select(domain, gene, complex, omega_tree)
) %>%
  rename(omega = omega_tree) %>%
  filter(!is.na(complex)) %>%  # Remove genes without complex assignment
  mutate(
    complex = factor(complex, levels = c("I", "II", "III", "IV", "V")),
    domain = factor(domain, levels = c("mt", "nu"))
  )

# Data summary
cat("\nFinal dataset summary:\n")
summary_table <- df %>%
  group_by(complex, domain) %>%
  summarise(
    count = n(),
    mean_omega = round(mean(omega, na.rm = TRUE), 4),
    median_omega = round(median(omega, na.rm = TRUE), 4),
    .groups = "drop"
  )
print(summary_table)

# Cross-tabulation
cat("\nCounts by complex x domain:\n")
print(table(df$complex, df$domain, useNA = "ifany"))

# --- Statistical Analysis ---
# Perform basic statistical tests if there are enough data points
if (nrow(df) > 10) {
  cat("\nPerforming statistical analysis...\n")
  
  # Test for overall differences between domains
  if (length(unique(df$domain)) > 1) {
    domain_test <- wilcox.test(omega ~ domain, data = df)
    cat("Wilcoxon test (mt vs nu overall): p =", format.pval(domain_test$p.value), "\n")
  }
  
  # Test for differences within each complex
  complex_tests <- df %>%
    group_by(complex) %>%
    filter(n_distinct(domain) > 1, n() >= 6) %>%  # Need both domains and enough samples
    do(test = tryCatch({
      wilcox.test(omega ~ domain, data = .)
    }, error = function(e) NULL)) %>%
    filter(!sapply(test, is.null)) %>%
    mutate(p_value = sapply(test, function(x) x$p.value))
  
  if (nrow(complex_tests) > 0) {
    cat("Within-complex comparisons (mt vs nu):\n")
    for (i in 1:nrow(complex_tests)) {
      cat("  Complex", complex_tests$complex[i], ": p =", 
          format.pval(complex_tests$p_value[i]), "\n")
    }
  }
}

# --- Enhanced Plotting ---
cat("Creating enhanced plot...\n")

# Calculate summary statistics for potential annotations
plot_stats <- df %>%
  group_by(complex, domain) %>%
  summarise(
    n = n(),
    median_omega = median(omega, na.rm = TRUE),
    q75 = quantile(omega, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Create the plot with enhanced aesthetics
pd <- position_dodge(width = 0.75)
pj <- position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75, seed = 42)

p <- ggplot(df, aes(x = complex, y = omega, fill = domain)) +
  # Add boxplots
  geom_boxplot(
    position = pd, 
    width = 0.6, 
    outlier.shape = NA, 
    na.rm = TRUE,
    alpha = 0.8,
    color = "black",
    linewidth = 0.5
  ) +
  # Add individual data points
  geom_point(
    aes(group = domain),
    position = pj, 
    size = 2, 
    alpha = 0.7, 
    na.rm = TRUE,
    stroke = 0.3,
    color = "white"
  ) +
  # Enhanced color scheme
  scale_fill_manual(
    values = c("mt" = "#E31A1C", "nu" = "#1F78B4"), 
    labels = c("mt" = "Mitochondrial", "nu" = "Nuclear"),
    drop = FALSE
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1)),
    labels = function(x) sprintf("%.2f", x)
  ) +
  # Labels and theme
  labs(
    x = "OXPHOS Complex", 
    y = expression(omega~"(dN/dS)"), 
    fill = "Genome",
    title = "Evolutionary Rates Across OXPHOS Complexes",
    subtitle = paste("Comparison of dN/dS ratios between mitochondrial and nuclear genes")
  ) +
  theme_classic(base_size = 12) +
  theme(
    # Axis styling
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    
    # Title styling
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    
    # Legend styling
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.margin = margin(5, 5, 5, 5),
    
    # Panel styling
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.2),
    
    # Plot margins
    plot.margin = margin(15, 15, 15, 15)
  )

# Add sample size annotations
p <- p + 
  geom_text(
    data = plot_stats,
    aes(x = complex, y = -Inf, label = paste0("n=", n)),
    position = position_dodge(width = 0.75),
    vjust = -0.5,
    size = 2.8,
    color = "gray40",
    inherit.aes = FALSE
  )

# Save the plot with high quality
tryCatch({
  ggsave(fn_out, p, width = 10, height = 7, dpi = 300, bg = "white")
  cat("Successfully saved plot to:", fn_out, "\n")
}, error = function(e) {
  stop("Error saving plot: ", e$message)
})

# --- Generate Summary Report ---
cat("\n" , rep("=", 60), "\n")
cat("ANALYSIS SUMMARY\n")
cat(rep("=", 60), "\n")

cat("Input files:\n")
cat("  MLC data:", fn_mlc, "\n")
cat("  BED file:", fn_bed, "\n")
cat("Output file:", fn_out, "\n\n")

cat("Data overview:\n")
cat("  Total genes analyzed:", nrow(df), "\n")
cat("  Mitochondrial genes:", sum(df$domain == "mt"), "\n")
cat("  Nuclear genes:", sum(df$domain == "nu"), "\n")
cat("  Complexes represented:", paste(levels(df$complex)[levels(df$complex) %in% df$complex], collapse = ", "), "\n\n")

# Summary statistics by domain
domain_summary <- df %>%
  group_by(domain) %>%
  summarise(
    count = n(),
    mean_omega = round(mean(omega, na.rm = TRUE), 4),
    median_omega = round(median(omega, na.rm = TRUE), 4),
    sd_omega = round(sd(omega, na.rm = TRUE), 4),
    min_omega = round(min(omega, na.rm = TRUE), 4),
    max_omega = round(max(omega, na.rm = TRUE), 4),
    .groups = "drop"
  )

cat("Summary by domain:\n")
print(domain_summary)

# Summary statistics by complex
complex_summary <- df %>%
  group_by(complex, domain) %>%
  summarise(
    count = n(),
    mean_omega = round(mean(omega, na.rm = TRUE), 4),
    median_omega = round(median(omega, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(complex, domain)

cat("\nSummary by complex and domain:\n")
print(complex_summary)

# Report any issues
if (length(missing_nu_genes) > 0) {
  cat("\nWarning: The following nuclear genes were excluded due to missing complex annotation:\n")
  cat("  ", paste(missing_nu_genes, collapse = ", "), "\n")
}

cat("\nAnalysis completed successfully!\n")
cat(rep("=", 60), "\n")