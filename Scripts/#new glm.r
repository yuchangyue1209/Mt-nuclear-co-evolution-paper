#new glm
cat > /work/cyu/ldx_all_subunits/run_ldpruned_mtlineage_quasibinomial_GLM_OXPHOS72.R <<'RSCRIPT'

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ============================================================
# Input files
# ============================================================
infile <- "/work/cyu/ldx_all_subunits/ld/deltaAF_long.noAMO_noLB.ldPruned_kept.keepMarine.R2_0.2.tsv.gz"
lineage_file <- "/work/cyu/mt_lineage_for_glm.tsv"
gene_list_file <- "/mnt/spareHD_2/nu_287/q2_parallelism/oxphos72_genes.list"

outdir <- "/work/cyu/ldx_all_subunits/ld/mtlineage_quasibinomial_GLM_OXPHOS72"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parameters
# ============================================================
min_depth <- 20
min_pops_total <- 6
min_lineages <- 2
min_pops_per_lineage <- 2

effect_threshold <- 0.20
q_threshold <- 0.05

# ============================================================
# Read OXPHOS 72 gene list
# ============================================================
ox72 <- fread(gene_list_file, header = FALSE)
setnames(ox72, "V1", "gene")
ox72[, gene := tolower(trimws(gene))]
ox72 <- unique(ox72)

cat("[INFO] OXPHOS genes loaded:", nrow(ox72), "\n")
print(head(ox72))

# ============================================================
# Read LD-pruned deltaAF table
# ============================================================
cat("[INFO] Reading LD-pruned deltaAF table...\n")
dat <- fread(infile)

cat("[INFO] Input rows before OXPHOS filter:", nrow(dat), "\n")
cat("[INFO] Columns:\n")
print(colnames(dat))

# Standardize names
if ("Pop" %in% names(dat)) setnames(dat, "Pop", "pop")
if ("Gene" %in% names(dat)) setnames(dat, "Gene", "gene")
if ("CHR" %in% names(dat)) setnames(dat, "CHR", "chr")
if ("POS" %in% names(dat)) setnames(dat, "POS", "pos")

dat[, pop := toupper(as.character(pop))]
dat[, gene := tolower(as.character(gene))]
dat[, chr := as.character(chr)]
dat[, pos := as.integer(pos)]

required_cols <- c("chr", "pos", "gene", "pop", "af", "depth")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# SNP ID
if (!"snp" %in% names(dat)) {
  dat[, snp := paste(chr, pos, sep = ":")]
}

dat[, snp_id := paste(gene, chr, pos, sep = "__")]

# ============================================================
# Filter to OXPHOS 72 genes only
# ============================================================
dat <- dat[gene %in% ox72$gene]

cat("[INFO] Rows after OXPHOS72 filter:", nrow(dat), "\n")
cat("[INFO] OXPHOS genes present in data:", uniqueN(dat$gene), "\n")
cat("[INFO] Missing OXPHOS genes from LD-pruned table:\n")
print(setdiff(ox72$gene, unique(dat$gene)))

if (nrow(dat) == 0) {
  stop("No rows left after OXPHOS72 filter. Check gene names.")
}

# Save filtered input
filtered_input <- file.path(outdir, "ldpruned_OXPHOS72_input_before_GLM.tsv.gz")
fwrite(dat, filtered_input, sep = "\t")
cat("[OK] wrote:", filtered_input, "\n")

# ============================================================
# Read mitochondrial lineage file
# ============================================================
lineage <- fread(lineage_file)
setnames(lineage, names(lineage), tolower(names(lineage)))

if (!all(c("pop", "mt_lineage") %in% names(lineage))) {
  stop("lineage_file must contain columns: pop, mt_lineage")
}

lineage[, pop := toupper(as.character(pop))]
lineage[, mt_lineage := as.factor(mt_lineage)]

dat <- merge(dat, lineage, by = "pop", all.x = TRUE)

if (any(is.na(dat$mt_lineage))) {
  cat("[WARNING] Some populations have no mt_lineage:\n")
  print(unique(dat[is.na(mt_lineage), pop]))
  dat <- dat[!is.na(mt_lineage)]
}

# ============================================================
# Construct allele counts
# ============================================================
dat[, af := as.numeric(af)]
dat[, depth := as.numeric(depth)]

dat <- dat[!is.na(af) & !is.na(depth)]
dat <- dat[depth >= min_depth]
dat <- dat[af >= 0 & af <= 1]

# focal count approximated from AF and depth
dat[, focal_count := round(af * depth)]
dat[, focal_count := pmax(0, pmin(focal_count, depth))]
dat[, other_count := depth - focal_count]

dat <- dat[!is.na(focal_count) & !is.na(other_count)]
dat <- dat[depth > 0]

cat("[INFO] Rows after depth/count filtering:", nrow(dat), "\n")
cat("[INFO] Unique OXPHOS72 SNPs:", uniqueN(dat$snp_id), "\n")
cat("[INFO] Unique OXPHOS72 genes:", uniqueN(dat$gene), "\n")
cat("[INFO] Populations:\n")
print(sort(unique(dat$pop)))
cat("[INFO] mt_lineages:\n")
print(table(dat$mt_lineage))

model_input_out <- file.path(outdir, "ldpruned_OXPHOS72_GLM_model_input_with_counts.tsv.gz")
fwrite(dat, model_input_out, sep = "\t")
cat("[OK] wrote:", model_input_out, "\n")

# ============================================================
# Function: test one SNP
# ============================================================
test_one_snp <- function(d) {

  d <- as.data.table(d)

  n_pops <- uniqueN(d$pop)
  n_lineages <- uniqueN(d$mt_lineage)

  if (n_pops < min_pops_total) return(NULL)
  if (n_lineages < min_lineages) return(NULL)

  lineage_counts <- d[, .N, by = mt_lineage]
  if (any(lineage_counts$N < min_pops_per_lineage)) return(NULL)

  if (length(unique(d$focal_count)) < 2 && length(unique(d$other_count)) < 2) {
    return(NULL)
  }

  fit0 <- try(
    glm(
      cbind(focal_count, other_count) ~ 1,
      family = quasibinomial,
      data = d
    ),
    silent = TRUE
  )

  fit1 <- try(
    glm(
      cbind(focal_count, other_count) ~ mt_lineage,
      family = quasibinomial,
      data = d
    ),
    silent = TRUE
  )

  if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) return(NULL)

  a <- try(anova(fit0, fit1, test = "F"), silent = TRUE)
  if (inherits(a, "try-error")) return(NULL)

  pval <- suppressWarnings(a$`Pr(>F)`[2])
  if (is.na(pval)) return(NULL)

  # Effect size: depth-weighted AF by mitochondrial lineage
  lineage_af <- d[, .(
    lineage_focal_count = sum(focal_count, na.rm = TRUE),
    lineage_depth = sum(depth, na.rm = TRUE),
    n_pops_lineage = uniqueN(pop)
  ), by = mt_lineage]

  lineage_af[, weighted_af := lineage_focal_count / lineage_depth]

  max_af <- max(lineage_af$weighted_af, na.rm = TRUE)
  min_af <- min(lineage_af$weighted_af, na.rm = TRUE)

  max_lineage_AF_diff <- max_af - min_af
  max_lineage <- as.character(lineage_af[which.max(weighted_af), mt_lineage])
  min_lineage <- as.character(lineage_af[which.min(weighted_af), mt_lineage])

  pop_af <- d[, .(
    pop_af = sum(focal_count) / sum(depth)
  ), by = .(pop, mt_lineage)]

  max_pop_AF_diff <- max(pop_af$pop_af, na.rm = TRUE) - min(pop_af$pop_af, na.rm = TRUE)

  data.table(
    gene = unique(d$gene)[1],
    chr = unique(d$chr)[1],
    pos = unique(d$pos)[1],
    snp_id = unique(d$snp_id)[1],
    n_pops = n_pops,
    n_lineages = n_lineages,
    p_value = pval,
    dispersion = summary(fit1)$dispersion,
    max_lineage_AF_diff = max_lineage_AF_diff,
    max_pop_AF_diff = max_pop_AF_diff,
    max_lineage = max_lineage,
    min_lineage = min_lineage,
    max_lineage_AF = max_af,
    min_lineage_AF = min_af
  )
}

# ============================================================
# Run GLM per SNP
# ============================================================
cat("[INFO] Running quasibinomial GLM per OXPHOS72 SNP...\n")

snp_list <- split(dat, dat$snp_id)

res_list <- vector("list", length(snp_list))
names(res_list) <- names(snp_list)

counter <- 0

for (id in names(snp_list)) {
  counter <- counter + 1

  if (counter %% 1000 == 0) {
    cat("[INFO] Tested", counter, "SNPs...\n")
  }

  res_list[[id]] <- test_one_snp(snp_list[[id]])
}

res <- rbindlist(res_list, fill = TRUE)

cat("[INFO] OXPHOS72 SNPs with valid models:", nrow(res), "\n")

if (nrow(res) == 0) {
  stop("No SNPs passed filters or model fitting.")
}

# Multiple testing correction
res[, q_value := p.adjust(p_value, method = "BH")]
res[, neglog10_q := -log10(q_value)]
res[is.infinite(neglog10_q), neglog10_q := max(neglog10_q[is.finite(neglog10_q)], na.rm = TRUE) + 1]

res[, significant := q_value < q_threshold]
res[, large_effect := max_lineage_AF_diff >= effect_threshold]
res[, candidate := significant & large_effect]

# Save full result
glm_out <- file.path(outdir, "ldpruned_OXPHOS72_mtlineage_quasibinomial_GLM_results.tsv")
fwrite(res[order(q_value, -max_lineage_AF_diff)], glm_out, sep = "\t")
cat("[OK] wrote:", glm_out, "\n")

# Save candidates
cand <- res[candidate == TRUE][order(q_value, -max_lineage_AF_diff)]

cand_out <- file.path(outdir, paste0(
  "candidate_OXPHOS72_SNPs_q", q_threshold,
  "_maxLineageAFdiff", effect_threshold,
  ".tsv"
))
fwrite(cand, cand_out, sep = "\t")
cat("[OK] wrote:", cand_out, "\n")

# ============================================================
# Gene-level summary
# ============================================================
gene_sum <- res[, .(
  n_tested_snps = .N,
  n_significant_snps = sum(significant, na.rm = TRUE),
  n_large_effect_snps = sum(large_effect, na.rm = TRUE),
  n_candidate_snps = sum(candidate, na.rm = TRUE),
  min_q_value = min(q_value, na.rm = TRUE),
  max_lineage_AF_diff = max(max_lineage_AF_diff, na.rm = TRUE),
  top_snp = snp_id[which.min(q_value)][1],
  top_snp_pos = pos[which.min(q_value)][1],
  top_snp_q = min(q_value, na.rm = TRUE)
), by = gene]

gene_sum <- gene_sum[order(-n_candidate_snps, min_q_value)]

gene_out <- file.path(outdir, "gene_level_summary_OXPHOS72_mtlineage_GLM.tsv")
fwrite(gene_sum, gene_out, sep = "\t")
cat("[OK] wrote:", gene_out, "\n")

# ============================================================
# Print summary
# ============================================================
cat("\n================ OXPHOS72 GLM summary ================\n")
cat("Input rows after OXPHOS72/depth filtering:", nrow(dat), "\n")
cat("Unique OXPHOS72 LD-pruned SNPs:", uniqueN(dat$snp_id), "\n")
cat("OXPHOS72 SNPs tested:", nrow(res), "\n")
cat("Significant SNPs q <", q_threshold, ":", sum(res$significant), "\n")
cat("Large-effect SNPs max lineage AF diff >=", effect_threshold, ":", sum(res$large_effect), "\n")
cat("Candidate SNPs q <", q_threshold, "and effect >=", effect_threshold, ":", sum(res$candidate), "\n")

cat("\nTop OXPHOS72 candidate SNPs:\n")
print(head(cand, 30))

cat("\nTop OXPHOS72 genes:\n")
print(head(gene_sum, 30))

# ============================================================
# Plot 1: effect-size volcano
# ============================================================
res[, plot_group := "not candidate"]
res[significant == TRUE & large_effect == FALSE, plot_group := "FDR only"]
res[significant == FALSE & large_effect == TRUE, plot_group := "large effect only"]
res[candidate == TRUE, plot_group := "candidate"]

top_label <- res[candidate == TRUE][order(q_value, -max_lineage_AF_diff)][1:min(.N, 25)]

res[, label_plot := ""]
res[snp_id %in% top_label$snp_id, label_plot := gene]

p1 <- ggplot(res, aes(x = max_lineage_AF_diff, y = neglog10_q)) +
  geom_point(aes(color = plot_group), alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = effect_threshold, linetype = "dashed") +
  geom_hline(yintercept = -log10(q_threshold), linetype = "dashed") +
  geom_text_repel(
    aes(label = label_plot),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.4,
    min.segment.length = 0
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Maximum allele-frequency difference among mitochondrial lineages",
    y = expression(-log[10]("FDR q-value")),
    color = NULL,
    title = "OXPHOS72 mitochondrial lineage-associated nuclear allele-frequency shifts"
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave(
  file.path(outdir, "Fig4A_OXPHOS72_mtlineage_GLM_effectsize_volcano.png"),
  p1,
  width = 7.2,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "Fig4A_OXPHOS72_mtlineage_GLM_effectsize_volcano.pdf"),
  p1,
  width = 7.2,
  height = 5.5
)

# ============================================================
# Plot 2: candidate genes
# ============================================================
gene_plot <- gene_sum[n_candidate_snps > 0]

if (nrow(gene_plot) > 0) {

  gene_plot[, gene := factor(gene, levels = rev(gene))]

  p2 <- ggplot(gene_plot, aes(x = gene, y = n_candidate_snps)) +
    geom_col(width = 0.75) +
    coord_flip() +
    theme_classic(base_size = 14) +
    labs(
      x = "Gene",
      y = "Number of candidate SNPs",
      title = "Candidate OXPHOS72 genes associated with mitochondrial lineage"
    ) +
    theme(
      axis.text.y = element_text(face = "italic"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )

  ggsave(
    file.path(outdir, "candidate_gene_counts_OXPHOS72_mtlineage_GLM.png"),
    p2,
    width = 6.5,
    height = max(4, 0.25 * nrow(gene_plot)),
    dpi = 300
  )

  ggsave(
    file.path(outdir, "candidate_gene_counts_OXPHOS72_mtlineage_GLM.pdf"),
    p2,
    width = 6.5,
    height = max(4, 0.25 * nrow(gene_plot))
  )
}

# ============================================================
# Plot 3: top candidate AF by lineage
# ============================================================
if (nrow(cand) > 0) {

  top_snps <- cand[1:min(.N, 12), snp_id]

  plot_dat <- dat[snp_id %in% top_snps]

  plot_sum <- plot_dat[, .(
    lineage_focal_count = sum(focal_count),
    lineage_depth = sum(depth),
    weighted_af = sum(focal_count) / sum(depth),
    n_pops = uniqueN(pop)
  ), by = .(gene, chr, pos, snp_id, mt_lineage)]

  plot_sum[, snp_label := paste0(gene, "\n", chr, ":", pos)]

  p3 <- ggplot(plot_sum, aes(x = mt_lineage, y = weighted_af, group = snp_label)) +
    geom_point(size = 2.2) +
    geom_line(alpha = 0.5) +
    facet_wrap(~ snp_label, scales = "free_y") +
    theme_classic(base_size = 13) +
    labs(
      x = "Mitochondrial lineage",
      y = "Depth-weighted focal allele frequency",
      title = "Top OXPHOS72 candidate SNPs by mitochondrial lineage"
    ) +
    theme(
      strip.text = element_text(face = "italic", size = 10),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    file.path(outdir, "top_candidate_OXPHOS72_SNPs_AF_by_mtlineage.png"),
    p3,
    width = 9,
    height = 7,
    dpi = 300
  )

  ggsave(
    file.path(outdir, "top_candidate_OXPHOS72_SNPs_AF_by_mtlineage.pdf"),
    p3,
    width = 9,
    height = 7
  )
}

cat("\nDONE\n")

RSCRIPT

Rscript /work/cyu/ldx_all_subunits/run_ldpruned_mtlineage_quasibinomial_GLM_OXPHOS72.R





cat > /work/cyu/ldx_all_subunits/run_OXPHOS72_separate_geo_mtlineage_GLM.R <<'RSCRIPT'

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
})

# ============================================================
# Input
# ============================================================

INFILE <- "/work/cyu/ldx_all_subunits/ld/mtlineage_quasibinomial_GLM_OXPHOS72/ldpruned_OXPHOS72_GLM_model_input_with_counts.tsv.gz"

OUTDIR <- "/work/cyu/ldx_all_subunits/ld/mtlineage_quasibinomial_GLM_OXPHOS72/separate_geo_mtlineage_GLM"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parameters
# ============================================================

Q_THRESHOLD <- 0.05
GEO_EFFECT_THRESHOLD <- 0.20
MT_EFFECT_THRESHOLD <- 0.20

MIN_POPS_TOTAL <- 6
MIN_GROUPS <- 2
MIN_POPS_PER_GROUP <- 2

# ============================================================
# Style
# ============================================================

COL_CAND <- c(
  "Not candidate" = "grey70",
  "FDR only" = "#4C78A8",
  "Large effect only" = "#E45756",
  "Candidate" = "#C77CFF"
)

COL_OVERLAP <- c(
  "Geography only" = "#00BFC4",
  "mt-lineage only" = "#F8766D",
  "Both" = "#C77CFF",
  "Neither" = "grey80"
)

theme_fig <- theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm")
  )

# ============================================================
# Helpers
# ============================================================

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl = TRUE)
  x
}

safe_anova_p <- function(fit0, fit1){
  a <- try(anova(fit0, fit1, test = "F"), silent = TRUE)
  if (inherits(a, "try-error") || nrow(a) < 2) return(NA_real_)
  as.numeric(a$`Pr(>F)`[2])
}

safe_anova_F <- function(fit0, fit1){
  a <- try(anova(fit0, fit1, test = "F"), silent = TRUE)
  if (inherits(a, "try-error") || nrow(a) < 2) return(NA_real_)
  as.numeric(a$F[2])
}

# ============================================================
# Read data
# ============================================================

cat("[read]", INFILE, "\n")
dat <- fread(INFILE)

dat[, pop := normalize_pop(pop)]
dat[, gene := tolower(as.character(gene))]
dat[, mt_lineage := as.factor(mt_lineage)]

# If region is missing, infer from population names
if (!"region" %in% names(dat)) {
  AK_POPS <- c("FG", "LG", "SL", "SR", "TL", "WB", "WK", "WT")
  BC_POPS <- c("BEA", "BOOT", "ECHO", "GOS", "JOE", "LAW", "MUC", "PYE", "ROB", "SWA", "THE")

  dat[pop %in% AK_POPS, region := "AK"]
  dat[pop %in% BC_POPS, region := "BC"]
}

dat <- dat[region %in% c("AK", "BC")]
dat[, region := factor(region, levels = c("AK", "BC"))]

needed <- c(
  "gene", "chr", "pos", "snp_id", "pop",
  "region", "mt_lineage",
  "focal_count", "other_count", "depth"
)

miss <- setdiff(needed, names(dat))
if (length(miss) > 0) {
  stop("Missing columns: ", paste(miss, collapse = ", "))
}

dat <- dat[
  is.finite(focal_count) &
    is.finite(other_count) &
    is.finite(depth) &
    depth > 0 &
    !is.na(region) &
    !is.na(mt_lineage)
]

cat("[info] rows:", nrow(dat), "\n")
cat("[info] SNPs:", uniqueN(dat$snp_id), "\n")
cat("[info] genes:", uniqueN(dat$gene), "\n")
cat("[info] populations:", uniqueN(dat$pop), "\n")

cat("\n[info] population table:\n")
print(unique(dat[, .(pop, region, mt_lineage)])[order(region, mt_lineage, pop)])

fwrite(
  unique(dat[, .(pop, region, mt_lineage)])[order(region, mt_lineage, pop)],
  file.path(OUTDIR, "pop_region_mtlineage_used.tsv"),
  sep = "\t"
)

# ============================================================
# Generic one-SNP GLM
# ============================================================

test_one_snp <- function(d, group_col){

  d <- as.data.table(d)

  d <- d[
    !is.na(get(group_col)) &
      is.finite(focal_count) &
      is.finite(other_count)
  ]

  n_pops <- uniqueN(d$pop)
  n_groups <- uniqueN(d[[group_col]])

  out_na <- data.table(
    gene = unique(d$gene)[1],
    chr = unique(d$chr)[1],
    pos = unique(d$pos)[1],
    snp_id = unique(d$snp_id)[1],
    n_pops = n_pops,
    n_groups = n_groups,
    test = group_col,
    F_value = NA_real_,
    p_value = NA_real_,
    dispersion = NA_real_,
    max_group_AF_diff = NA_real_,
    max_group = NA_character_,
    min_group = NA_character_,
    max_group_AF = NA_real_,
    min_group_AF = NA_real_
  )

  if (n_pops < MIN_POPS_TOTAL) return(out_na)
  if (n_groups < MIN_GROUPS) return(out_na)

  group_counts <- d[, .(n_pops_group = uniqueN(pop)), by = group_col]
  if (any(group_counts$n_pops_group < MIN_POPS_PER_GROUP)) return(out_na)

  d[, group_tmp := factor(get(group_col))]

  fit0 <- try(
    glm(
      cbind(focal_count, other_count) ~ 1,
      family = quasibinomial,
      data = d
    ),
    silent = TRUE
  )

  fit1 <- try(
    glm(
      cbind(focal_count, other_count) ~ group_tmp,
      family = quasibinomial,
      data = d
    ),
    silent = TRUE
  )

  if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) {
    return(out_na)
  }

  pval <- safe_anova_p(fit0, fit1)
  Fval <- safe_anova_F(fit0, fit1)

  if (!is.finite(pval)) return(out_na)

  group_af <- d[, .(
    focal_sum = sum(focal_count, na.rm = TRUE),
    depth_sum = sum(depth, na.rm = TRUE),
    n_pops_group = uniqueN(pop)
  ), by = group_col]

  group_af <- group_af[depth_sum > 0]
  group_af[, weighted_AF := focal_sum / depth_sum]

  max_af <- max(group_af$weighted_AF, na.rm = TRUE)
  min_af <- min(group_af$weighted_AF, na.rm = TRUE)

  max_group <- as.character(group_af[which.max(weighted_AF)][[group_col]])
  min_group <- as.character(group_af[which.min(weighted_AF)][[group_col]])

  out <- data.table(
    gene = unique(d$gene)[1],
    chr = unique(d$chr)[1],
    pos = unique(d$pos)[1],
    snp_id = unique(d$snp_id)[1],
    n_pops = n_pops,
    n_groups = n_groups,
    test = group_col,
    F_value = Fval,
    p_value = pval,
    dispersion = summary(fit1)$dispersion,
    max_group_AF_diff = max_af - min_af,
    max_group = max_group,
    min_group = min_group,
    max_group_AF = max_af,
    min_group_AF = min_af
  )

  # explicit AK and BC AF for geography
  if (group_col == "region") {
    af_region <- d[, .(
      weighted_AF = sum(focal_count) / sum(depth)
    ), by = region]

    ak_af <- af_region[region == "AK", weighted_AF]
    bc_af <- af_region[region == "BC", weighted_AF]

    out[, AK_AF := ifelse(length(ak_af) == 1, ak_af, NA_real_)]
    out[, BC_AF := ifelse(length(bc_af) == 1, bc_af, NA_real_)]
    out[, geo_AF_diff_BC_minus_AK := BC_AF - AK_AF]
    out[, abs_geo_AF_diff := abs(geo_AF_diff_BC_minus_AK)]
  }

  out
}

run_glm_by_group <- function(dat, group_col){

  cat("\n[run] GLM:", group_col, "\n")

  snps <- unique(dat$snp_id)
  setkey(dat, snp_id)

  res_list <- vector("list", length(snps))

  for (i in seq_along(snps)) {

    if (i %% 1000 == 0) {
      cat("[progress]", group_col, i, "/", length(snps), "\n")
    }

    s <- snps[i]
    d <- dat[list(s)]

    res_list[[i]] <- test_one_snp(d, group_col)
  }

  res <- rbindlist(res_list, fill = TRUE)

  res <- res[is.finite(p_value) & p_value > 0 & p_value <= 1]
  res[, q_value := p.adjust(p_value, method = "BH")]
  res[, neglog10_q := -log10(q_value)]

  res[]
}

# ============================================================
# Run separate GLMs
# ============================================================

res_geo <- run_glm_by_group(dat, "region")
res_mt  <- run_glm_by_group(dat, "mt_lineage")

# ============================================================
# Candidate definitions
# ============================================================

res_geo[, significant := q_value < Q_THRESHOLD]
res_geo[, large_effect := abs_geo_AF_diff >= GEO_EFFECT_THRESHOLD]
res_geo[, candidate := significant & large_effect]

res_mt[, significant := q_value < Q_THRESHOLD]
res_mt[, large_effect := max_group_AF_diff >= MT_EFFECT_THRESHOLD]
res_mt[, candidate := significant & large_effect]

# Rename columns
setnames(
  res_geo,
  old = c(
    "F_value", "p_value", "q_value", "neglog10_q",
    "max_group_AF_diff", "max_group", "min_group",
    "max_group_AF", "min_group_AF",
    "significant", "large_effect", "candidate"
  ),
  new = c(
    "F_geo", "p_geo", "q_geo", "neglog10_q_geo",
    "max_geo_AF_diff", "max_geo_group", "min_geo_group",
    "max_geo_AF", "min_geo_AF",
    "significant_geo", "large_effect_geo", "candidate_geo"
  )
)

setnames(
  res_mt,
  old = c(
    "F_value", "p_value", "q_value", "neglog10_q",
    "max_group_AF_diff", "max_group", "min_group",
    "max_group_AF", "min_group_AF",
    "significant", "large_effect", "candidate"
  ),
  new = c(
    "F_mt", "p_mt", "q_mt", "neglog10_q_mt",
    "max_lineage_AF_diff", "max_lineage", "min_lineage",
    "max_lineage_AF", "min_lineage_AF",
    "significant_mt", "large_effect_mt", "candidate_mt"
  )
)

# ============================================================
# Save outputs
# ============================================================

geo_out <- file.path(OUTDIR, "OXPHOS72_geography_quasibinomial_GLM_results.tsv")
mt_out  <- file.path(OUTDIR, "OXPHOS72_mtlineage_quasibinomial_GLM_results.tsv")

fwrite(res_geo[order(q_geo, -abs_geo_AF_diff)], geo_out, sep = "\t")
fwrite(res_mt[order(q_mt, -max_lineage_AF_diff)], mt_out, sep = "\t")

cat("[write]", geo_out, "\n")
cat("[write]", mt_out, "\n")

cand_geo_out <- file.path(OUTDIR, "candidate_OXPHOS72_geography_q0.05_absGeoAFdiff0.2.tsv")
cand_mt_out  <- file.path(OUTDIR, "candidate_OXPHOS72_mtlineage_q0.05_maxLineageAFdiff0.2.tsv")

fwrite(res_geo[candidate_geo == TRUE][order(q_geo, -abs_geo_AF_diff)], cand_geo_out, sep = "\t")
fwrite(res_mt[candidate_mt == TRUE][order(q_mt, -max_lineage_AF_diff)], cand_mt_out, sep = "\t")

# ============================================================
# Merge overlap
# ============================================================

geo_keep <- res_geo[, .(
  snp_id, gene, chr, pos,
  p_geo, q_geo, neglog10_q_geo,
  AK_AF, BC_AF,
  geo_AF_diff_BC_minus_AK,
  abs_geo_AF_diff,
  significant_geo,
  large_effect_geo,
  candidate_geo
)]

mt_keep <- res_mt[, .(
  snp_id,
  p_mt, q_mt, neglog10_q_mt,
  max_lineage_AF_diff,
  max_lineage,
  min_lineage,
  max_lineage_AF,
  min_lineage_AF,
  significant_mt,
  large_effect_mt,
  candidate_mt
)]

cmp <- merge(geo_keep, mt_keep, by = "snp_id", all = TRUE)

cmp[, candidate_geo := fifelse(is.na(candidate_geo), FALSE, candidate_geo)]
cmp[, candidate_mt := fifelse(is.na(candidate_mt), FALSE, candidate_mt)]

cmp[, category := fifelse(
  candidate_geo & candidate_mt,
  "Both",
  fifelse(
    candidate_geo & !candidate_mt,
    "Geography only",
    fifelse(
      !candidate_geo & candidate_mt,
      "mt-lineage only",
      "Neither"
    )
  )
)]

cmp[, category := factor(
  category,
  levels = c("Geography only", "mt-lineage only", "Both", "Neither")
)]

cmp_out <- file.path(OUTDIR, "OXPHOS72_geography_vs_mtlineage_candidate_overlap.tsv")
fwrite(cmp[order(category, q_geo, q_mt)], cmp_out, sep = "\t")
cat("[write]", cmp_out, "\n")

overlap_summary <- cmp[, .N, by = category][order(category)]
fwrite(overlap_summary, file.path(OUTDIR, "candidate_overlap_summary.tsv"), sep = "\t")

# ============================================================
# Gene summaries
# ============================================================

gene_geo <- res_geo[, .(
  n_tested_snps = .N,
  n_significant_geo = sum(significant_geo, na.rm = TRUE),
  n_large_effect_geo = sum(large_effect_geo, na.rm = TRUE),
  n_candidate_geo = sum(candidate_geo, na.rm = TRUE),
  min_q_geo = min(q_geo, na.rm = TRUE),
  max_abs_geo_AF_diff = max(abs_geo_AF_diff, na.rm = TRUE),
  top_geo_snp = snp_id[which.min(q_geo)][1],
  top_geo_pos = pos[which.min(q_geo)][1]
), by = gene][order(-n_candidate_geo, min_q_geo)]

gene_mt <- res_mt[, .(
  n_tested_snps = .N,
  n_significant_mt = sum(significant_mt, na.rm = TRUE),
  n_large_effect_mt = sum(large_effect_mt, na.rm = TRUE),
  n_candidate_mt = sum(candidate_mt, na.rm = TRUE),
  min_q_mt = min(q_mt, na.rm = TRUE),
  max_lineage_AF_diff = max(max_lineage_AF_diff, na.rm = TRUE),
  top_mt_snp = snp_id[which.min(q_mt)][1],
  top_mt_pos = pos[which.min(q_mt)][1]
), by = gene][order(-n_candidate_mt, min_q_mt)]

fwrite(gene_geo, file.path(OUTDIR, "gene_summary_OXPHOS72_geography_GLM.tsv"), sep = "\t")
fwrite(gene_mt, file.path(OUTDIR, "gene_summary_OXPHOS72_mtlineage_GLM.tsv"), sep = "\t")

# ============================================================
# Summary print
# ============================================================

cat("\n================ Separate OXPHOS72 GLM summary ================\n")

cat("\nGeography model:\n")
cat("SNPs tested:", nrow(res_geo), "\n")
cat("Significant q < 0.05:", sum(res_geo$significant_geo, na.rm = TRUE), "\n")
cat("Large effect abs AK-BC AF diff >= 0.2:", sum(res_geo$large_effect_geo, na.rm = TRUE), "\n")
cat("Candidates:", sum(res_geo$candidate_geo, na.rm = TRUE), "\n")

cat("\nmt-lineage model:\n")
cat("SNPs tested:", nrow(res_mt), "\n")
cat("Significant q < 0.05:", sum(res_mt$significant_mt, na.rm = TRUE), "\n")
cat("Large effect max lineage AF diff >= 0.2:", sum(res_mt$large_effect_mt, na.rm = TRUE), "\n")
cat("Candidates:", sum(res_mt$candidate_mt, na.rm = TRUE), "\n")

cat("\nCandidate overlap:\n")
print(overlap_summary)

cat("\nTop geography candidates:\n")
print(head(res_geo[candidate_geo == TRUE][order(q_geo, -abs_geo_AF_diff)], 20))

cat("\nTop mt-lineage candidates:\n")
print(head(res_mt[candidate_mt == TRUE][order(q_mt, -max_lineage_AF_diff)], 20))

# ============================================================
# Plot A1: geography volcano
# ============================================================

res_geo[, plot_group := "Not candidate"]
res_geo[significant_geo == TRUE & large_effect_geo == FALSE, plot_group := "FDR only"]
res_geo[significant_geo == FALSE & large_effect_geo == TRUE, plot_group := "Large effect only"]
res_geo[candidate_geo == TRUE, plot_group := "Candidate"]

res_geo[, plot_group := factor(
  plot_group,
  levels = c("Not candidate", "FDR only", "Large effect only", "Candidate")
)]

TOP_GEO <- res_geo[candidate_geo == TRUE][order(q_geo, -abs_geo_AF_diff)]
TOP_GEO <- TOP_GEO[1:min(.N, 15)]

res_geo[, neglog10_q_geo_plot := pmin(neglog10_q_geo, 35)]
TOP_GEO[, neglog10_q_geo_plot := pmin(neglog10_q_geo, 35)]

p_geo <- ggplot(res_geo, aes(x = abs_geo_AF_diff, y = neglog10_q_geo_plot)) +
  geom_point(aes(color = plot_group), alpha = 0.75, size = 1.7) +
  geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = GEO_EFFECT_THRESHOLD, linetype = "dashed", linewidth = 0.4) +
  geom_text_repel(
    data = TOP_GEO,
    aes(label = gene),
    size = 3.8,
    color = "black",
    fontface = "italic",
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.3,
    max.overlaps = Inf,
    force = 2
  ) +
  scale_color_manual(values = COL_CAND, name = NULL) +
  labs(
    title = "Geography",
    x = "Absolute AK–BC allele-frequency difference",
    y = expression(-log[10]("FDR q-value"))
  ) +
  theme_fig +
  theme(legend.position = "right")

ggsave(
  file.path(OUTDIR, "Fig_geography_OXPHOS72_GLM_volcano.png"),
  p_geo,
  width = 7.2,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig_geography_OXPHOS72_GLM_volcano.pdf"),
  p_geo,
  width = 7.2,
  height = 5.5
)

# ============================================================
# Plot A2: mt-lineage volcano
# ============================================================

res_mt[, plot_group := "Not candidate"]
res_mt[significant_mt == TRUE & large_effect_mt == FALSE, plot_group := "FDR only"]
res_mt[significant_mt == FALSE & large_effect_mt == TRUE, plot_group := "Large effect only"]
res_mt[candidate_mt == TRUE, plot_group := "Candidate"]

res_mt[, plot_group := factor(
  plot_group,
  levels = c("Not candidate", "FDR only", "Large effect only", "Candidate")
)]

TOP_MT <- res_mt[candidate_mt == TRUE][order(q_mt, -max_lineage_AF_diff)]
TOP_MT <- TOP_MT[1:min(.N, 15)]

res_mt[, neglog10_q_mt_plot := pmin(neglog10_q_mt, 35)]
TOP_MT[, neglog10_q_mt_plot := pmin(neglog10_q_mt, 35)]

p_mt <- ggplot(res_mt, aes(x = max_lineage_AF_diff, y = neglog10_q_mt_plot)) +
  geom_point(aes(color = plot_group), alpha = 0.75, size = 1.7) +
  geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = MT_EFFECT_THRESHOLD, linetype = "dashed", linewidth = 0.4) +
  geom_text_repel(
    data = TOP_MT,
    aes(label = gene),
    size = 3.8,
    color = "black",
    fontface = "italic",
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.3,
    max.overlaps = Inf,
    force = 2
  ) +
  scale_color_manual(values = COL_CAND, name = NULL) +
  labs(
    title = "Mitochondrial lineage",
    x = "Maximum allele-frequency difference among mt lineages",
    y = expression(-log[10]("FDR q-value"))
  ) +
  theme_fig +
  theme(legend.position = "right")

ggsave(
  file.path(OUTDIR, "Fig_mtlineage_OXPHOS72_GLM_volcano.png"),
  p_mt,
  width = 7.2,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig_mtlineage_OXPHOS72_GLM_volcano.pdf"),
  p_mt,
  width = 7.2,
  height = 5.5
)

# ============================================================
# Plot: candidate overlap Venn-style
# ============================================================

total_oxphos <- uniqueN(cmp$snp_id)
geo_total <- sum(cmp$candidate_geo, na.rm = TRUE)
mt_total <- sum(cmp$candidate_mt, na.rm = TRUE)
overlap <- sum(cmp$candidate_geo & cmp$candidate_mt, na.rm = TRUE)

geo_only <- sum(cmp$candidate_geo & !cmp$candidate_mt, na.rm = TRUE)
mt_only <- sum(!cmp$candidate_geo & cmp$candidate_mt, na.rm = TRUE)

VENN_LABELS <- data.table(
  label = c(
    paste0(total_oxphos, "\nOXPHOS72 SNPs"),
    paste0(geo_only, "\nGeography only"),
    paste0(mt_only, "\nmt-lineage only"),
    paste0(overlap, "\nBoth")
  ),
  x = c(0.50, 0.67, 0.33, 0.50),
  y = c(0.78, 0.43, 0.43, 0.43)
)

p_venn <- ggplot() +
  geom_ellipse(
    aes(x0 = 0.50, y0 = 0.52, a = 0.43, b = 0.31, angle = 0),
    fill = "grey92",
    alpha = 0.45,
    color = "grey40",
    linewidth = 0.8
  ) +
  geom_ellipse(
    aes(x0 = 0.39, y0 = 0.43, a = 0.24, b = 0.18, angle = 0),
    fill = "#F8766D",
    alpha = 0.35,
    color = "grey35",
    linewidth = 0.8
  ) +
  geom_ellipse(
    aes(x0 = 0.61, y0 = 0.43, a = 0.24, b = 0.18, angle = 0),
    fill = "#00BFC4",
    alpha = 0.35,
    color = "grey35",
    linewidth = 0.8
  ) +
  geom_text(
    data = VENN_LABELS,
    aes(x = x, y = y, label = label),
    size = 4.6,
    fontface = "bold",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 0.22,
    y = 0.20,
    label = "mt-lineage",
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.78,
    y = 0.20,
    label = "Geography",
    size = 5,
    fontface = "bold"
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0.12, 0.90)) +
  theme_void()

ggsave(
  file.path(OUTDIR, "Fig_candidate_overlap_venn_style_OXPHOS72.png"),
  p_venn,
  width = 5.8,
  height = 4.6,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig_candidate_overlap_venn_style_OXPHOS72.pdf"),
  p_venn,
  width = 5.8,
  height = 4.6
)

cat("\nDONE\n")
cat("Output dir:\n", OUTDIR, "\n")

RSCRIPT

Rscript /work/cyu/ldx_all_subunits/run_OXPHOS72_separate_geo_mtlineage_GLM.R