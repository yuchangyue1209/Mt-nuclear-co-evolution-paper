#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(patchwork)
  library(readxl)
})

# ============================================================
# OXPHOS72 geography vs mt-lineage GLM
#
# Updated according to Jesse:
#   1. Geography = latitude + longitude, not AK/BC region.
#   2. First mt-lineage analysis = all retained populations combined.
#   3. Region-specific mt-lineage analyses are follow-up tables only.
#   4. Main figure = panels A and B only.
#   5. Former panel C is saved as a separate figure.
#
# Inputs:
#   LD-pruned deltaAF table from withAMO_noLB + DP20 + R2 0.2
# ============================================================

# ============================================================
# Input / output
# ============================================================

INFILE <- "/work/cyu/ldx_all_subunits_withAMO_noLB_DP20/ld_pruned_R2_0.2/deltaAF_long.withAMO_noLB.ldPruned_kept.keepMarine.DP20.R2_0.2.tsv.gz"

GENE_LIST_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/oxphos72_genes.withAMO_noLB.list"

POOLINFO <- "/work/cyu/Poolinfo.csv"

CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"

OUTDIR <- "/work/cyu/ldx_all_subunits_withAMO_noLB_DP20/ld_pruned_R2_0.2/geoLatLong_mtlineage_GLM_OXPHOS72"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parameters
# ============================================================

Q_THRESHOLD <- 0.05
GEO_EFFECT_THRESHOLD <- 0.20
MT_EFFECT_THRESHOLD <- 0.20

MIN_DEPTH <- 20

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
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?(_DEDUP)?$", "\\2", x, perl = TRUE)
  x
}

safe_read_poolinfo <- function(f){
  cat("[read] Poolinfo:", f, "\n")

  # file may be Excel-formatted even if named .csv
  magic <- readBin(f, what = "raw", n = 4)
  is_xlsx <- identical(as.integer(magic[1:2]), as.integer(charToRaw("PK")))

  if (is_xlsx) {
    as.data.table(read_excel(f))
  } else {
    fread(f)
  }
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
# Read OXPHOS72 gene list
# ============================================================

ox72 <- fread(GENE_LIST_FILE, header = FALSE)
setnames(ox72, "V1", "gene")
ox72[, gene := tolower(trimws(gene))]
ox72 <- unique(ox72)

cat("[INFO] OXPHOS genes loaded:", nrow(ox72), "\n")

# ============================================================
# Read LD-pruned deltaAF table
# ============================================================

cat("[read] LD-pruned deltaAF:", INFILE, "\n")
dat <- fread(INFILE)

if ("Pop" %in% names(dat)) setnames(dat, "Pop", "pop")
if ("Gene" %in% names(dat)) setnames(dat, "Gene", "gene")
if ("CHR" %in% names(dat)) setnames(dat, "CHR", "chr")
if ("POS" %in% names(dat)) setnames(dat, "POS", "pos")

dat[, pop := normalize_pop(pop)]
dat[, gene := tolower(as.character(gene))]
dat[, chr := as.character(chr)]
dat[, pos := as.integer(pos)]

required_cols <- c("chr", "pos", "gene", "pop", "af", "depth")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!"snp" %in% names(dat)) {
  dat[, snp := paste(chr, pos, gene, sep = ":")]
}

dat[, snp_id := paste(gene, chr, pos, sep = "__")]

# OXPHOS72 only
dat <- dat[gene %in% ox72$gene]

# depth/count filtering
dat[, af := as.numeric(af)]
dat[, depth := as.numeric(depth)]

dat <- dat[
  is.finite(af) &
    is.finite(depth) &
    depth >= MIN_DEPTH &
    af >= 0 &
    af <= 1
]

dat[, focal_count := round(af * depth)]
dat[, focal_count := pmax(0, pmin(focal_count, depth))]
dat[, other_count := depth - focal_count]

cat("[INFO] rows after OXPHOS72/depth filtering:", nrow(dat), "\n")
cat("[INFO] SNPs:", uniqueN(dat$snp_id), "\n")
cat("[INFO] genes:", uniqueN(dat$gene), "\n")

# ============================================================
# Read geography metadata
# ============================================================

info <- safe_read_poolinfo(POOLINFO)

setnames(info, old = "Population", new = "pop", skip_absent = TRUE)

required_info <- c("pop", "Latitude", "Longitude")
missing_info <- setdiff(required_info, names(info))
if (length(missing_info) > 0) {
  stop("Poolinfo missing columns: ", paste(missing_info, collapse = ", "))
}

info[, pop := normalize_pop(pop)]
info[, Latitude := as.numeric(Latitude)]
info[, Longitude := as.numeric(Longitude)]

if (!"Region" %in% names(info)) {
  info[, Region := NA_character_]
}
if (!"Habitat" %in% names(info)) {
  info[, Habitat := NA_character_]
}
if (!"Watershed" %in% names(info)) {
  info[, Watershed := NA_character_]
}

info <- unique(
  info[, .(pop, Region, Latitude, Longitude, Habitat, Watershed)],
  by = "pop"
)

# ============================================================
# Read mtCluster metadata
# ============================================================

cat("[read] mtCluster:", CLUSTER_FILE, "\n")
cl <- fread(CLUSTER_FILE, sep = "\t", header = TRUE, fill = TRUE, blank.lines.skip = TRUE)

need_cl <- c("pop", "mtCluster")
missing_cl <- setdiff(need_cl, names(cl))
if (length(missing_cl) > 0) {
  stop("Cluster file missing columns: ", paste(missing_cl, collapse = ", "))
}

cl <- cl[!is.na(pop) & pop != "" & !is.na(mtCluster) & mtCluster != ""]
cl[, pop := normalize_pop(pop)]
cl <- unique(cl[, .(pop, mtCluster)], by = "pop")
cl[, mt_lineage := factor(mtCluster)]

# ============================================================
# Merge metadata
# ============================================================

dat <- merge(dat, info, by = "pop", all.x = TRUE)
dat <- merge(dat, cl[, .(pop, mt_lineage)], by = "pop", all.x = TRUE)

# Keep only rows with geography and mt-lineage.
# AMO will be retained only if it has mt_lineage in mtCluster_manual.
dat <- dat[
  is.finite(Latitude) &
    is.finite(Longitude) &
    !is.na(mt_lineage)
]

# Infer region if missing
AK_POPS <- c("FG", "LG", "SL", "SR", "TL", "WB", "WK", "WT")
BC_POPS <- c("AMO", "BEA", "BOOT", "ECHO", "GOS", "JOE", "LAW", "MUC", "PYE", "ROB", "SWA", "THE")

dat[pop %in% AK_POPS, region := "AK"]
dat[pop %in% BC_POPS, region := "BC"]
dat <- dat[region %in% c("AK", "BC")]
dat[, region := factor(region, levels = c("AK", "BC"))]

# scale geography across all included freshwater populations
dat[, Latitude_z := as.numeric(scale(Latitude))]
dat[, Longitude_z := as.numeric(scale(Longitude))]

cat("\n[INFO] population table used:\n")
pop_used <- unique(dat[, .(pop, region, Latitude, Longitude, mt_lineage)])[order(region, mt_lineage, pop)]
print(pop_used)

fwrite(
  pop_used,
  file.path(OUTDIR, "pop_geo_mtlineage_used.tsv"),
  sep = "\t"
)

model_input_out <- file.path(OUTDIR, "OXPHOS72_model_input_with_counts_geo_mtlineage.tsv.gz")
fwrite(dat, model_input_out, sep = "\t")
cat("[OK] wrote:", model_input_out, "\n")

# ============================================================
# GLM functions
# ============================================================

test_geo_one_snp <- function(d){

  d <- as.data.table(d)

  d <- d[
    is.finite(focal_count) &
      is.finite(other_count) &
      is.finite(depth) &
      is.finite(Latitude_z) &
      is.finite(Longitude_z)
  ]

  n_pops <- uniqueN(d$pop)

  out_na <- data.table(
    gene = unique(d$gene)[1],
    chr = unique(d$chr)[1],
    pos = unique(d$pos)[1],
    snp_id = unique(d$snp_id)[1],
    n_pops = n_pops,
    test = "geography_lat_long",
    F_geo = NA_real_,
    p_geo = NA_real_,
    dispersion_geo = NA_real_,
    geo_fitted_AF_range = NA_real_,
    max_pop_AF_diff = NA_real_
  )

  if (n_pops < MIN_POPS_TOTAL) return(out_na)

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
      cbind(focal_count, other_count) ~ Latitude_z + Longitude_z,
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

  # geography effect size:
  # range of fitted AF across observed population coordinates
  d[, fitted_AF_geo := as.numeric(predict(fit1, type = "response"))]
  geo_range <- max(d$fitted_AF_geo, na.rm = TRUE) - min(d$fitted_AF_geo, na.rm = TRUE)

  pop_af <- d[, .(
    pop_af = sum(focal_count, na.rm = TRUE) / sum(depth, na.rm = TRUE)
  ), by = pop]

  max_pop_diff <- max(pop_af$pop_af, na.rm = TRUE) - min(pop_af$pop_af, na.rm = TRUE)

  data.table(
    gene = unique(d$gene)[1],
    chr = unique(d$chr)[1],
    pos = unique(d$pos)[1],
    snp_id = unique(d$snp_id)[1],
    n_pops = n_pops,
    test = "geography_lat_long",
    F_geo = Fval,
    p_geo = pval,
    dispersion_geo = summary(fit1)$dispersion,
    geo_fitted_AF_range = geo_range,
    max_pop_AF_diff = max_pop_diff
  )
}

test_mt_one_snp <- function(d){

  d <- as.data.table(d)

  d <- d[
    !is.na(mt_lineage) &
      is.finite(focal_count) &
      is.finite(other_count) &
      is.finite(depth)
  ]

  n_pops <- uniqueN(d$pop)
  n_lineages <- uniqueN(d$mt_lineage)

  out_na <- data.table(
    gene = unique(d$gene)[1],
    chr = unique(d$chr)[1],
    pos = unique(d$pos)[1],
    snp_id = unique(d$snp_id)[1],
    n_pops = n_pops,
    n_lineages = n_lineages,
    test = "mt_lineage_all_pops",
    F_mt = NA_real_,
    p_mt = NA_real_,
    dispersion_mt = NA_real_,
    max_lineage_AF_diff = NA_real_,
    max_lineage = NA_character_,
    min_lineage = NA_character_,
    max_lineage_AF = NA_real_,
    min_lineage_AF = NA_real_
  )

  if (n_pops < MIN_POPS_TOTAL) return(out_na)
  if (n_lineages < MIN_GROUPS) return(out_na)

  lineage_counts <- d[, .(n_pops_lineage = uniqueN(pop)), by = mt_lineage]
  if (any(lineage_counts$n_pops_lineage < MIN_POPS_PER_GROUP)) return(out_na)

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

  if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) {
    return(out_na)
  }

  pval <- safe_anova_p(fit0, fit1)
  Fval <- safe_anova_F(fit0, fit1)

  if (!is.finite(pval)) return(out_na)

  lineage_af <- d[, .(
    lineage_focal_count = sum(focal_count, na.rm = TRUE),
    lineage_depth = sum(depth, na.rm = TRUE),
    n_pops_lineage = uniqueN(pop)
  ), by = mt_lineage]

  lineage_af <- lineage_af[lineage_depth > 0]
  lineage_af[, weighted_AF := lineage_focal_count / lineage_depth]

  max_af <- max(lineage_af$weighted_AF, na.rm = TRUE)
  min_af <- min(lineage_af$weighted_AF, na.rm = TRUE)

  max_lineage <- as.character(lineage_af[which.max(weighted_AF), mt_lineage])
  min_lineage <- as.character(lineage_af[which.min(weighted_AF), mt_lineage])

  data.table(
    gene = unique(d$gene)[1],
    chr = unique(d$chr)[1],
    pos = unique(d$pos)[1],
    snp_id = unique(d$snp_id)[1],
    n_pops = n_pops,
    n_lineages = n_lineages,
    test = "mt_lineage_all_pops",
    F_mt = Fval,
    p_mt = pval,
    dispersion_mt = summary(fit1)$dispersion,
    max_lineage_AF_diff = max_af - min_af,
    max_lineage = max_lineage,
    min_lineage = min_lineage,
    max_lineage_AF = max_af,
    min_lineage_AF = min_af
  )
}

run_by_snp <- function(dat, FUN, label){

  cat("\n[run]", label, "\n")

  snps <- unique(dat$snp_id)
  setkey(dat, snp_id)

  out <- vector("list", length(snps))

  for (i in seq_along(snps)) {
    if (i %% 1000 == 0) {
      cat("[progress]", label, i, "/", length(snps), "\n")
    }

    s <- snps[i]
    out[[i]] <- FUN(dat[list(s)])
  }

  rbindlist(out, fill = TRUE)
}

# ============================================================
# Main analyses:
# A. Geography = latitude + longitude
# B. mt-lineage = all populations combined
# ============================================================

res_geo <- run_by_snp(dat, test_geo_one_snp, "geography latitude + longitude")
res_mt  <- run_by_snp(dat, test_mt_one_snp, "mt-lineage all populations combined")

res_geo <- res_geo[is.finite(p_geo) & p_geo > 0 & p_geo <= 1]
res_mt  <- res_mt[is.finite(p_mt) & p_mt > 0 & p_mt <= 1]

res_geo[, q_geo := p.adjust(p_geo, method = "BH")]
res_geo[, neglog10_q_geo := -log10(q_geo)]
res_geo[, significant_geo := q_geo < Q_THRESHOLD]
res_geo[, large_effect_geo := geo_fitted_AF_range >= GEO_EFFECT_THRESHOLD]
res_geo[, candidate_geo := significant_geo & large_effect_geo]

res_mt[, q_mt := p.adjust(p_mt, method = "BH")]
res_mt[, neglog10_q_mt := -log10(q_mt)]
res_mt[, significant_mt := q_mt < Q_THRESHOLD]
res_mt[, large_effect_mt := max_lineage_AF_diff >= MT_EFFECT_THRESHOLD]
res_mt[, candidate_mt := significant_mt & large_effect_mt]

# Save main tables
geo_out <- file.path(OUTDIR, "OXPHOS72_geography_LatLong_quasibinomial_GLM_results.tsv")
mt_out  <- file.path(OUTDIR, "OXPHOS72_mtlineage_allPops_quasibinomial_GLM_results.tsv")

fwrite(res_geo[order(q_geo, -geo_fitted_AF_range)], geo_out, sep = "\t")
fwrite(res_mt[order(q_mt, -max_lineage_AF_diff)], mt_out, sep = "\t")

cat("[OK] wrote:", geo_out, "\n")
cat("[OK] wrote:", mt_out, "\n")

# Candidate tables
fwrite(
  res_geo[candidate_geo == TRUE][order(q_geo, -geo_fitted_AF_range)],
  file.path(OUTDIR, "candidate_OXPHOS72_geography_LatLong_q0.05_effect0.2.tsv"),
  sep = "\t"
)

fwrite(
  res_mt[candidate_mt == TRUE][order(q_mt, -max_lineage_AF_diff)],
  file.path(OUTDIR, "candidate_OXPHOS72_mtlineage_allPops_q0.05_effect0.2.tsv"),
  sep = "\t"
)

# ============================================================
# Merge overlap for separate panel C
# ============================================================

geo_keep <- res_geo[, .(
  snp_id, gene, chr, pos,
  p_geo, q_geo, neglog10_q_geo,
  geo_fitted_AF_range,
  max_pop_AF_diff,
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

cmp_out <- file.path(OUTDIR, "OXPHOS72_geographyLatLong_vs_mtlineage_candidate_overlap.tsv")
fwrite(cmp[order(category, q_geo, q_mt)], cmp_out, sep = "\t")

overlap_summary <- cmp[, .N, by = category][order(category)]
fwrite(overlap_summary, file.path(OUTDIR, "candidate_overlap_summary.tsv"), sep = "\t")

# ============================================================
# Gene-level summaries
# ============================================================

gene_geo <- res_geo[, .(
  n_tested_snps = .N,
  n_significant_geo = sum(significant_geo, na.rm = TRUE),
  n_large_effect_geo = sum(large_effect_geo, na.rm = TRUE),
  n_candidate_geo = sum(candidate_geo, na.rm = TRUE),
  min_q_geo = min(q_geo, na.rm = TRUE),
  max_geo_fitted_AF_range = max(geo_fitted_AF_range, na.rm = TRUE),
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

fwrite(gene_geo, file.path(OUTDIR, "gene_summary_OXPHOS72_geography_LatLong_GLM.tsv"), sep = "\t")
fwrite(gene_mt, file.path(OUTDIR, "gene_summary_OXPHOS72_mtlineage_allPops_GLM.tsv"), sep = "\t")

# ============================================================
# Follow-up: region-specific mt-lineage models
# Not part of main figure
# ============================================================

run_region_mt <- function(REGION_NAME){

  sub <- dat[region == REGION_NAME]

  if (nrow(sub) == 0) {
    return(NULL)
  }

  res <- run_by_snp(sub, test_mt_one_snp, paste0("mt-lineage within ", REGION_NAME))
  res <- res[is.finite(p_mt) & p_mt > 0 & p_mt <= 1]

  if (nrow(res) == 0) {
    return(NULL)
  }

  res[, region_followup := REGION_NAME]
  res[, q_mt := p.adjust(p_mt, method = "BH")]
  res[, neglog10_q_mt := -log10(q_mt)]
  res[, significant_mt := q_mt < Q_THRESHOLD]
  res[, large_effect_mt := max_lineage_AF_diff >= MT_EFFECT_THRESHOLD]
  res[, candidate_mt := significant_mt & large_effect_mt]

  res[]
}

res_mt_AK <- run_region_mt("AK")
res_mt_BC <- run_region_mt("BC")

res_region <- rbindlist(list(res_mt_AK, res_mt_BC), fill = TRUE)

if (nrow(res_region) > 0) {
  fwrite(
    res_region[order(region_followup, q_mt, -max_lineage_AF_diff)],
    file.path(OUTDIR, "followup_regionSpecific_mtlineage_GLM_AK_BC.tsv"),
    sep = "\t"
  )

  region_sum <- res_region[, .(
    n_tested_snps = .N,
    n_significant = sum(significant_mt, na.rm = TRUE),
    n_large_effect = sum(large_effect_mt, na.rm = TRUE),
    n_candidate = sum(candidate_mt, na.rm = TRUE)
  ), by = region_followup]

  fwrite(
    region_sum,
    file.path(OUTDIR, "followup_regionSpecific_mtlineage_GLM_summary.tsv"),
    sep = "\t"
  )
}

# ============================================================
# Summary print
# ============================================================

cat("\n================ Main OXPHOS72 GLM summary ================\n")

cat("\nGeography model: Latitude + Longitude\n")
cat("SNPs tested:", nrow(res_geo), "\n")
cat("Significant q < 0.05:", sum(res_geo$significant_geo, na.rm = TRUE), "\n")
cat("Large effect fitted AF range >= 0.2:", sum(res_geo$large_effect_geo, na.rm = TRUE), "\n")
cat("Candidates:", sum(res_geo$candidate_geo, na.rm = TRUE), "\n")

cat("\nmt-lineage model: all populations combined\n")
cat("SNPs tested:", nrow(res_mt), "\n")
cat("Significant q < 0.05:", sum(res_mt$significant_mt, na.rm = TRUE), "\n")
cat("Large effect max lineage AF diff >= 0.2:", sum(res_mt$large_effect_mt, na.rm = TRUE), "\n")
cat("Candidates:", sum(res_mt$candidate_mt, na.rm = TRUE), "\n")

cat("\nCandidate overlap:\n")
print(overlap_summary)

if (exists("region_sum")) {
  cat("\nRegion-specific follow-up:\n")
  print(region_sum)
}

# ============================================================
# Plot A: geography volcano
# ============================================================

res_geo[, plot_group := "Not candidate"]
res_geo[significant_geo == TRUE & large_effect_geo == FALSE, plot_group := "FDR only"]
res_geo[significant_geo == FALSE & large_effect_geo == TRUE, plot_group := "Large effect only"]
res_geo[candidate_geo == TRUE, plot_group := "Candidate"]

res_geo[, plot_group := factor(
  plot_group,
  levels = c("Not candidate", "FDR only", "Large effect only", "Candidate")
)]

TOP_GEO <- res_geo[candidate_geo == TRUE][order(q_geo, -geo_fitted_AF_range)]
TOP_GEO <- TOP_GEO[1:min(.N, 15)]

res_geo[, neglog10_q_geo_plot := pmin(neglog10_q_geo, 35)]
TOP_GEO[, neglog10_q_geo_plot := pmin(neglog10_q_geo, 35)]

p_geo <- ggplot(res_geo, aes(x = geo_fitted_AF_range, y = neglog10_q_geo_plot)) +
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
    x = "Fitted allele-frequency range across latitude/longitude",
    y = expression(-log[10]("FDR q-value"))
  ) +
  theme_fig +
  theme(legend.position = "right")

ggsave(
  file.path(OUTDIR, "Fig4A_geography_LatLong_OXPHOS72_GLM_volcano.png"),
  p_geo,
  width = 7.2,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig4A_geography_LatLong_OXPHOS72_GLM_volcano.pdf"),
  p_geo,
  width = 7.2,
  height = 5.5
)

# ============================================================
# Plot B: mt-lineage volcano, all populations combined
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
  file.path(OUTDIR, "Fig4B_mtlineage_allPops_OXPHOS72_GLM_volcano.png"),
  p_mt,
  width = 7.2,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig4B_mtlineage_allPops_OXPHOS72_GLM_volcano.pdf"),
  p_mt,
  width = 7.2,
  height = 5.5
)

# ============================================================
# Main figure AB only
# ============================================================

figAB <- p_geo + p_mt + plot_annotation(tag_levels = "A")

ggsave(
  file.path(OUTDIR, "Fig4AB_geoLatLong_mtlineage_main_ABonly.png"),
  figAB,
  width = 14.5,
  height = 5.8,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig4AB_geoLatLong_mtlineage_main_ABonly.pdf"),
  figAB,
  width = 14.5,
  height = 5.8
)

# ============================================================
# Former panel C saved separately
# Candidate overlap Venn-style
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

p_overlap <- ggplot() +
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
  file.path(OUTDIR, "Fig4C_candidate_overlap_separate_geoLatLong_mtlineage.png"),
  p_overlap,
  width = 5.8,
  height = 4.6,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig4C_candidate_overlap_separate_geoLatLong_mtlineage.pdf"),
  p_overlap,
  width = 5.8,
  height = 4.6
)

cat("\nDONE\n")
cat("Output dir:\n", OUTDIR, "\n")
