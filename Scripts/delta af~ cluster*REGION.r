delta af~ cluster*REGION
/mnt/spareHD_2/nu_287/q2_parallelism/09_mtClusterRegionInteraction_manualCluster_deltaAF_SNPlevel_LM.R
#!/usr/bin/env Rscript
# ============================================================
# 09_mtClusterRegionInteraction_manualCluster_deltaAF_SNPlevel_LM.R
#
# SNP-level regression across AK + BC together:
#   deltaAF ~ mtCluster(manual) * region + treePC1 + treePC2
#
# Main test of interest:
#   interaction term = mtCluster:region
# via nested-model F-test:
#   reduced: deltaAF ~ mtCluster + region + treePC1 + treePC2
#   full   : deltaAF ~ mtCluster * region + treePC1 + treePC2
#
# Inputs:
#   DELTA_FILE:
#     /mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz
#   CLUSTER_FILE:
#     /mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv
#
# Output:
#   /mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtClusterRegionInteraction_manualCluster_noAMO/
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ----------------------------
# Config
# ----------------------------
DELTA_FILE <- Sys.getenv(
  "DELTA",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
)

CLUSTER_FILE <- Sys.getenv(
  "CLUSTER",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
)

OUT_DIR <- Sys.getenv(
  "OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtClusterRegionInteraction_manualCluster_noAMO"
)

MIN_N_TOTAL <- as.integer(Sys.getenv("MIN_N", unset = "8"))
MIN_N_AK    <- as.integer(Sys.getenv("MIN_N_AK", unset = "3"))
MIN_N_BC    <- as.integer(Sys.getenv("MIN_N_BC", unset = "3"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl = TRUE)
  x
}

safe_zread <- function(f){
  if(!file.exists(f)) stop("File not found: ", f)
  fread(cmd = paste("zcat", shQuote(f)))
}

cluster_region_mean_str <- function(dt){
  tmp <- dt[, .(
    mu = mean(deltaAF, na.rm = TRUE),
    n  = .N
  ), by = .(region, cluster)][order(region, cluster)]

  if(nrow(tmp) == 0) return(NA_character_)

  paste0(
    tmp$region, ":", tmp$cluster, ":",
    sprintf("%.6g", tmp$mu), "(n=", tmp$n, ")",
    collapse = ";"
  )
}

safe_anova_p <- function(fit1, fit2){
  a <- try(anova(fit1, fit2), silent = TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2) return(NA_real_)
  as.numeric(a$`Pr(>F)`[2])
}

safe_anova_F <- function(fit1, fit2){
  a <- try(anova(fit1, fit2), silent = TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2) return(NA_real_)
  as.numeric(a$F[2])
}

# ----------------------------
# Core function
# ----------------------------
fit_interaction_model <- function(dt){

  dt <- copy(dt)
  dt <- dt[
    is.finite(deltaAF) &
    is.finite(treePC1) &
    is.finite(treePC2) &
    !is.na(region) &
    !is.na(cluster)
  ]

  n_all <- nrow(dt)
  n_AK  <- sum(dt$region == "AK")
  n_BC  <- sum(dt$region == "BC")

  out_na <- list(
    n = n_all,
    n_AK = n_AK,
    n_BC = n_BC,
    n_clusters = uniqueN(dt$cluster),
    F_interaction = NA_real_,
    p_interaction = NA_real_,
    df1 = NA_real_,
    df2 = NA_real_
  )

  if(n_all < MIN_N_TOTAL) return(out_na)
  if(n_AK < MIN_N_AK || n_BC < MIN_N_BC) return(out_na)
  if(length(unique(dt$region)) < 2) return(out_na)
  if(length(unique(dt$cluster)) < 2) return(out_na)

  dt[, region := factor(region, levels = c("AK", "BC"))]
  dt[, cluster := factor(cluster)]

  fit_add <- try(
    lm(deltaAF ~ cluster + region + treePC1 + treePC2, data = dt),
    silent = TRUE
  )

  fit_full <- try(
    lm(deltaAF ~ cluster * region + treePC1 + treePC2, data = dt),
    silent = TRUE
  )

  if(inherits(fit_add, "try-error") || inherits(fit_full, "try-error")){
    return(out_na)
  }

  a <- try(anova(fit_add, fit_full), silent = TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2){
    return(out_na)
  }

  list(
    n = n_all,
    n_AK = n_AK,
    n_BC = n_BC,
    n_clusters = uniqueN(dt$cluster),
    F_interaction = as.numeric(a$F[2]),
    p_interaction = as.numeric(a$`Pr(>F)`[2]),
    df1 = as.numeric(a$Df[2]),
    df2 = as.numeric(a$Res.Df[2])
  )
}

# ----------------------------
# Load deltaAF
# ----------------------------
cat("[read] deltaAF_long: ", DELTA_FILE, "\n", sep = "")
DEL <- safe_zread(DELTA_FILE)

req <- c("region","snp","chr","pos","gene","pop","deltaAF","treePC1","treePC2")
miss <- setdiff(req, names(DEL))
if(length(miss) > 0){
  cat("\n[ERROR] delta file missing columns:\n  ", paste(miss, collapse = ", "), "\n", sep = "")
  cat("\n[FOUND columns]\n")
  print(names(DEL))
  stop("delta file format mismatch.")
}

DEL[, pop := normalize_pop(pop)]
DEL <- DEL[pop != "AMO"]
DEL <- DEL[region %in% c("AK", "BC")]
DEL <- DEL[is.finite(deltaAF) & is.finite(treePC1) & is.finite(treePC2)]

cat("[info] rows after filters: ", nrow(DEL), "\n", sep = "")
cat("[info] pops by region:\n")
print(unique(DEL[, .(pop, region)])[ , .N, by = region][order(region)])

# ----------------------------
# Load manual cluster
# ----------------------------
cat("[read] mtCluster_manual: ", CLUSTER_FILE, "\n", sep = "")
CL <- fread(CLUSTER_FILE, sep = "\t", header = TRUE, fill = TRUE)

CL <- CL[!(is.na(pop) | pop == "")]
need_cl <- c("pop","mtCluster")
miss_cl <- setdiff(need_cl, names(CL))
if(length(miss_cl) > 0){
  stop("Cluster file missing columns: ", paste(miss_cl, collapse = ", "))
}

CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL[, .(pop, mtCluster)], by = "pop")
CL[, cluster := factor(mtCluster)]

cat("[info] clusters loaded:\n")
print(CL[, .N, by = cluster][order(-N, cluster)])

DEL <- merge(DEL, CL[, .(pop, cluster)], by = "pop", all.x = TRUE)

cat("[info] pops missing cluster: ", sum(is.na(DEL$cluster)), "\n", sep = "")
if(sum(is.na(DEL$cluster)) > 0){
  cat("[WARN] pops without cluster:\n")
  print(unique(DEL[is.na(cluster), pop]))
}

DEL <- DEL[!is.na(cluster)]
DEL[, region := factor(region, levels = c("AK", "BC"))]

# save pop-cluster map
pop_cluster_out <- unique(DEL[, .(pop, region, cluster)], by = "pop")
fwrite(
  pop_cluster_out[order(region, cluster, pop)],
  file.path(OUT_DIR, "mtCluster_by_pop.manual.noAMO.tsv"),
  sep = "\t"
)
cat("[write] mtCluster_by_pop.manual.noAMO.tsv\n")

# ----------------------------
# Per-SNP run
# ----------------------------
cat("\n[info] rows by region:\n")
print(DEL[, .N, by = region][order(region)])

snps <- unique(DEL$snp)
cat("[info] total SNPs: ", length(snps), "\n", sep = "")

setkey(DEL, snp)

res_list <- vector("list", length(snps))
idx <- 0L

for(s in snps){
  dt <- DEL[list(s)]

  chr_  <- dt$chr[1]
  pos_  <- dt$pos[1]
  gene_ <- dt$gene[1]

  out <- fit_interaction_model(dt)

  idx <- idx + 1L
  res_list[[idx]] <- data.table(
    snp = s,
    chr = chr_,
    pos = pos_,
    gene = gene_,
    n = out$n,
    n_AK = out$n_AK,
    n_BC = out$n_BC,
    n_clusters = out$n_clusters,
    F_interaction = out$F_interaction,
    p_interaction = out$p_interaction,
    df1 = out$df1,
    df2 = out$df2,
    cluster_region_means = cluster_region_mean_str(dt)
  )
}

RES <- rbindlist(res_list[seq_len(idx)], use.names = TRUE, fill = TRUE)

out_main <- file.path(OUT_DIR, "LM_perSNP_mtClusterRegionInteraction_manualCluster_plus_treePC12.tsv.gz")
fwrite(RES, out_main, sep = "\t")
cat("[write] ", out_main, "\n", sep = "")

# ----------------------------
# Summary
# ----------------------------
SUM <- RES[, .(
  n_snps = .N,
  prop_p05 = mean(p_interaction < 0.05, na.rm = TRUE),
  prop_p10 = mean(p_interaction < 0.10, na.rm = TRUE),
  min_p = suppressWarnings(min(p_interaction, na.rm = TRUE)),
  med_p = suppressWarnings(median(p_interaction, na.rm = TRUE))
)]

sum_out <- file.path(OUT_DIR, "LM_interaction_summary.tsv")
fwrite(SUM, sum_out, sep = "\t")
cat("[write] LM_interaction_summary.tsv\n\n")

cat("=== summary (cluster:region interaction) ===\n")
print(SUM)

# ----------------------------
# Top SNPs
# ----------------------------
cat("\n=== Top 20 SNPs by p_interaction ===\n")
TOP <- RES[is.finite(p_interaction)][order(p_interaction)][1:min(20, .N)]
print(TOP[, .(
  snp, gene, chr, pos, n, n_AK, n_BC, n_clusters,
  F_interaction, p_interaction
)])

# ----------------------------
# Plots
# ----------------------------
P1 <- ggplot(RES[is.finite(p_interaction)], aes(p_interaction)) +
  geom_histogram(bins = 50) +
  theme_classic(base_size = 14) +
  labs(
    x = "p-value for cluster:region interaction",
    y = "count",
    title = "deltaAF ~ manual mtCluster * region + treePC1 + treePC2"
  )

ggsave(
  file.path(OUT_DIR, "Fig_pHist_mtClusterRegionInteraction_manual.png"),
  P1, width = 7, height = 5.5, dpi = 300
)

QQ <- RES[is.finite(p_interaction) & p_interaction > 0 & p_interaction <= 1]
if(nrow(QQ) > 0){
  QQ <- copy(QQ)
  QQ[, obs := -log10(sort(p_interaction))]
  QQ[, exp := -log10(ppoints(.N))]

  P2 <- ggplot(QQ, aes(exp, obs)) +
    geom_point(alpha = 0.35, size = 0.8) +
    geom_abline(slope = 1, intercept = 0) +
    theme_classic(base_size = 14) +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = "QQ plot: manual mtCluster × region interaction"
    )

  ggsave(
    file.path(OUT_DIR, "Fig_QQ_mtClusterRegionInteraction_manual.png"),
    P2, width = 6.2, height = 4.8, dpi = 300
  )
}

cat("\n[OK] done. Outputs in:\n  ", OUT_DIR, "\n", sep = "")




/mnt/spareHD_2/nu_287/q2_parallelism/09_mtClusterRegion_variancePartition.R
#!/usr/bin/env Rscript
# ============================================================
# 09_mtClusterManual_region_variancePartition.R
#
# Per-SNP variance partition across AK + BC:
#
#   M0    : deltaAF ~ treePC1 + treePC2
#   M_reg : deltaAF ~ region  + treePC1 + treePC2
#   M_clu : deltaAF ~ cluster + treePC1 + treePC2
#   M_full: deltaAF ~ region + cluster + treePC1 + treePC2
#
# Goal:
#   compare external driver (region) vs internal driver (manual mtCluster)
#
# Inputs:
#   DELTA_FILE:
#     /mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz
#   CLUSTER_FILE:
#     /mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv
#
# Output:
#   /mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_variancePartition_manualCluster_noAMO/
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ----------------------------
# Config
# ----------------------------
DELTA_FILE <- Sys.getenv(
  "DELTA",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
)

CLUSTER_FILE <- Sys.getenv(
  "CLUSTER",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
)

OUT_DIR <- Sys.getenv(
  "OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_variancePartition_manualCluster_noAMO"
)

MIN_N_TOTAL <- as.integer(Sys.getenv("MIN_N", unset = "8"))
MIN_N_AK    <- as.integer(Sys.getenv("MIN_N_AK", unset = "3"))
MIN_N_BC    <- as.integer(Sys.getenv("MIN_N_BC", unset = "3"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl = TRUE)
  x
}

safe_zread <- function(f){
  if(!file.exists(f)) stop("File not found: ", f)
  fread(cmd = paste("zcat", shQuote(f)))
}

safe_r2 <- function(fit){
  out <- try(summary(fit)$r.squared, silent = TRUE)
  if(inherits(out, "try-error")) return(NA_real_)
  as.numeric(out)
}

safe_adj_r2 <- function(fit){
  out <- try(summary(fit)$adj.r.squared, silent = TRUE)
  if(inherits(out, "try-error")) return(NA_real_)
  as.numeric(out)
}

safe_anova_p <- function(fit1, fit2){
  a <- try(anova(fit1, fit2), silent = TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2) return(NA_real_)
  as.numeric(a$`Pr(>F)`[2])
}

safe_anova_F <- function(fit1, fit2){
  a <- try(anova(fit1, fit2), silent = TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2) return(NA_real_)
  as.numeric(a$F[2])
}

cluster_region_mean_str <- function(dt){
  tmp <- dt[, .(
    mu = mean(deltaAF, na.rm = TRUE),
    n  = .N
  ), by = .(region, cluster)][order(region, cluster)]

  if(nrow(tmp) == 0) return(NA_character_)

  paste0(
    tmp$region, ":", tmp$cluster, ":",
    sprintf("%.6g", tmp$mu), "(n=", tmp$n, ")",
    collapse = ";"
  )
}

# ----------------------------
# Core function
# ----------------------------
fit_partition_models <- function(dt){

  dt <- copy(dt)
  dt <- dt[
    is.finite(deltaAF) &
    is.finite(treePC1) &
    is.finite(treePC2) &
    !is.na(region) &
    !is.na(cluster)
  ]

  n_all <- nrow(dt)
  n_AK  <- sum(dt$region == "AK")
  n_BC  <- sum(dt$region == "BC")

  out_na <- list(
    n = n_all,
    n_AK = n_AK,
    n_BC = n_BC,
    n_clusters = uniqueN(dt$cluster),

    R2_M0 = NA_real_,
    R2_reg = NA_real_,
    R2_clu = NA_real_,
    R2_full = NA_real_,

    adjR2_M0 = NA_real_,
    adjR2_reg = NA_real_,
    adjR2_clu = NA_real_,
    adjR2_full = NA_real_,

    contrib_region = NA_real_,
    contrib_cluster = NA_real_,
    shared_like = NA_real_,

    F_reg_vs_M0 = NA_real_,
    p_reg_vs_M0 = NA_real_,

    F_clu_vs_M0 = NA_real_,
    p_clu_vs_M0 = NA_real_,

    F_region_unique = NA_real_,
    p_region_unique = NA_real_,

    F_cluster_unique = NA_real_,
    p_cluster_unique = NA_real_
  )

  if(n_all < MIN_N_TOTAL) return(out_na)
  if(n_AK < MIN_N_AK || n_BC < MIN_N_BC) return(out_na)
  if(length(unique(dt$region)) < 2) return(out_na)
  if(length(unique(dt$cluster)) < 2) return(out_na)

  dt[, region := factor(region, levels = c("AK", "BC"))]
  dt[, cluster := factor(cluster)]

  M0 <- try(lm(deltaAF ~ treePC1 + treePC2, data = dt), silent = TRUE)
  M_reg <- try(lm(deltaAF ~ region + treePC1 + treePC2, data = dt), silent = TRUE)
  M_clu <- try(lm(deltaAF ~ cluster + treePC1 + treePC2, data = dt), silent = TRUE)
  M_full <- try(lm(deltaAF ~ region + cluster + treePC1 + treePC2, data = dt), silent = TRUE)

  if(inherits(M0, "try-error") ||
     inherits(M_reg, "try-error") ||
     inherits(M_clu, "try-error") ||
     inherits(M_full, "try-error")){
    return(out_na)
  }

  R2_M0   <- safe_r2(M0)
  R2_reg  <- safe_r2(M_reg)
  R2_clu  <- safe_r2(M_clu)
  R2_full <- safe_r2(M_full)

  adjR2_M0   <- safe_adj_r2(M0)
  adjR2_reg  <- safe_adj_r2(M_reg)
  adjR2_clu  <- safe_adj_r2(M_clu)
  adjR2_full <- safe_adj_r2(M_full)

  contrib_region  <- R2_full - R2_clu
  contrib_cluster <- R2_full - R2_reg
  shared_like     <- R2_full - R2_M0 - contrib_region - contrib_cluster

  F_reg_vs_M0 <- safe_anova_F(M0, M_reg)
  p_reg_vs_M0 <- safe_anova_p(M0, M_reg)

  F_clu_vs_M0 <- safe_anova_F(M0, M_clu)
  p_clu_vs_M0 <- safe_anova_p(M0, M_clu)

  F_region_unique <- safe_anova_F(M_clu, M_full)
  p_region_unique <- safe_anova_p(M_clu, M_full)

  F_cluster_unique <- safe_anova_F(M_reg, M_full)
  p_cluster_unique <- safe_anova_p(M_reg, M_full)

  list(
    n = n_all,
    n_AK = n_AK,
    n_BC = n_BC,
    n_clusters = uniqueN(dt$cluster),

    R2_M0 = R2_M0,
    R2_reg = R2_reg,
    R2_clu = R2_clu,
    R2_full = R2_full,

    adjR2_M0 = adjR2_M0,
    adjR2_reg = adjR2_reg,
    adjR2_clu = adjR2_clu,
    adjR2_full = adjR2_full,

    contrib_region = contrib_region,
    contrib_cluster = contrib_cluster,
    shared_like = shared_like,

    F_reg_vs_M0 = F_reg_vs_M0,
    p_reg_vs_M0 = p_reg_vs_M0,

    F_clu_vs_M0 = F_clu_vs_M0,
    p_clu_vs_M0 = p_clu_vs_M0,

    F_region_unique = F_region_unique,
    p_region_unique = p_region_unique,

    F_cluster_unique = F_cluster_unique,
    p_cluster_unique = p_cluster_unique
  )
}

# ----------------------------
# Load deltaAF
# ----------------------------
cat("[read] deltaAF_long: ", DELTA_FILE, "\n", sep = "")
DEL <- safe_zread(DELTA_FILE)

req <- c("region","snp","chr","pos","gene","pop","deltaAF","treePC1","treePC2")
miss <- setdiff(req, names(DEL))
if(length(miss) > 0){
  cat("\n[ERROR] delta file missing columns:\n  ", paste(miss, collapse = ", "), "\n", sep = "")
  cat("\n[FOUND columns]\n")
  print(names(DEL))
  stop("delta file format mismatch.")
}

DEL[, pop := normalize_pop(pop)]
DEL <- DEL[pop != "AMO"]
DEL <- DEL[region %in% c("AK","BC")]
DEL <- DEL[is.finite(deltaAF) & is.finite(treePC1) & is.finite(treePC2)]

cat("[info] rows after filters: ", nrow(DEL), "\n", sep = "")
cat("[info] pops by region:\n")
print(unique(DEL[, .(pop, region)])[ , .N, by = region][order(region)])

# ----------------------------
# Load manual cluster
# ----------------------------
cat("[read] mtCluster_manual: ", CLUSTER_FILE, "\n", sep = "")
CL <- fread(CLUSTER_FILE, sep = "\t", header = TRUE, fill = TRUE)

CL <- CL[!(is.na(pop) | pop == "")]
need_cl <- c("pop","mtCluster")
miss_cl <- setdiff(need_cl, names(CL))
if(length(miss_cl) > 0){
  stop("Cluster file missing columns: ", paste(miss_cl, collapse = ", "))
}

CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL[, .(pop, mtCluster)], by = "pop")
CL[, cluster := factor(mtCluster)]

cat("[info] clusters loaded:\n")
print(CL[, .N, by = cluster][order(-N, cluster)])

DEL <- merge(DEL, CL[, .(pop, cluster)], by = "pop", all.x = TRUE)

cat("[info] pops missing cluster: ", sum(is.na(DEL$cluster)), "\n", sep = "")
if(sum(is.na(DEL$cluster)) > 0){
  cat("[WARN] pops without cluster:\n")
  print(unique(DEL[is.na(cluster), pop]))
}

DEL <- DEL[!is.na(cluster)]
DEL[, region := factor(region, levels = c("AK", "BC"))]

# save pop->cluster map
pop_cluster_out <- unique(DEL[, .(pop, region, cluster)], by = "pop")
fwrite(
  pop_cluster_out[order(region, cluster, pop)],
  file.path(OUT_DIR, "mtCluster_by_pop.manual.noAMO.tsv"),
  sep = "\t"
)
cat("[write] mtCluster_by_pop.manual.noAMO.tsv\n")

# ----------------------------
# Per-SNP run
# ----------------------------
cat("\n[info] rows by region:\n")
print(DEL[, .N, by = region][order(region)])

snps <- unique(DEL$snp)
cat("[info] total SNPs: ", length(snps), "\n", sep = "")

setkey(DEL, snp)

res_list <- vector("list", length(snps))
idx <- 0L

for(s in snps){
  dt <- DEL[list(s)]

  chr_  <- dt$chr[1]
  pos_  <- dt$pos[1]
  gene_ <- dt$gene[1]

  out <- fit_partition_models(dt)

  idx <- idx + 1L
  res_list[[idx]] <- data.table(
    snp = s,
    chr = chr_,
    pos = pos_,
    gene = gene_,
    n = out$n,
    n_AK = out$n_AK,
    n_BC = out$n_BC,
    n_clusters = out$n_clusters,

    R2_M0 = out$R2_M0,
    R2_reg = out$R2_reg,
    R2_clu = out$R2_clu,
    R2_full = out$R2_full,

    adjR2_M0 = out$adjR2_M0,
    adjR2_reg = out$adjR2_reg,
    adjR2_clu = out$adjR2_clu,
    adjR2_full = out$adjR2_full,

    contrib_region = out$contrib_region,
    contrib_cluster = out$contrib_cluster,
    shared_like = out$shared_like,

    F_reg_vs_M0 = out$F_reg_vs_M0,
    p_reg_vs_M0 = out$p_reg_vs_M0,

    F_clu_vs_M0 = out$F_clu_vs_M0,
    p_clu_vs_M0 = out$p_clu_vs_M0,

    F_region_unique = out$F_region_unique,
    p_region_unique = out$p_region_unique,

    F_cluster_unique = out$F_cluster_unique,
    p_cluster_unique = out$p_cluster_unique,

    cluster_region_means = cluster_region_mean_str(dt)
  )
}

RES <- rbindlist(res_list[seq_len(idx)], use.names = TRUE, fill = TRUE)

out_main <- file.path(OUT_DIR, "LM_perSNP_variancePartition_manualCluster_region_plus_treePC12.tsv.gz")
fwrite(RES, out_main, sep = "\t")
cat("[write] ", out_main, "\n", sep = "")

# ----------------------------
# Summary
# ----------------------------
SUM <- RES[, .(
  n_snps = .N,

  mean_R2_M0 = mean(R2_M0, na.rm = TRUE),
  mean_R2_reg = mean(R2_reg, na.rm = TRUE),
  mean_R2_clu = mean(R2_clu, na.rm = TRUE),
  mean_R2_full = mean(R2_full, na.rm = TRUE),

  mean_contrib_region = mean(contrib_region, na.rm = TRUE),
  mean_contrib_cluster = mean(contrib_cluster, na.rm = TRUE),

  median_contrib_region = median(contrib_region, na.rm = TRUE),
  median_contrib_cluster = median(contrib_cluster, na.rm = TRUE),

  prop_region_total_p05 = mean(p_reg_vs_M0 < 0.05, na.rm = TRUE),
  prop_cluster_total_p05 = mean(p_clu_vs_M0 < 0.05, na.rm = TRUE),

  prop_region_unique_p05 = mean(p_region_unique < 0.05, na.rm = TRUE),
  prop_cluster_unique_p05 = mean(p_cluster_unique < 0.05, na.rm = TRUE)
)]

sum_out <- file.path(OUT_DIR, "LM_variancePartition_summary.tsv")
fwrite(SUM, sum_out, sep = "\t")
cat("[write] LM_variancePartition_summary.tsv\n\n")

cat("=== summary ===\n")
print(SUM)

# ----------------------------
# Top SNPs
# ----------------------------
cat("\n=== Top 20 SNPs by unique cluster contribution ===\n")
TOP_CLU <- RES[is.finite(contrib_cluster)][order(-contrib_cluster)][1:min(20, .N)]
print(TOP_CLU[, .(
  snp, gene, chr, pos, n, n_AK, n_BC,
  contrib_cluster, contrib_region, R2_full, p_cluster_unique, p_region_unique
)])

cat("\n=== Top 20 SNPs by unique region contribution ===\n")
TOP_REG <- RES[is.finite(contrib_region)][order(-contrib_region)][1:min(20, .N)]
print(TOP_REG[, .(
  snp, gene, chr, pos, n, n_AK, n_BC,
  contrib_region, contrib_cluster, R2_full, p_region_unique, p_cluster_unique
)])

# ----------------------------
# Plots
# ----------------------------
P1 <- ggplot(
  RES[is.finite(contrib_region) & is.finite(contrib_cluster)],
  aes(contrib_region, contrib_cluster)
) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Unique contribution of region (R2_full - R2_clu)",
    y = "Unique contribution of cluster (R2_full - R2_reg)",
    title = "Per-SNP variance partition: region vs manual mtCluster"
  )

ggsave(
  file.path(OUT_DIR, "Fig_scatter_contribRegion_vs_contribCluster.png"),
  P1, width = 6.8, height = 5.8, dpi = 300
)

BOXDT <- rbind(
  RES[, .(component = "region", value = contrib_region)],
  RES[, .(component = "cluster", value = contrib_cluster)]
)
BOXDT <- BOXDT[is.finite(value)]

P2 <- ggplot(BOXDT, aes(component, value)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "Unique contribution (ΔR²)",
    title = "Distribution of unique contributions across SNPs"
  )

ggsave(
  file.path(OUT_DIR, "Fig_boxplot_uniqueContrib_region_vs_cluster.png"),
  P2, width = 5.8, height = 5.2, dpi = 300
)

# ----------------------------
# Simple win counts
# ----------------------------
WIN <- RES[is.finite(contrib_region) & is.finite(contrib_cluster), .(
  n_region_gt_cluster = sum(contrib_region > contrib_cluster),
  n_cluster_gt_region = sum(contrib_cluster > contrib_region),
  n_equal = sum(contrib_cluster == contrib_region),
  prop_region_gt_cluster = mean(contrib_region > contrib_cluster),
  prop_cluster_gt_region = mean(contrib_cluster > contrib_region)
)]

fwrite(WIN, file.path(OUT_DIR, "LM_variancePartition_winCounts.tsv"), sep = "\t")
cat("[write] LM_variancePartition_winCounts.tsv\n")

cat("\n=== win counts ===\n")
print(WIN)

cat("\n[OK] done. Outputs in:\n  ", OUT_DIR, "\n", sep = "")




#0.2 ld la+long
mean_unique_cluster median_unique_geo median_unique_cluster
                 <num>             <num>                 <num>
1:          0.09142093           0.20492            0.04321348
   prop_geo_p05_total prop_cluster_p05_total prop_geo_p05_unique
                <num>                  <num>               <num>
1:          0.4456812             0.03369428           0.4345724
   prop_cluster_p05_unique prop_geo_larger prop_cluster_larger
                     <num>           <num>               <num>
1:              0.02538403       0.6774112           0.3225888

[done] outputs in:cat
mean_unique_cluster: command not found
-bash: syntax error near unexpected token `<'
1:: command not found
prop_geo_p05_total: command not found
-bash: syntax error near unexpected token `<'
1:: command not found
prop_cluster_p05_unique: command not found
-bash: syntax error near unexpected token `<'
1:: command not found
[done]: command not found
(r_env) cyu@stickleback:~$ cat /mnt/spareHD_2/nu_287/q2_parallelism/09_geo_mtCluster_variancePartition_LDpruned_r2_0.2.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Geographic + mtCluster variance partition
# LD-pruned r2 = 0.2 version
#
# Models per SNP:
#   M0      : deltaAF ~ treePC1 + treePC2
#   M_geo   : deltaAF ~ Latitude_z + Longitude_z + treePC1 + treePC2
#   M_clu   : deltaAF ~ cluster + treePC1 + treePC2
#   M_full  : deltaAF ~ Latitude_z + Longitude_z + cluster + treePC1 + treePC2
#
# Main questions:
#   1. Does geography explain deltaAF?
#   2. Does mtCluster explain additional variance after geography?
# ============================================================

# ----------------------------
# Config
# ----------------------------
DELTA_FILE <- "/work/cyu/ldx_all_subunits/ld/bayenv_style_r2_0.2/deltaAF_long.noAMO.bayenvStyle_r2_0.2_ldPruned_kept.tsv.gz"

POOLINFO <- "/work/cyu/Poolinfo.csv"

CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"

OUT_DIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_geo_mtCluster_variancePartition_LDpruned_r2_0.2"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_N_TOTAL <- 8

# ----------------------------
# Helpers
# ----------------------------
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl = TRUE)
  x
}

safe_r2 <- function(fit){
  x <- try(summary(fit)$r.squared, silent = TRUE)
  if(inherits(x, "try-error")) return(NA_real_)
  as.numeric(x)
}

safe_adj_r2 <- function(fit){
  x <- try(summary(fit)$adj.r.squared, silent = TRUE)
  if(inherits(x, "try-error")) return(NA_real_)
  as.numeric(x)
}

safe_anova_p <- function(fit1, fit2){
  a <- try(anova(fit1, fit2), silent = TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2) return(NA_real_)
  as.numeric(a$`Pr(>F)`[2])
}

safe_anova_F <- function(fit1, fit2){
  a <- try(anova(fit1, fit2), silent = TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2) return(NA_real_)
  as.numeric(a$F[2])
}

fit_one_snp <- function(dt){

  dt <- copy(dt)

  dt <- dt[
    is.finite(deltaAF) &
      is.finite(treePC1) &
      is.finite(treePC2) &
      is.finite(Latitude_z) &
      is.finite(Longitude_z) &
      !is.na(cluster)
  ]

  out_na <- list(
    n = nrow(dt),
    n_clusters = uniqueN(dt$cluster),

    R2_M0 = NA_real_,
    R2_geo = NA_real_,
    R2_clu = NA_real_,
    R2_full = NA_real_,

    adjR2_M0 = NA_real_,
    adjR2_geo = NA_real_,
    adjR2_clu = NA_real_,
    adjR2_full = NA_real_,

    contrib_geo = NA_real_,
    contrib_cluster = NA_real_,
    shared_like = NA_real_,

    F_geo_vs_M0 = NA_real_,
    p_geo_vs_M0 = NA_real_,

    F_clu_vs_M0 = NA_real_,
    p_clu_vs_M0 = NA_real_,

    F_geo_unique = NA_real_,
    p_geo_unique = NA_real_,

    F_cluster_unique = NA_real_,
    p_cluster_unique = NA_real_
  )

  if(nrow(dt) < MIN_N_TOTAL) return(out_na)
  if(uniqueN(dt$cluster) < 2) return(out_na)

  dt[, cluster := factor(cluster)]

  M0 <- try(lm(deltaAF ~ treePC1 + treePC2, data = dt), silent = TRUE)

  M_geo <- try(
    lm(deltaAF ~ Latitude_z + Longitude_z + treePC1 + treePC2, data = dt),
    silent = TRUE
  )

  M_clu <- try(
    lm(deltaAF ~ cluster + treePC1 + treePC2, data = dt),
    silent = TRUE
  )

  M_full <- try(
    lm(deltaAF ~ Latitude_z + Longitude_z + cluster + treePC1 + treePC2, data = dt),
    silent = TRUE
  )

  if(
    inherits(M0, "try-error") |
    inherits(M_geo, "try-error") |
    inherits(M_clu, "try-error") |
    inherits(M_full, "try-error")
  ){
    return(out_na)
  }

  R2_M0   <- safe_r2(M0)
  R2_geo  <- safe_r2(M_geo)
  R2_clu  <- safe_r2(M_clu)
  R2_full <- safe_r2(M_full)

  contrib_geo     <- R2_full - R2_clu
  contrib_cluster <- R2_full - R2_geo
  shared_like     <- R2_full - R2_M0 - contrib_geo - contrib_cluster

  list(
    n = nrow(dt),
    n_clusters = uniqueN(dt$cluster),

    R2_M0 = R2_M0,
    R2_geo = R2_geo,
    R2_clu = R2_clu,
    R2_full = R2_full,

    adjR2_M0 = safe_adj_r2(M0),
    adjR2_geo = safe_adj_r2(M_geo),
    adjR2_clu = safe_adj_r2(M_clu),
    adjR2_full = safe_adj_r2(M_full),

    contrib_geo = contrib_geo,
    contrib_cluster = contrib_cluster,
    shared_like = shared_like,

    F_geo_vs_M0 = safe_anova_F(M0, M_geo),
    p_geo_vs_M0 = safe_anova_p(M0, M_geo),

    F_clu_vs_M0 = safe_anova_F(M0, M_clu),
    p_clu_vs_M0 = safe_anova_p(M0, M_clu),

    F_geo_unique = safe_anova_F(M_clu, M_full),
    p_geo_unique = safe_anova_p(M_clu, M_full),

    F_cluster_unique = safe_anova_F(M_geo, M_full),
    p_cluster_unique = safe_anova_p(M_geo, M_full)
  )
}

# ----------------------------
# Load data
# ----------------------------
cat("[read] deltaAF:", DELTA_FILE, "\n")
DEL <- fread(cmd = paste("zcat", shQuote(DELTA_FILE)))

DEL[, pop := normalize_pop(pop)]
DEL <- DEL[pop != "AMO"]

cat("[read] Poolinfo Excel:", POOLINFO, "\n")

suppressPackageStartupMessages({
  library(readxl)
})

INFO <- as.data.table(read_excel(POOLINFO))

setnames(INFO, old = "Population", new = "pop", skip_absent = TRUE)
INFO[, pop := normalize_pop(pop)]
INFO[, Latitude := as.numeric(Latitude)]
INFO[, Longitude := as.numeric(Longitude)]

INFO <- unique(INFO[, .(pop, Region, Latitude, Longitude, Habitat, Watershed)], by = "pop")

cat("[read] cluster:", CLUSTER_FILE, "\n")
CL <- fread(CLUSTER_FILE, sep = "\t", header = TRUE, fill = TRUE, blank.lines.skip = TRUE)
CL <- CL[!is.na(pop) & pop != "" & !is.na(mtCluster) & mtCluster != ""]
CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL[, .(pop, mtCluster)], by = "pop")
CL[, cluster := factor(mtCluster)]

# merge
DEL <- merge(DEL, INFO, by = "pop", all.x = TRUE)
DEL <- merge(DEL, CL[, .(pop, cluster)], by = "pop", all.x = TRUE)

# filters
DEL <- DEL[
  is.finite(deltaAF) &
    is.finite(treePC1) &
    is.finite(treePC2) &
    is.finite(Latitude) &
    is.finite(Longitude) &
    !is.na(cluster)
]

# scale geography
DEL[, Latitude_z := as.numeric(scale(Latitude))]
DEL[, Longitude_z := as.numeric(scale(Longitude))]

cat("[info] rows:", nrow(DEL), "\n")
cat("[info] pops:", uniqueN(DEL$pop), "\n")
cat("[info] SNPs:", uniqueN(DEL$snp), "\n")
cat("[info] cluster table:\n")
print(unique(DEL[, .(pop, Region, Latitude, Longitude, cluster)])[order(Region, cluster, pop)])

# save merged metadata
fwrite(
  unique(DEL[, .(pop, Region, Latitude, Longitude, Latitude_z, Longitude_z, Habitat, Watershed, cluster)])[order(Region, cluster, pop)],
  file.path(OUT_DIR, "pop_geo_mtCluster_used.tsv"),
  sep = "\t"
)

# ----------------------------
# Per-SNP variance partition
# ----------------------------
setkey(DEL, snp)
snps <- unique(DEL$snp)

cat("[run] per-SNP models...\n")

res_list <- vector("list", length(snps))

for(i in seq_along(snps)){
  if(i %% 1000 == 0) cat("[progress]", i, "/", length(snps), "\n")

  s <- snps[i]
  dt <- DEL[list(s)]

  out <- fit_one_snp(dt)

  res_list[[i]] <- data.table(
    snp = s,
    chr = dt$chr[1],
    pos = dt$pos[1],
    gene = dt$gene[1],

    n = out$n,
    n_clusters = out$n_clusters,

    R2_M0 = out$R2_M0,
    R2_geo = out$R2_geo,
    R2_clu = out$R2_clu,
    R2_full = out$R2_full,

    adjR2_M0 = out$adjR2_M0,
    adjR2_geo = out$adjR2_geo,
    adjR2_clu = out$adjR2_clu,
    adjR2_full = out$adjR2_full,

    contrib_geo = out$contrib_geo,
    contrib_cluster = out$contrib_cluster,
    shared_like = out$shared_like,

    F_geo_vs_M0 = out$F_geo_vs_M0,
    p_geo_vs_M0 = out$p_geo_vs_M0,

    F_clu_vs_M0 = out$F_clu_vs_M0,
    p_clu_vs_M0 = out$p_clu_vs_M0,

    F_geo_unique = out$F_geo_unique,
    p_geo_unique = out$p_geo_unique,

    F_cluster_unique = out$F_cluster_unique,
    p_cluster_unique = out$p_cluster_unique
  )
}

RES <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

RES[, q_geo_unique := p.adjust(p_geo_unique, method = "BH")]
RES[, q_cluster_unique := p.adjust(p_cluster_unique, method = "BH")]

out_main <- file.path(OUT_DIR, "variancePartition_geo_mtCluster_LDpruned_r2_0.2.tsv.gz")
fwrite(RES, out_main, sep = "\t")

cat("[write]", out_main, "\n")

# ----------------------------
# Summary table
# ----------------------------
SUM <- RES[, .(
  n_snps = .N,

  mean_R2_geo = mean(R2_geo, na.rm = TRUE),
  mean_R2_cluster = mean(R2_clu, na.rm = TRUE),
  mean_R2_full = mean(R2_full, na.rm = TRUE),

  mean_unique_geo = mean(contrib_geo, na.rm = TRUE),
  mean_unique_cluster = mean(contrib_cluster, na.rm = TRUE),

  median_unique_geo = median(contrib_geo, na.rm = TRUE),
  median_unique_cluster = median(contrib_cluster, na.rm = TRUE),

  prop_geo_p05_total = mean(p_geo_vs_M0 < 0.05, na.rm = TRUE),
  prop_cluster_p05_total = mean(p_clu_vs_M0 < 0.05, na.rm = TRUE),

  prop_geo_p05_unique = mean(p_geo_unique < 0.05, na.rm = TRUE),
  prop_cluster_p05_unique = mean(p_cluster_unique < 0.05, na.rm = TRUE),

  prop_geo_larger = mean(contrib_geo > contrib_cluster, na.rm = TRUE),
  prop_cluster_larger = mean(contrib_cluster > contrib_geo, na.rm = TRUE)
)]

fwrite(SUM, file.path(OUT_DIR, "summary_geo_mtCluster_LDpruned_r2_0.2.tsv"), sep = "\t")

cat("\n=== summary ===\n")
print(SUM)

# ----------------------------
# Gene-level summary
# ----------------------------
GENE_SUM <- RES[, .(
  n_snps = .N,
  mean_unique_geo = mean(contrib_geo, na.rm = TRUE),
  mean_unique_cluster = mean(contrib_cluster, na.rm = TRUE),
  prop_cluster_unique_p05 = mean(p_cluster_unique < 0.05, na.rm = TRUE),
  prop_geo_unique_p05 = mean(p_geo_unique < 0.05, na.rm = TRUE),
  min_p_cluster_unique = min(p_cluster_unique, na.rm = TRUE),
  min_p_geo_unique = min(p_geo_unique, na.rm = TRUE)
), by = gene][order(-mean_unique_cluster)]

fwrite(GENE_SUM, file.path(OUT_DIR, "gene_summary_geo_mtCluster_LDpruned_r2_0.2.tsv"), sep = "\t")

cat("\n[done] outputs in:\n", OUT_DIR, "\n")




#venn
library(VennDiagram)
library(grid)

grid.newpage()

venn.plot <- draw.triple.venn(
  area1 = total_oxphos,
  area2 = geo_total,
  area3 = cluster_total,
  
  n12 = geo_total,
  n13 = cluster_total,
  n23 = overlap,
  n123 = overlap,
  
  category = c("", "", ""),   # 关闭默认label
  
  scaled = FALSE,
  
  fill = c("grey92", "#9ecae1", "#fcae91"),
  alpha = c(0.45, 0.70, 0.70),
  
  col = c("grey40", "grey35", "grey35"),
  lwd = 2,
  
  cex = 1.8,
  fontface = "bold"
)

grid.draw(venn.plot)

# =========================
# 手动加文字
# x,y 范围都是 0~1
# =========================

grid.text(
  "OXPHOS SNPs",
  x = 0.5,
  y = 0.82,
  gp = gpar(fontsize = 20, fontface = "bold")
)

grid.text(
  "Geography",
  x = 0.75,
  y = 0.40,
  gp = gpar(fontsize = 20, fontface = "bold")
)

grid.text(
  "mtCluster",
  x = 0.25,
  y = 0.40,
  gp = gpar(fontsize = 20, fontface = "bold")
)

