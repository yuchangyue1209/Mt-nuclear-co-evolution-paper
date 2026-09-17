#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# 287-gene PBS with Norway outgroup
# Drop LB from AK freshwater populations
# ============================================================

# -------------------------
# paths
# -------------------------
MT_FST_FILE <- "/work/cyu/poolseq/PPalign_output/fst_mt_noDloop_withNorway/mtgenome_noDloop_withNorway_fst.csv"
NU_FST_DIR  <- "/work/cyu/gene_fst_work_withNorway/fst_out_all_withNorway"

OUTDIR <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup_dropLB"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# population design
# -------------------------
# AK freshwater: LB removed
AK_FW <- c("FG", "LG", "SR", "SL", "TL", "WB", "WT", "WK")

# BC freshwater: keep AMO here; remove AMO later for OXPHOS72 / figure analyses if needed
BC_FW <- c("SWA", "THE", "JOE", "BEA", "MUC", "PYE", "AMO",
           "BOOT", "ECHO", "LAW", "GOS", "ROB")

design <- rbindlist(list(
  data.table(focal = AK_FW, region = "AK", marine = "RS",  outgroup = "Norway"),
  data.table(focal = BC_FW, region = "BC", marine = "SAY", outgroup = "Norway")
))

design <- unique(design)

cat("\n[Design populations]\n")
print(design)

# ============================================================
# helper functions
# ============================================================

fst_to_matrix <- function(file) {
  df <- fread(file)

  if (nrow(df) < 1 || ncol(df) < 4) {
    stop("Bad FST file: ", file)
  }

  row <- df[1]

  pair_cols <- colnames(df)[4:ncol(df)]
  vals <- as.numeric(unlist(row[, 4:ncol(df), with = FALSE]))

  pair_names <- sub("\\.fst$", "", pair_cols)
  sp <- tstrsplit(pair_names, ":", fixed = TRUE)

  pop1 <- sp[[1]]
  pop2 <- sp[[2]]

  pops <- sort(unique(c(pop1, pop2)))

  mat <- matrix(
    NA_real_,
    nrow = length(pops),
    ncol = length(pops),
    dimnames = list(pops, pops)
  )

  diag(mat) <- 0

  for (i in seq_along(vals)) {
    mat[pop1[i], pop2[i]] <- vals[i]
    mat[pop2[i], pop1[i]] <- vals[i]
  }

  mat
}

get_fst <- function(mat, a, b) {
  if (!(a %in% rownames(mat)) || !(b %in% colnames(mat))) {
    return(NA_real_)
  }
  mat[a, b]
}

calc_pbs_one <- function(mat, focal, marine, outgroup) {
  fst_AB <- get_fst(mat, focal, marine)
  fst_AC <- get_fst(mat, focal, outgroup)
  fst_BC <- get_fst(mat, marine, outgroup)

  fst_AB <- pmin(pmax(fst_AB, 0), 0.999999)
  fst_AC <- pmin(pmax(fst_AC, 0), 0.999999)
  fst_BC <- pmin(pmax(fst_BC, 0), 0.999999)

  if (!all(is.finite(c(fst_AB, fst_AC, fst_BC)))) {
    return(data.table(
      fst_focal_marine = fst_AB,
      fst_focal_outgroup = fst_AC,
      fst_marine_outgroup = fst_BC,
      T_focal_marine = NA_real_,
      T_focal_outgroup = NA_real_,
      T_marine_outgroup = NA_real_,
      PBS = NA_real_
    ))
  }

  T_AB <- -log(1 - fst_AB)
  T_AC <- -log(1 - fst_AC)
  T_BC <- -log(1 - fst_BC)

  PBS <- (T_AB + T_AC - T_BC) / 2

  data.table(
    fst_focal_marine = fst_AB,
    fst_focal_outgroup = fst_AC,
    fst_marine_outgroup = fst_BC,
    T_focal_marine = T_AB,
    T_focal_outgroup = T_AC,
    T_marine_outgroup = T_BC,
    PBS = PBS
  )
}

calc_pbs_table <- function(mat, design_dt) {
  rbindlist(lapply(seq_len(nrow(design_dt)), function(i) {
    d <- design_dt[i]

    out <- calc_pbs_one(
      mat = mat,
      focal = d$focal,
      marine = d$marine,
      outgroup = d$outgroup
    )

    cbind(d, out)
  }))
}

# ============================================================
# 1. mitochondrial PBS
# ============================================================

cat("\n[1] Reading mt FST with Norway...\n")
mt_mat <- fst_to_matrix(MT_FST_FILE)

design_use <- design[
  focal %in% rownames(mt_mat) &
    marine %in% rownames(mt_mat) &
    outgroup %in% rownames(mt_mat)
]

cat("[Info] Design rows used for mt PBS:", nrow(design_use), "\n")
print(design_use)

mt_pbs <- calc_pbs_table(mt_mat, design_use)
setnames(mt_pbs, "PBS", "mt_PBS")

fwrite(
  mt_pbs,
  file.path(OUTDIR, "mt_noDloop_PBS_NorwayOutgroup_dropLB_by_population.tsv"),
  sep = "\t"
)

# ============================================================
# 2. nuclear gene PBS
# ============================================================

files <- list.files(
  NU_FST_DIR,
  pattern = "_fst\\.csv$",
  full.names = TRUE
)

cat("\n[2] Nuclear gene FST files found:", length(files), "\n")

nu_list <- vector("list", length(files))

for (i in seq_along(files)) {
  f <- files[i]
  gene <- sub("_fst\\.csv$", "", basename(f))

  cat("[run] ", i, "/", length(files), "  ", gene, "\n", sep = "")

  nu_list[[i]] <- tryCatch({
    mat <- fst_to_matrix(f)

    tmp <- calc_pbs_table(mat, design_use)
    tmp[, gene := gene]

    tmp
  }, error = function(e) {
    data.table(
      focal = NA_character_,
      region = NA_character_,
      marine = NA_character_,
      outgroup = NA_character_,
      gene = gene,
      fst_focal_marine = NA_real_,
      fst_focal_outgroup = NA_real_,
      fst_marine_outgroup = NA_real_,
      T_focal_marine = NA_real_,
      T_focal_outgroup = NA_real_,
      T_marine_outgroup = NA_real_,
      PBS = NA_real_,
      error = conditionMessage(e)
    )
  })
}

nu_pbs <- rbindlist(nu_list, use.names = TRUE, fill = TRUE)
setnames(nu_pbs, "PBS", "nu_PBS")

fwrite(
  nu_pbs,
  file.path(OUTDIR, "nuclear_gene_PBS_NorwayOutgroup_dropLB_by_population.tsv"),
  sep = "\t"
)

# ============================================================
# 3. merge mt and nuclear PBS
# ============================================================

merged <- merge(
  nu_pbs,
  mt_pbs[, .(focal, region, marine, outgroup, mt_PBS)],
  by = c("focal", "region", "marine", "outgroup"),
  all.x = TRUE
)

# Safety check: LB should not exist
merged <- merged[focal != "LB"]

fwrite(
  merged,
  file.path(OUTDIR, "merged_mt_nuclear_PBS_NorwayOutgroup_dropLB_by_gene_population.tsv"),
  sep = "\t"
)

# ============================================================
# 4. gene-wise correlation: nuclear PBS vs mt PBS
# ============================================================

gene_res <- merged[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))

      data.table(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    } else {
      data.table(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    }
  },
  by = gene
]

gene_res[, p_bh := p.adjust(p, method = "BH")]
setorder(gene_res, p, -cor_r)

fwrite(
  gene_res,
  file.path(OUTDIR, "gene_correlation_nuPBS_vs_mtPBS_NorwayOutgroup_dropLB.tsv"),
  sep = "\t"
)

# ============================================================
# 5. region-specific gene correlation
# ============================================================

gene_region_res <- merged[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 4 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))

      data.table(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    } else {
      data.table(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    }
  },
  by = .(gene, region)
]

gene_region_res[, p_bh := p.adjust(p, method = "BH"), by = region]
setorder(gene_region_res, region, p, -cor_r)

fwrite(
  gene_region_res,
  file.path(OUTDIR, "gene_region_correlation_nuPBS_vs_mtPBS_NorwayOutgroup_dropLB.tsv"),
  sep = "\t"
)

# ============================================================
# 6. quick volcano plot
# ============================================================

plot_dt <- gene_res[is.finite(cor_r) & is.finite(p) & p > 0]

if (nrow(plot_dt) > 0) {
  plot_dt[, logp := -log10(p)]

  p1 <- ggplot(plot_dt, aes(x = cor_r, y = logp)) +
    geom_point(size = 2.2, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    theme_classic(base_size = 14) +
    labs(
      x = "Correlation between nuclear gene PBS and mtDNA PBS",
      y = expression(-log[10](p)),
      title = "287 genes: nuclear PBS vs mtDNA PBS",
      subtitle = "Norway outgroup; LB removed"
    )

  ggsave(
    file.path(OUTDIR, "Fig_287genes_nuPBS_vs_mtPBS_NorwayOutgroup_dropLB_volcano.png"),
    p1,
    width = 7,
    height = 5,
    dpi = 300
  )

  ggsave(
    file.path(OUTDIR, "Fig_287genes_nuPBS_vs_mtPBS_NorwayOutgroup_dropLB_volcano.pdf"),
    p1,
    width = 7,
    height = 5
  )
}

# ============================================================
# 7. final checks
# ============================================================

cat("\n[Check] Populations in merged PBS table:\n")
print(sort(unique(merged$focal)))

cat("\n[Check] Number of genes:\n")
print(uniqueN(merged$gene))

cat("\n[Check] Rows by region:\n")
print(merged[, .N, by = region])

cat("\n[Check] Valid nuclear PBS rows by region:\n")
print(merged[is.finite(nu_PBS), .N, by = region])

cat("\n[Check] Top 20 gene correlations:\n")
print(head(gene_res, 20))

cat("\n[OK] 287-gene PBS with LB removed is done.\n")
cat("[Output directory]\n")
cat(OUTDIR, "\n")
