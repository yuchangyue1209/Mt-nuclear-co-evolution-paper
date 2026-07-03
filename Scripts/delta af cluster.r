delta af cluster
#delta /mnt/spareHD_2/nu_287/q2_parallelism/
AF long (原料)
  └── af_long_final_72genes_subunit_with_si.tsv.gz
        |
        v
(1) 08_q2_deltaAF_SNPlevel_vs_mtPCk_plus_treePC12.R
        ├── deltaAF_long.noAMO.tsv.gz          [核心中间文件]
        ├── LM_perSNP_mitoPCk_plus_treePC12.tsv.gz
        └── LM_perSNP_treeOnly.tsv.gz
        |
        v
(2) 09_mtCluster_deltaAF_SNPlevel_LM_perm.R    [tree-based cluster + permutation]
        ├── LM_perSNP_mtCluster_plus_treePC12.tsv.gz
        ├── Perm_perSNP_mtCluster_plus_treePC12.tsv.gz
        └── Fig_*.png + Perm_global_summary.tsv
        |
        v
(3) 09_mtClusterManual_deltaAF_SNPlevel_perm0.R [manual cluster + driver_cluster]
        ├── LM_perSNP_mtCluster_manual_plus_treePC12.tsv.gz
        ├── LM_driverCluster_share_sigSNP.tsv
        └── Fig_*.png
        |
        v
(4) 10_gene_enrichment_from_cluster_LM.R       [BH + gene enrichment]
        └── gene_enrichment/*
              ├── LM_with_q_and_driver.tsv.gz
              ├── Gene_enrichment_overall.tsv.gz
              ├── GeneCluster_enrichment.tsv.gz
              └── Fig_topGenes_byCluster_*.png

(5) 06_q2_deltaAF_meanSD_vs_mtPC_noAMO.R       [pop-level meanΔ/sdΔ summary]
        ├── LM_meanDelta_vs_mitoPC1to5_plus_treePC12.tsv
        └── LM_sdDelta_vs_mitoPC1to5_plus_treePC12.tsv
             |
             v
      (你贴的 pop-level 整理脚本) -> LM_summary_all.sorted.tsv / minP / PC1_only

#/mnt/spareHD_2/nu_287/q2_parallelism/08_q2_deltaAF_SNPlevel_vs_mtPCk_plus_treePC12.R
#snp level
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ========= paths =========
af_file  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
pc_file  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"
cov_file <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

out_dir  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_delta <- file.path(out_dir, "deltaAF_long.noAMO.tsv.gz")
out_mito  <- file.path(out_dir, "LM_perSNP_mitoPCk_plus_treePC12.tsv.gz")
out_tree  <- file.path(out_dir, "LM_perSNP_treeOnly.tsv.gz")
out_log   <- file.path(out_dir, "LOG_summary.txt")

# ========= params =========
min_depth <- 10
AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")
AK_marine <- "RS"
BC_marine <- "SAY"
min_n_AK <- 5
min_n_BC <- 8

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}
infer_region <- function(pop){
  fifelse(pop %in% c(AK_fresh, AK_marine), "AK",
          fifelse(pop %in% c(BC_fresh, BC_marine), "BC", NA_character_))
}

# ========= fast lm helpers =========
fast_lm_1x2t <- function(y, x, t1, t2){
  ok <- is.finite(y) & is.finite(x) & is.finite(t1) & is.finite(t2)
  y <- y[ok]; x <- x[ok]; t1 <- t1[ok]; t2 <- t2[ok]
  n <- length(y)
  if(n < 5) return(NULL)
  X <- cbind(1, x, t1, t2)
  fit <- lm.fit(X, y)
  p <- ncol(X); df_res <- n - p
  if(df_res <= 0) return(NULL)

  rss <- sum(fit$residuals^2)
  tss <- sum((y - mean(y))^2)
  R2 <- ifelse(tss > 0, 1 - rss/tss, NA_real_)
  adjR2 <- ifelse(tss > 0, 1 - (1 - R2) * (n - 1) / df_res, NA_real_)

  XtX_inv <- tryCatch(solve(crossprod(X)), error=function(e) NULL)
  if(is.null(XtX_inv)) return(NULL)
  sigma2 <- rss / df_res
  se <- sqrt(diag(XtX_inv) * sigma2)

  beta <- fit$coefficients
  tval <- beta / se
  pval <- 2 * pt(abs(tval), df=df_res, lower.tail=FALSE)

  list(beta=beta[2], se=se[2], t=tval[2], p=pval[2], R2=R2, adjR2=adjR2, n=n)
}

fast_lm_tree_only <- function(y, t1, t2){
  ok <- is.finite(y) & is.finite(t1) & is.finite(t2)
  y <- y[ok]; t1 <- t1[ok]; t2 <- t2[ok]
  n <- length(y)
  if(n < 4) return(NULL)
  X <- cbind(1, t1, t2)
  fit <- lm.fit(X, y)
  p <- ncol(X); df_res <- n - p
  if(df_res <= 0) return(NULL)

  rss <- sum(fit$residuals^2)
  tss <- sum((y - mean(y))^2)
  R2 <- ifelse(tss > 0, 1 - rss/tss, NA_real_)
  adjR2 <- ifelse(tss > 0, 1 - (1 - R2) * (n - 1) / df_res, NA_real_)
  list(R2=R2, adjR2=adjR2, n=n)
}

# ========= read AF =========
AF <- fread(cmd = paste("zcat", shQuote(af_file)), sep = ",", header = TRUE)
stopifnot(all(c("chr","pos","gene","pop","af","depth") %in% names(AF)))
AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= min_depth]
AF[, snp := paste(chr, pos, gene, sep=":")]

keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]

# drop AMO everywhere
AF <- AF[pop != "AMO"]

# add region now
AF[, region := infer_region(pop)]
AF <- AF[!is.na(region)]

# ========= read mitoPC + treePC =========
PC <- fread(pc_file)
PC[, pop := toupper(pop)]
PC <- PC[pop != "AMO"]

COV <- fread(cov_file)
COV[, pop := toupper(pop)]
COV <- COV[pop != "AMO"]

stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))

pc_cols <- grep("^mitoPC[1-5]$", names(PC), value=TRUE)
if(length(pc_cols) < 1){
  stop("pc_file must contain mitoPC1..mitoPC5 at least one. Found: ",
       paste(names(PC), collapse=","))
}

PC  <- unique(PC[, c("pop", pc_cols), with=FALSE])
COV <- unique(COV[, .(pop, treePC1, treePC2)])

PRED <- merge(PC, COV, by="pop", all=FALSE)

# ========= compute marine AF per SNP within region =========
marine_ak <- AF[region=="AK" & pop==AK_marine, .(region, snp, marine_af=af)]
marine_bc <- AF[region=="BC" & pop==BC_marine, .(region, snp, marine_af=af)]
marine_dt <- rbindlist(list(marine_ak, marine_bc), use.names=TRUE, fill=TRUE)

marine_dt <- marine_dt[, .(marine_af = mean(marine_af, na.rm=TRUE)), by=.(region, snp)]

# ========= freshwater deltaAF long =========
fresh_dt <- AF[(region=="AK" & pop %in% AK_fresh) | (region=="BC" & pop %in% BC_fresh),
               .(region, snp, pop, af)]

DT <- merge(fresh_dt, marine_dt, by=c("region","snp"), all.x=TRUE)
DT <- DT[is.finite(af) & is.finite(marine_af)]
DT[, deltaAF := af - marine_af]

# attach SNP annotation
SNP_ANN <- unique(AF[, .(snp, chr, pos, gene)])
DT <- merge(DT, SNP_ANN, by="snp", all.x=TRUE)

# attach predictors
DT <- merge(DT, PRED, by="pop", all=FALSE)

# ---- HARD GUARANTEE region exists ----
if(!("region" %in% names(DT))){
  # if merge created region.x/region.y, rescue it
  if("region.x" %in% names(DT)) setnames(DT, "region.x", "region")
  if("region.y" %in% names(DT) && !("region" %in% names(DT))) setnames(DT, "region.y", "region")
}
if(!("region" %in% names(DT))){
  DT[, region := infer_region(pop)]
}
DT <- DT[!is.na(region)]

# sanity print
cat("[DEBUG] names(DT):\n")
print(names(DT))

stopifnot(!any(duplicated(names(DT))))

# region-specific minimum n per SNP
DT[, n_in_snp := .N, by=.(region, snp)]
DT <- DT[(region=="AK" & n_in_snp >= min_n_AK) | (region=="BC" & n_in_snp >= min_n_BC)]
DT[, n_in_snp := NULL]

# ========= write deltaAF long (SAFE column selection) =========
cols_out <- c("region","snp","chr","pos","gene","pop","deltaAF","af","marine_af","treePC1","treePC2", pc_cols)
cols_out <- cols_out[cols_out %in% names(DT)]
fwrite(DT[, ..cols_out], out_delta, sep="\t", compress="gzip")

# ========= per-SNP regressions =========
tree_res <- DT[, {
  fit <- fast_lm_tree_only(deltaAF, treePC1, treePC2)
  if(is.null(fit)) return(NULL)
  .(chr=chr[1], pos=pos[1], gene=gene[1], R2=fit$R2, adjR2=fit$adjR2, n=fit$n)
}, by=.(region, snp)]
fwrite(tree_res, out_tree, sep="\t", compress="gzip")

mito_res_list <- vector("list", length(pc_cols))
names(mito_res_list) <- pc_cols
for(pc in pc_cols){
  tmp <- DT[, {
    fit <- fast_lm_1x2t(deltaAF, get(pc), treePC1, treePC2)
    if(is.null(fit)) return(NULL)
    .(chr=chr[1], pos=pos[1], gene=gene[1],
      mitoPC=pc,
      beta_mitoPC=fit$beta,
      se_mitoPC=fit$se,
      t_mitoPC=fit$t,
      p_mitoPC=fit$p,
      R2=fit$R2,
      adjR2=fit$adjR2,
      n=fit$n)
  }, by=.(region, snp)]
  mito_res_list[[pc]] <- tmp
}
mito_res <- rbindlist(mito_res_list, use.names=TRUE, fill=TRUE)
setorder(mito_res, region, mitoPC, p_mitoPC)
fwrite(mito_res, out_mito, sep="\t", compress="gzip")

# ========= log =========
log_txt <- c(
  sprintf("[INFO] AF rows after depth/pop filter: %d", nrow(AF)),
  sprintf("[INFO] deltaAF observations: %d", nrow(DT)),
  sprintf("[INFO] mitoPC cols used: %s", paste(pc_cols, collapse=",")),
  sprintf("[INFO] tree-only rows: %d", nrow(tree_res)),
  sprintf("[INFO] mito rows: %d", nrow(mito_res)),
  "",
  "=== populations used per region ===",
  capture.output(DT[, .N, by=.(region, pop)][order(region, -N)])
)
writeLines(log_txt, out_log)

cat("[OK] wrote:\n")
cat(" - ", out_delta, "\n", sep="")
cat(" - ", out_tree,  "\n", sep="")
cat(" - ", out_mito,  "\n", sep="")
cat(" - ", out_log,   "\n", sep="")

Rscript /mnt/spareHD_2/nu_287/q2_parallelism/08_q2_deltaAF_SNPlevel_vs_mtPCk_plus_treePC12.R





library(data.table)

f <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/LM_perSNP_mitoPCk_plus_treePC12.tsv.gz"
DT <- fread(cmd=paste("zcat", shQuote(f)))

# 只保留有p的行
DT <- DT[p_mitoPC != "" & !is.na(p_mitoPC)]
DT[, p := as.numeric(p_mitoPC)]
DT[, b := as.numeric(beta_mitoPC)]

# summary：每个 region × mitoPC
S <- DT[, .(
  n = .N,
  min_p = min(p, na.rm=TRUE),
  med_p = median(p, na.rm=TRUE),
  prop_p05 = mean(p < 0.05, na.rm=TRUE),
  prop_p10 = mean(p < 0.10, na.rm=TRUE),
  mean_beta = mean(b, na.rm=TRUE),
  prop_beta_pos = mean(b > 0, na.rm=TRUE)
), by = .(region, mitoPC)]

setorder(S, region, -prop_p05, med_p)
print(S)

# 给你每个 region 的“最强PC”（按 p<0.05 富集）
TOP <- S[, .SD[1], by=region]
cat("\n=== top PC per region (by prop_p05) ===\n")
print(TOP)




#pop level
cat 06_q2_deltaAF_meanSD_vs_mtPC_noAMO.R 
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =======================
# Inputs
# =======================
af_file  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"  # CSV gz: chr,pos,gene,pop,af,depth
pc_file  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"          # TSV: pop, mitoPC1..mitoPC5 (noAMO)
cov_file <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"                                      # TSV: pop, treePC1, treePC2, ...

out_dir  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_vs_mtPC_noAMO"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_pop_sum <- file.path(out_dir, "deltaAF_pop_summary_mean_sd.tsv")
out_lm_mean <- file.path(out_dir, "LM_meanDelta_vs_mitoPC1to5_plus_treePC12.tsv")
out_lm_sd   <- file.path(out_dir, "LM_sdDelta_vs_mitoPC1to5_plus_treePC12.tsv")
fig_mean_pc1 <- file.path(out_dir, "Fig_meanDelta_vs_mitoPC1.png")
fig_sd_pc1   <- file.path(out_dir, "Fig_sdDelta_vs_mitoPC1.png")

# =======================
# Params
# =======================
min_depth <- 10

AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")  # 9
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS",
               "BOOT","ECHO","LAW","GOS","ROB")               # 12  ✅ AMO removed
AK_marine <- "RS"
BC_marine <- "SAY"

# =======================
# Helpers
# =======================
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

stopifnot(file.exists(af_file), file.exists(pc_file), file.exists(cov_file))

# =======================
# Read AF (CSV.gz)
# =======================
AF <- fread(cmd = paste("zcat", shQuote(af_file)), sep = ",", header = TRUE)
req_af <- c("chr","pos","gene","pop","af","depth")
if(!all(req_af %in% names(AF))){
  stop("AF file missing required columns. Need: ", paste(req_af, collapse=", "),
       "\nFound: ", paste(names(AF), collapse=", "))
}
AF[, pop := normalize_pop(pop)]
AF[, gene := tolower(gene)]
AF <- AF[depth >= min_depth]

# only keep relevant pops (no AMO here)
keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]

AF[, snp := paste(chr, pos, gene, sep=":")]

# =======================
# Build marine lookup per SNP
# =======================
marine_RS  <- AF[pop == AK_marine, .(snp, af_marine = af)]
marine_SAY <- AF[pop == BC_marine, .(snp, af_marine = af)]

# =======================
# Compute per-pop ΔAF within each region
#   ΔAF(pop) = AF_pop - AF_marine(region)
# =======================
calc_delta <- function(AF_dt, fresh_pops, marine_dt, region_label){
  sub <- AF_dt[pop %in% fresh_pops, .(snp, pop, af)]
  m <- merge(sub, marine_dt, by="snp", all.x=TRUE)
  m[, delta := af - af_marine]
  # summarize per population
  out <- m[is.finite(delta), .(
    region = region_label,
    n_snps = .N,
    mean_delta = mean(delta, na.rm=TRUE),
    sd_delta   = sd(delta,   na.rm=TRUE)
  ), by=.(pop)]
  out[]
}

AK_pop <- calc_delta(AF, AK_fresh, marine_RS,  "AK")
BC_pop <- calc_delta(AF, BC_fresh, marine_SAY, "BC")
POP_SUM <- rbindlist(list(AK_pop, BC_pop), use.names=TRUE, fill=TRUE)

# =======================
# Read mitoPC + treePC
# =======================
PC <- fread(pc_file)
PC[, pop := toupper(pop)]
need_pcs <- paste0("mitoPC", 1:5)
if(!all(need_pcs %in% names(PC))){
  stop("pc_file must contain mitoPC1..mitoPC5. Found: ", paste(names(PC), collapse=", "))
}

COV <- fread(cov_file)
COV[, pop := toupper(pop)]
if(!all(c("pop","treePC1","treePC2") %in% names(COV))){
  stop("cov_file must contain pop, treePC1, treePC2 at minimum. Found: ", paste(names(COV), collapse=", "))
}

# merge
df <- merge(POP_SUM, PC[, c("pop", need_pcs), with=FALSE], by="pop", all=FALSE)
df <- merge(df, COV[, .(pop, treePC1, treePC2)], by="pop", all=FALSE)

# sanity: should be AK=9, BC=12 (unless some pops lack AF/PC/COV)
cat("[info] populations after merge:\n")
print(df[, .N, by=region])

fwrite(df, out_pop_sum, sep="\t")

# =======================
# LM runner
# =======================
fit_lm_grid <- function(dat, response_col){
  out <- list()
  for(reg in sort(unique(dat$region))){
    dd <- dat[region == reg]
    for(k in 1:5){
      pc <- paste0("mitoPC", k)
      fml <- as.formula(paste0(response_col, " ~ ", pc, " + treePC1 + treePC2"))
      m <- lm(fml, data=dd)

      sm <- summary(m)
      co <- coef(sm)

      if(!(pc %in% rownames(co))) next

      out[[length(out)+1]] <- data.table(
        region = reg,
        response = response_col,
        mitoPC = pc,
        beta_mitoPC = unname(co[pc, "Estimate"]),
        se_mitoPC   = unname(co[pc, "Std. Error"]),
        t_mitoPC    = unname(co[pc, "t value"]),
        p_mitoPC    = unname(co[pc, "Pr(>|t|)"]),
        adjR2       = unname(sm$adj.r.squared),
        R2          = unname(sm$r.squared),
        n           = nrow(dd)
      )
    }
  }
  rbindlist(out, use.names=TRUE, fill=TRUE)
}

LM_mean <- fit_lm_grid(df, "mean_delta")
LM_sd   <- fit_lm_grid(df, "sd_delta")

setorder(LM_mean, region, mitoPC)
setorder(LM_sd,   region, mitoPC)

fwrite(LM_mean, out_lm_mean, sep="\t")
fwrite(LM_sd,   out_lm_sd,   sep="\t")

# =======================
# Figures: PC1 only (simple, for quick viewing)
# =======================
p_mean <- ggplot(df, aes(x=mitoPC1, y=mean_delta)) +
  geom_point(size=2, alpha=0.85) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~region, scales="free_y") +
  labs(x="mitoPC1 (noAMO)", y="mean ΔAF (pop vs regional marine)", title="Mean ΔAF vs mitoPC1 (noAMO; +treePC1-2 in LM)") +
  theme_classic(base_size = 14)

p_sd <- ggplot(df, aes(x=mitoPC1, y=sd_delta)) +
  geom_point(size=2, alpha=0.85) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~region, scales="free_y") +
  labs(x="mitoPC1 (noAMO)", y="SD(ΔAF) (pop vs regional marine)", title="SD(ΔAF) vs mitoPC1 (noAMO; +treePC1-2 in LM)") +
  theme_classic(base_size = 14)

ggsave(fig_mean_pc1, p_mean, width=7.2, height=4.5, dpi=300)
ggsave(fig_sd_pc1,   p_sd,   width=7.2, height=4.5, dpi=300)

cat("\n[OK] wrote:\n")
cat(" - ", out_pop_sum, "\n", sep="")
cat(" - ", out_lm_mean, "\n", sep="")
cat(" - ", out_lm_sd, "\n", sep="")
cat(" - ", fig_mean_pc1, "\n", sep="")
cat(" - ", fig_sd_pc1, "\n", sep="")






#ΔAF ~ mtCluster + treePC + permutation（tree-based cluster）
MT_TREE=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_noDloop_noAMO.iqtree.treefile \
N_PERM=200 \
Rscript /mnt/spareHD_2/nu_287/q2_parallelism/09_mtCluster_deltaAF_SNPlevel_LM_perm.R


#delta/af  ～  cluster treepc
#!/usr/bin/env Rscript
# ============================================================
# 09_mtCluster_deltaAF_SNPlevel_LM_perm.R
#
# SNP-level regression (within region AK/BC):
#   deltaAF ~ mtCluster + treePC1 + treePC2
# Test mtCluster term by nested-model F-test (anova(reduced, full)).
#
# Permutation test:
#   permute mtCluster labels among POPs *within each region*,
#   recompute F for each SNP, get empirical p_perm.
#
# Inputs (noAMO):
#   DELTA (long table): deltaAF_long.noAMO.tsv.gz
#     must contain: region, snp, chr, pos, gene, pop, deltaAF, treePC1, treePC2
#   MT_TREE: mt_noDloop_noAMO*.treefile  (non-ultrametric OK)
#
# Outputs:
#   OUTDIR/
#     mtCluster_by_pop.noAMO.tsv
#     LM_perSNP_mtCluster_plus_treePC12.tsv.gz
#     LM_cluster_summary_by_region.tsv
#     Perm_perSNP_mtCluster_plus_treePC12.tsv.gz   (if N_PERM>0)
#     Fig_pHist_mtCluster_LM.png
#     Fig_QQ_mtCluster_LM.png
#     Fig_pHist_mtCluster_perm.png (if N_PERM>0)
#     Perm_global_summary.tsv      (global permutation summary)
#
# Run:
#   Rscript 09_mtCluster_deltaAF_SNPlevel_LM_perm.R
# or:
#   MT_TREE=/path/to/treefile DELTA=/path/to/delta OUTDIR=/path/to/out \
#   K=4 N_PERM=200 SEED=1 Rscript 09_mtCluster_deltaAF_SNPlevel_LM_perm.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
})

# ----------------------------
# Config (env override)
# ----------------------------
DELTA_FILE <- Sys.getenv("DELTA",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
)

OUT_DIR <- Sys.getenv("OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_noAMO"
)

K_CLUST <- as.integer(Sys.getenv("K", unset = "4"))
N_PERM  <- as.integer(Sys.getenv("N_PERM", unset = "200"))
SEED    <- as.integer(Sys.getenv("SEED", unset = "1"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

# ----------------------------
# Helpers
# ----------------------------
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

safe_zread <- function(f){
  if(!file.exists(f)) stop("File not found: ", f)
  fread(cmd = paste("zcat", shQuote(f)))
}

# Find a treefile if not explicitly provided
find_treefile <- function(){
  # If user provides MT_TREE, use it.
  tf <- Sys.getenv("MT_TREE", unset = "")
  if(nzchar(tf) && file.exists(tf)) return(tf)

  # Otherwise search common locations (no recursive arg; base Sys.glob doesn't support it)
  candidates <- unique(c(
    Sys.glob("/work/cyu/poolseq/PPalign_output/**/mt*noAMO*.treefile"),
    Sys.glob("/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/*noAMO*.treefile"),
    Sys.glob("/mnt/spareHD_2/nu_287/_assoc72_subunit/*noAMO*.treefile"),
    Sys.glob("/mnt/spareHD_2/mt_gene_tree/*noAMO*.treefile")
  ))

  # Sys.glob doesn't support "**" portably; the first pattern may return none.
  # If none found, do shallow glob on a few directories:
  if(length(candidates) == 0){
    candidates <- unique(c(
      Sys.glob("/work/cyu/poolseq/PPalign_output/*/*noAMO*.treefile"),
      Sys.glob("/work/cyu/poolseq/PPalign_output/*/*/*noAMO*.treefile"),
      Sys.glob("/work/cyu/poolseq/PPalign_output/*/*/*/*noAMO*.treefile"),
      Sys.glob("/work/cyu/poolseq/PPalign_output/*/*/*/*/*noAMO*.treefile"),
      Sys.glob("/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/*noAMO*.treefile"),
      Sys.glob("/mnt/spareHD_2/nu_287/_assoc72_subunit/*noAMO*.treefile"),
      Sys.glob("/mnt/spareHD_2/mt_gene_tree/*noAMO*.treefile")
    ))
  }

  candidates <- candidates[file.exists(candidates)]
  if(length(candidates) == 0) return(NA_character_)

  # Prefer iqtree treefile if available
  pref <- candidates[grepl("iqtree", candidates, ignore.case=TRUE)]
  if(length(pref) > 0) return(pref[1])
  candidates[1]
}

stop_msg_tree <- function(){
  cat("\n[ERROR] Cannot find mt noAMO treefile.\n")
  cat("Set it explicitly like:\n")
  cat("  MT_TREE=/full/path/to/mt_noDloop_noAMO.iqtree.treefile Rscript 09_...\n\n")
  stop("Tree file not found.")
}

# Compute F-test p-value for cluster term via nested model
lm_cluster_Ftest <- function(dt){
  dt <- dt[is.finite(deltaAF) & !is.na(cluster) & is.finite(treePC1) & is.finite(treePC2)]
  if(nrow(dt) < 6) return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))
  if(length(unique(dt$cluster)) < 2) return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))

  fit_full <- try(lm(deltaAF ~ cluster + treePC1 + treePC2, data=dt), silent=TRUE)
  fit_red  <- try(lm(deltaAF ~          treePC1 + treePC2, data=dt), silent=TRUE)
  if(inherits(fit_full, "try-error") || inherits(fit_red, "try-error")) {
    return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))
  }

  a <- try(anova(fit_red, fit_full), silent=TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2) {
    return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))
  }

  list(
    p   = as.numeric(a$`Pr(>F)`[2]),
    F   = as.numeric(a$F[2]),
    df1 = as.numeric(a$Df[2]),
    df2 = as.numeric(a$Res.Df[2]),
    n   = nrow(dt)
  )
}

# Permute cluster labels among POPs within a region
permute_clusters_within_region <- function(pop_cluster_dt, region_name){
  sub <- pop_cluster_dt[region == region_name]
  if(nrow(sub) == 0) return(data.table(pop=character(), cluster_perm=factor()))
  sub[, cluster_perm := sample(cluster)]
  sub[, .(pop, cluster_perm)]
}

# ----------------------------
# Load deltaAF long table
# ----------------------------
cat("[read] deltaAF_long: ", DELTA_FILE, "\n", sep="")
DEL <- safe_zread(DELTA_FILE)

req <- c("region","snp","chr","pos","gene","pop","deltaAF","treePC1","treePC2")
miss <- setdiff(req, names(DEL))
if(length(miss) > 0){
  cat("\n[ERROR] delta file missing columns:\n  ", paste(miss, collapse=", "), "\n", sep="")
  cat("\n[FOUND columns]\n")
  print(names(DEL))
  stop("delta file format mismatch.")
}

DEL[, pop := normalize_pop(pop)]
DEL <- DEL[pop != "AMO"]  # double safety

# Filter NA
DEL <- DEL[is.finite(deltaAF) & is.finite(treePC1) & is.finite(treePC2)]
cat("[info] rows after NA filter: ", nrow(DEL), "\n", sep="")

# ----------------------------
# Load mt tree -> mt clusters (non-ultrametric OK)
# ----------------------------
TREE_FILE <- find_treefile()
if(!is.character(TREE_FILE) || is.na(TREE_FILE) || !file.exists(TREE_FILE)) stop_msg_tree()

cat("[read tree] ", TREE_FILE, "\n", sep="")
tr <- read.tree(TREE_FILE)
cat("[tree] ntip = ", length(tr$tip.label), "\n", sep="")
cat("[tree] first tips: ", paste(head(tr$tip.label, 10), collapse=", "), "\n", sep="")

# Use patristic distances; no need ultrametric
D <- cophenetic(tr)  # matrix tip x tip
hc <- hclust(as.dist(D), method="average")
cl <- cutree(hc, k=K_CLUST)

mt_cluster <- data.table(
  pop = normalize_pop(names(cl)),
  mt_cluster = as.integer(cl)
)
mt_cluster <- unique(mt_cluster, by="pop")

# Keep only pops that appear in DEL (otherwise clusters for other tips unused)
pops_in_DEL <- unique(DEL$pop)
mt_cluster <- mt_cluster[pop %in% pops_in_DEL]

# Merge cluster
DEL <- merge(DEL, mt_cluster, by="pop", all.x=TRUE)

cat("\n[info] rows deltaAF: ", nrow(DEL), "\n", sep="")
cat("[info] unique pops: ", uniqueN(DEL$pop), "\n", sep="")
cat("[info] pops missing mt_cluster: ", sum(is.na(DEL$mt_cluster)), "\n", sep="")
if(sum(is.na(DEL$mt_cluster)) > 0){
  cat("[WARN] pops without mt_cluster (first 30):\n")
  print(head(unique(DEL[is.na(mt_cluster), pop]), 30))
}

# Drop missing cluster rows
DEL <- DEL[!is.na(mt_cluster)]
DEL[, cluster := factor(paste0("C", mt_cluster))]

# Save pop->cluster table (with region)
pop_region <- unique(DEL[, .(pop, region)], by="pop")
pop_cluster_out <- merge(pop_region, unique(DEL[, .(pop, mt_cluster)], by="pop"), by="pop", all.x=TRUE)
fwrite(pop_cluster_out[order(region, mt_cluster, pop)],
       file.path(OUT_DIR, "mtCluster_by_pop.noAMO.tsv"), sep="\t")
cat("[write] mtCluster_by_pop.noAMO.tsv\n")

# ----------------------------
# SNP-level LM by region
# ----------------------------
cat("\n[info] rows by region:\n")
print(DEL[, .N, by=region][order(region)])

setkey(DEL, region, snp)

regions <- sort(unique(DEL$region))
snps_by_region <- lapply(regions, function(r) unique(DEL[region==r, snp]))
names(snps_by_region) <- regions

cat("\n[info] #SNPs per region:\n")
for(r in regions){
  cat("  ", r, ": ", length(snps_by_region[[r]]), "\n", sep="")
}

lm_list <- vector("list", 1000)
idx <- 0L

for(r in regions){
  snps <- snps_by_region[[r]]
  for(s in snps){
    dt <- DEL[list(r, s)]
    chr_  <- dt$chr[1]
    pos_  <- dt$pos[1]
    gene_ <- dt$gene[1]

    out <- lm_cluster_Ftest(dt)

    cm <- dt[, .(mu=mean(deltaAF, na.rm=TRUE), n=.N), by=cluster][order(cluster)]
    cm_str <- paste0(cm$cluster, ":", sprintf("%.6g", cm$mu), "(n=", cm$n, ")", collapse=";")

    idx <- idx + 1L
    lm_list[[idx]] <- data.table(
      region=r, snp=s, chr=chr_, pos=pos_, gene=gene_,
      n=out$n, n_clusters=length(unique(dt$cluster)),
      F_cluster=out$F, df1=out$df1, df2=out$df2, p_cluster=out$p,
      cluster_means=cm_str
    )
  }
}

LM <- rbindlist(lm_list[seq_len(idx)], use.names=TRUE, fill=TRUE)

lm_out <- file.path(OUT_DIR, "LM_perSNP_mtCluster_plus_treePC12.tsv.gz")
fwrite(LM, lm_out, sep="\t")
cat("\n[write] ", lm_out, "\n", sep="")

SUM <- LM[, .(
  n_snps = .N,
  prop_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
  prop_p10 = mean(p_cluster < 0.10, na.rm=TRUE),
  min_p = suppressWarnings(min(p_cluster, na.rm=TRUE)),
  med_p = suppressWarnings(median(p_cluster, na.rm=TRUE))
), by=region][order(region)]

fwrite(SUM, file.path(OUT_DIR, "LM_cluster_summary_by_region.tsv"), sep="\t")
cat("[write] LM_cluster_summary_by_region.tsv\n\n")
cat("=== summary (cluster term) ===\n")
print(SUM)

# ----------------------------
# Permutation test
# ----------------------------
perm_out_file <- file.path(OUT_DIR, "Perm_perSNP_mtCluster_plus_treePC12.tsv.gz")
perm_global_file <- file.path(OUT_DIR, "Perm_global_summary.tsv")

if(N_PERM <= 0){
  cat("\n[perm] N_PERM=0 -> skip permutation.\n")
} else {

  cat("\n[perm] running permutations: N_PERM=", N_PERM, "\n", sep="")

  # pop -> region -> cluster table
  pop_cluster_dt <- unique(DEL[, .(pop, region, cluster)], by=c("pop","region"))
  # observed F
  OBS <- LM[, .(region, snp, F_obs=F_cluster)]
  setkey(OBS, region, snp)
  OBS[, ge_count := 0L]
  OBS[, n_perm_used := 0L]

  # Minimal X
  BASE <- DEL[, .(region, snp, pop, deltaAF, treePC1, treePC2)]
  setkey(BASE, region, snp)

  # progress
  for(b in seq_len(N_PERM)){
    if(b %% 25 == 0) cat("[perm] ", b, "/", N_PERM, "\n", sep="")

    perm_maps <- lapply(regions, function(r){
      pm <- permute_clusters_within_region(pop_cluster_dt, r)
      pm[, region := r]
      pm
    })
    PM <- rbindlist(perm_maps, use.names=TRUE, fill=TRUE)

    X <- merge(BASE, PM, by=c("region","pop"), all=FALSE)
    X[, cluster_perm := factor(cluster_perm)]
    setkey(X, region, snp)

    for(r in regions){
      snps <- snps_by_region[[r]]
      for(s in snps){
        dt <- X[list(r, s)]
        if(nrow(dt) < 6) next
        if(length(unique(dt$cluster_perm)) < 2) next

        dt2 <- dt
        dt2[, cluster := cluster_perm]
        out <- lm_cluster_Ftest(dt2)
        if(!is.finite(out$F)) next

        f_obs <- OBS[list(r, s), F_obs]
        if(!is.finite(f_obs)) next

        # update counts
        krow <- which(OBS$region==r & OBS$snp==s)
        if(length(krow)==1){
          OBS$n_perm_used[krow] <- OBS$n_perm_used[krow] + 1L
          if(out$F >= f_obs) OBS$ge_count[krow] <- OBS$ge_count[krow] + 1L
        }
      }
    }
  }

  OBS[, p_perm := (1 + ge_count) / (1 + n_perm_used)]

  PERM <- merge(LM, OBS[, .(region, snp, n_perm_used, p_perm)], by=c("region","snp"), all.x=TRUE)
  fwrite(PERM, perm_out_file, sep="\t")
  cat("\n[write] ", perm_out_file, "\n", sep="")

  # Global perm summary per region:
  G <- PERM[, .(
    n_snps=.N,
    prop_LM_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
    prop_perm_p05 = mean(p_perm < 0.05, na.rm=TRUE),
    med_perm_p = suppressWarnings(median(p_perm, na.rm=TRUE)),
    n_perm_used_med = suppressWarnings(median(n_perm_used, na.rm=TRUE))
  ), by=region][order(region)]
  fwrite(G, perm_global_file, sep="\t")
  cat("[write] Perm_global_summary.tsv\n")

  cat("\n=== permutation summary ===\n")
  print(G)
}

# ----------------------------
# Plots
# ----------------------------
p_hist <- ggplot(LM[is.finite(p_cluster)], aes(p_cluster)) +
  geom_histogram(bins=50) +
  facet_wrap(~region, ncol=1, scales="free_y") +
  labs(x="LM p-value for mtCluster term (F-test)", y="count",
       title="SNP-level association: ΔAF ~ mtCluster + treePC1 + treePC2") +
  theme_classic(base_size=14)
ggsave(file.path(OUT_DIR, "Fig_pHist_mtCluster_LM.png"), p_hist, width=7.2, height=6.5, dpi=300)

LMq <- LM[is.finite(p_cluster) & p_cluster > 0 & p_cluster <= 1]
LMq[, obs := -log10(sort(p_cluster))]
LMq[, exp := -log10(ppoints(.N))]
p_qq <- ggplot(LMq, aes(exp, obs)) +
  geom_point(alpha=0.4, size=0.8) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~region) +
  labs(x="Expected -log10(p)", y="Observed -log10(p)",
       title="QQ plot: mtCluster term p-values") +
  theme_classic(base_size=14)
ggsave(file.path(OUT_DIR, "Fig_QQ_mtCluster_LM.png"), p_qq, width=7.2, height=4.2, dpi=300)

if(N_PERM > 0 && file.exists(perm_out_file)){
  PERM <- safe_zread(perm_out_file)
  if("p_perm" %in% names(PERM)){
    p_hist2 <- ggplot(PERM[is.finite(p_perm)], aes(p_perm)) +
      geom_histogram(bins=50) +
      facet_wrap(~region, ncol=1, scales="free_y") +
      labs(x="Permutation p-value for mtCluster term", y="count",
           title=paste0("Permutation test (N_PERM=", N_PERM, ")")) +
      theme_classic(base_size=14)
    ggsave(file.path(OUT_DIR, "Fig_pHist_mtCluster_perm.png"), p_hist2, width=7.2, height=6.5, dpi=300)
  }
}

# ----------------------------
# Top hits
# ----------------------------
cat("\n=== Top 20 SNPs by LM p_cluster (per region) ===\n")
LM_top <- LM[is.finite(p_cluster)][order(p_cluster)][, head(.SD, 20), by=region]
print(LM_top[, .(region, snp, gene, chr, pos, n, n_clusters, F_cluster, p_cluster)])

cat("\n[OK] done. Outputs in:\n  ", OUT_DIR, "\n", sep="")





#ΔAF ~ manual mtCluster + treePC（perm0 + driver_cluster）
#perm0 delta af～ cluster+treepc
#!/usr/bin/env Rscript
# ============================================================
# /mnt/spareHD_2/nu_287/q2_parallelism/09_mtClusterManual_deltaAF_SNPlevel_perm0.R
# deltaAF ~ mtCluster(manual) + treePC1 + treePC2
# per SNP, within region
# NO permutation (perm=0)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

DELTA_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.withAMO.tsv"
OUTDIR       <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

safe_zread <- function(f){
  fread(cmd = paste("zcat", shQuote(f)))
}

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

lm_cluster_ftest <- function(dt){
  dt <- dt[
    is.finite(deltaAF) &
      is.finite(treePC1) &
      is.finite(treePC2) &
      !is.na(mtCluster)
  ]
  if (nrow(dt) < 6) return(list(p=NA_real_, F=NA_real_, n=nrow(dt)))
  if (length(unique(dt$mtCluster)) < 2) return(list(p=NA_real_, F=NA_real_, n=nrow(dt)))

  fit_full <- lm(deltaAF ~ mtCluster + treePC1 + treePC2, data = dt)
  fit_red  <- lm(deltaAF ~             treePC1 + treePC2, data = dt)
  a <- anova(fit_red, fit_full)

  list(
    p = as.numeric(a$`Pr(>F)`[2]),
    F = as.numeric(a$F[2]),
    n = nrow(dt)
  )
}

# ----------------------------
# Load data
# ----------------------------
cat("[read] deltaAF_long\n")
DEL <- safe_zread(DELTA_FILE)
DEL[, pop := normalize_pop(pop)]
DEL <- DEL[pop != "AMO"]

cat("[read] mtCluster_manual\n")
# fill=TRUE: allow blank lines / ragged lines
CL <- fread(CLUSTER_FILE, sep="\t", header=TRUE, fill=TRUE)

# drop truly empty rows
CL <- CL[!(is.na(pop) | pop=="")]

# enforce two columns
need <- c("pop","mtCluster")
miss <- setdiff(need, names(CL))
if(length(miss) > 0){
  stop("Cluster file missing columns: ", paste(miss, collapse=", "))
}

CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL, by="pop")

# merge
DEL <- merge(DEL, CL, by = "pop", all.x = TRUE)
DEL <- DEL[!is.na(mtCluster)]
DEL[, mtCluster := factor(mtCluster)]

cat("[info] pops per cluster:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# ----------------------------
# Per SNP LM
# ----------------------------
setkey(DEL, region, snp)
regions <- sort(unique(DEL$region))

res <- vector("list", 100000)
k <- 0L

for (r in regions) {
  snps <- unique(DEL[region == r, snp])
  for (s in snps) {
    dt <- DEL[list(r, s)]

    out <- lm_cluster_ftest(dt)

    cm <- dt[, .(
      mean_deltaAF = mean(deltaAF),
      n = .N
    ), by = mtCluster]

    cm[, abs_mean := abs(mean_deltaAF)]
    driver_i <- which.max(cm$abs_mean)
    driver <- if(length(driver_i)==1) as.character(cm$mtCluster[driver_i]) else NA_character_
    driver_mu <- if(length(driver_i)==1) cm$mean_deltaAF[driver_i] else NA_real_

    k <- k + 1L
    res[[k]] <- data.table(
      region = r,
      snp = s,
      gene = dt$gene[1],
      chr = dt$chr[1],
      pos = dt$pos[1],
      n = out$n,
      F_cluster = out$F,
      p_cluster = out$p,
      driver_cluster = driver,
      driver_mu = driver_mu
    )
  }
}

LM <- rbindlist(res[seq_len(k)])

fwrite(LM,
       file.path(OUTDIR, "LM_perSNP_mtCluster_manual_plus_treePC12.tsv.gz"),
       sep = "\t")

# ----------------------------
# Summary
# ----------------------------
SUM <- LM[, .(
  n_snps = .N,
  prop_p05 = mean(p_cluster < 0.05, na.rm = TRUE),
  prop_p10 = mean(p_cluster < 0.10, na.rm = TRUE),
  min_p = suppressWarnings(min(p_cluster, na.rm = TRUE)),
  med_p = suppressWarnings(median(p_cluster, na.rm = TRUE))
), by = region][order(region)]

fwrite(SUM,
       file.path(OUTDIR, "LM_cluster_summary_by_region.tsv"),
       sep = "\t")

cat("\n=== summary ===\n")
print(SUM)

# ----------------------------
# Which cluster drives signal?
# ----------------------------
SIG <- LM[is.finite(p_cluster) & p_cluster < 0.05]
DRIVER <- SIG[, .N, by = .(region, driver_cluster)]
DRIVER[, prop := N / sum(N), by = region]
setorder(DRIVER, region, -prop)

fwrite(DRIVER,
       file.path(OUTDIR, "LM_driverCluster_share_sigSNP.tsv"),
       sep = "\t")

cat("\n=== driver cluster share among significant SNPs (p<0.05) ===\n")
print(DRIVER)

# ----------------------------
# Plots
# ----------------------------
p1 <- ggplot(LM[is.finite(p_cluster)], aes(p_cluster)) +
  geom_histogram(bins = 50) +
  facet_wrap(~region, ncol = 1) +
  theme_classic(base_size = 14) +
  labs(x="p-value (mtCluster F-test)", y="count",
       title="deltaAF ~ mtCluster(manual) + treePC1 + treePC2")

ggsave(file.path(OUTDIR, "Fig_pHist_mtCluster_manual.png"),
       p1, width = 7, height = 6, dpi = 300)

p2 <- ggplot(DRIVER, aes(driver_cluster, prop)) +
  geom_col() +
  facet_wrap(~region, ncol = 1, scales="free_x") +
  theme_classic(base_size = 14) +
  labs(x="driver cluster (largest |mean deltaAF|)", y="share among significant SNPs",
       title="Which mtCluster drives deltaAF signal?")

ggsave(file.path(OUTDIR, "Fig_driverCluster_share_sigSNP.png"),
       p2, width = 8, height = 5, dpi = 300)

cat("\n[OK] done. Outputs in:\n", OUTDIR, "\n", sep="")




#bh
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# Gene enrichment from mtCluster-driven SNPs (BH after LM)
#
# Inputs (from your manual-cluster 09 script):
#   OUTDIR contains:
#     - LM_perSNP_mtCluster_plus_treePC12.tsv.gz
#   (optional) If LM file doesn't contain driver_cluster, we will
#   infer it from cluster means if present; but best is: LM already
#   has driver_cluster column from your 09 script output.
#
# Outputs:
#   - LM_with_q_and_driver.tsv.gz
#   - Gene_enrichment_overall.tsv.gz
#   - GeneCluster_enrichment.tsv.gz
#   - Fig_topGenes_byCluster_<region>.png
# ============================================================

# ----------------------------
# Config (env override)
# ----------------------------
IN_LM <- Sys.getenv("IN_LM",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0/LM_perSNP_mtCluster_plus_treePC12.tsv.gz"
)

OUTDIR <- Sys.getenv("OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0/gene_enrichment"
)

# BH threshold for calling SNPs significant
Q_CUT <- as.numeric(Sys.getenv("Q_CUT", unset="0.10"))   # try 0.05 later

# Minimum SNPs per gene in the "universe" to test (avoid tiny genes)
MIN_SNP_GENE <- as.integer(Sys.getenv("MIN_SNP_GENE", unset="20"))

# Only consider these regions (auto if blank)
REGIONS <- Sys.getenv("REGIONS", unset="")  # e.g. "AK,BC" or "" for all

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Read LM results
# ----------------------------
cat("[read] LM: ", IN_LM, "\n", sep="")
LM <- fread(IN_LM)

need <- c("region","snp","gene","p_cluster")
miss <- setdiff(need, names(LM))
if(length(miss)>0){
  cat("[ERROR] Missing required columns in LM: ", paste(miss, collapse=", "), "\n", sep="")
  cat("[FOUND] columns:\n")
  print(names(LM))
  stop("LM input mismatch.")
}

# region filter
if(nchar(REGIONS)>0){
  keep <- trimws(unlist(strsplit(REGIONS, ",")))
  LM <- LM[region %in% keep]
}

# BH per region
LM[, q_cluster := p.adjust(p_cluster, method="BH"), by=region]

# driver_cluster must exist for cluster-enrichment
# If your LM has driver_cluster column, use it.
# If not, try to infer from cluster_means string if present.
if(!("driver_cluster" %in% names(LM))){
  if("cluster_means" %in% names(LM)){
    # cluster_means like: "C1:0.12(n=3);C2:-0.05(n=2)..."
    # Infer driver as max |mean - overall mean| approx by max |mean|
    infer_driver <- function(x){
      if(is.na(x) || x=="") return(NA_character_)
      parts <- unlist(strsplit(x, ";", fixed=TRUE))
      # extract "Ck:mean"
      cl <- sub("^([^:]+):.*$", "\\1", parts)
      mu <- sub("^[^:]+:([^\\(]+)\\(.*$", "\\1", parts)
      mu <- suppressWarnings(as.numeric(mu))
      if(all(!is.finite(mu))) return(NA_character_)
      cl[ which.max(abs(mu)) ]
    }
    LM[, driver_cluster := vapply(cluster_means, infer_driver, character(1))]
    cat("[warn] driver_cluster not found; inferred from cluster_means.\n")
  } else {
    cat("[ERROR] driver_cluster not found and cannot infer (no cluster_means).\n")
    stop("Please rerun 09 so LM contains driver_cluster (recommended).")
  }
}

# Save LM with q + driver
out_lm2 <- file.path(OUTDIR, "LM_with_q_and_driver.tsv.gz")
fwrite(LM, out_lm2, sep="\t")
cat("[write] ", out_lm2, "\n", sep="")

# ----------------------------
# Define universe + significant sets
# ----------------------------
LM <- LM[is.finite(p_cluster) & is.finite(q_cluster)]
LM[, is_sig := (q_cluster < Q_CUT)]

cat("\n[info] SNP counts by region:\n")
print(LM[, .(
  n_snps=.N,
  n_sig=sum(is_sig, na.rm=TRUE),
  prop_sig=mean(is_sig, na.rm=TRUE),
  min_q=min(q_cluster, na.rm=TRUE),
  med_q=median(q_cluster, na.rm=TRUE)
), by=region][order(region)])

# ----------------------------
# 1) Overall gene enrichment (sig vs non-sig), ignoring cluster
# Fisher test per gene within region
# ----------------------------
cat("\n[run] gene overall enrichment (sig vs non-sig)\n")

# precompute region totals
REG_TOT <- LM[, .(
  N_total=.N,
  N_sig=sum(is_sig)
), by=region]

# per gene counts
G <- LM[, .(
  n_gene = .N,
  n_sig_gene = sum(is_sig)
), by=.(region, gene)]

# filter tiny genes
G <- G[n_gene >= MIN_SNP_GENE]

# Fisher per gene
fisher_one <- function(a, b, c, d){
  # matrix [[a,b],[c,d]]; return OR and p
  m <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
  ft <- fisher.test(m)
  list(or=unname(ft$estimate), p=ft$p.value)
}

G[, `:=`(or_overall=NA_real_, p_overall=NA_real_)]

for(r in unique(G$region)){
  tot <- REG_TOT[region==r]
  # a = sig in gene
  # b = non-sig in gene
  # c = sig not in gene
  # d = non-sig not in gene
  idx <- which(G$region==r)
  for(i in idx){
    a <- G$n_sig_gene[i]
    b <- G$n_gene[i] - a
    c <- tot$N_sig - a
    d <- (tot$N_total - tot$N_sig) - b
    out <- fisher_one(a,b,c,d)
    G$or_overall[i] <- out$or
    G$p_overall[i]  <- out$p
  }
}

G[, q_overall := p.adjust(p_overall, method="BH"), by=region]

out_g_overall <- file.path(OUTDIR, "Gene_enrichment_overall.tsv.gz")
fwrite(G[order(region, q_overall, p_overall)],
       out_g_overall, sep="\t")
cat("[write] ", out_g_overall, "\n", sep="")

# ----------------------------
# 2) Gene × driver_cluster enrichment among significant SNPs
# Question: within significant SNPs, is driver_cluster overrepresented for a gene?
#
# For each region, gene, cluster:
#   a = # sig SNPs in gene with driver_cluster==cl
#   b = # sig SNPs in gene with driver_cluster!=cl
#   c = # sig SNPs not in gene with driver_cluster==cl
#   d = # sig SNPs not in gene with driver_cluster!=cl
# Fisher exact test.
# ----------------------------
cat("\n[run] gene × driver_cluster enrichment (within significant SNPs)\n")

SIG <- LM[is_sig == TRUE & !is.na(driver_cluster)]

# if no sig SNPs, stop gracefully
if(nrow(SIG)==0){
  cat("[WARN] No significant SNPs at q<", Q_CUT, ". Try Q_CUT=0.10 or 0.20.\n", sep="")
  quit(save="no", status=0)
}

# totals among sig SNPs
SIG_TOT <- SIG[, .N, by=.(region)]        # total sig per region
SIG_CL  <- SIG[, .N, by=.(region, driver_cluster)]  # sig per region×cluster

# per gene×cluster counts among sig
GC <- SIG[, .N, by=.(region, gene, driver_cluster)]
setnames(GC, "N", "n_sig_gene_cl")

# add gene total sig
GS <- SIG[, .N, by=.(region, gene)]
setnames(GS, "N", "n_sig_gene")
GC <- merge(GC, GS, by=c("region","gene"), all.x=TRUE)

# add region totals and region×cluster totals
GC <- merge(GC, SIG_TOT, by="region", all.x=TRUE)
setnames(GC, "N", "n_sig_region")
GC <- merge(GC, SIG_CL, by=c("region","driver_cluster"), all.x=TRUE)
setnames(GC, "N", "n_sig_region_cl")

# add universe size per gene (all SNPs) to allow filtering by MIN_SNP_GENE
GC <- merge(GC, G[, .(region, gene, n_gene)], by=c("region","gene"), all.x=TRUE)

GC <- GC[!is.na(n_gene) & n_gene >= MIN_SNP_GENE]

# fisher for each row
GC[, `:=`(or_fisher=NA_real_, p_fisher=NA_real_)]

for(i in seq_len(nrow(GC))){
  a <- GC$n_sig_gene_cl[i]
  b <- GC$n_sig_gene[i] - a
  c <- GC$n_sig_region_cl[i] - a
  d <- (GC$n_sig_region[i] - GC$n_sig_region_cl[i]) - b
  out <- fisher_one(a,b,c,d)
  GC$or_fisher[i] <- out$or
  GC$p_fisher[i]  <- out$p
}

GC[, q_fisher := p.adjust(p_fisher, method="BH"), by=.(region, driver_cluster)]

out_gc <- file.path(OUTDIR, "GeneCluster_enrichment.tsv.gz")
fwrite(GC[order(region, driver_cluster, q_fisher, p_fisher)],
       out_gc, sep="\t")
cat("[write] ", out_gc, "\n", sep="")

# ----------------------------
# Plot: top enriched genes per cluster per region
# ----------------------------
plot_top <- function(dt, region_name){
  dt <- dt[region == region_name]
  if(nrow(dt)==0) return(NULL)

  # keep only strong signals for plotting (raw p for ranking, but show q)
  # You can tweak these:
  dt2 <- dt[is.finite(p_fisher) & is.finite(or_fisher)]
  if(nrow(dt2)==0) return(NULL)

  # pick top 10 genes per cluster by smallest q (then p)
  dt2 <- dt2[order(q_fisher, p_fisher)]
  top <- dt2[, head(.SD, 10), by=driver_cluster]

  top[, gene := factor(gene, levels=rev(unique(gene[order(driver_cluster, q_fisher, p_fisher)])))]

  p <- ggplot(top, aes(x=gene, y=log2(or_fisher))) +
    geom_col() +
    coord_flip() +
    facet_wrap(~driver_cluster, scales="free_y") +
    labs(
      title=paste0("Top gene enrichments by mt driver_cluster (", region_name,
                   ")  |  significant SNPs: q_cluster<", Q_CUT),
      x="gene",
      y="log2(Fisher OR) for cluster among significant SNPs"
    ) +
    theme_classic(base_size=13)

  p
}

for(r in unique(GC$region)){
  p <- plot_top(GC, r)
  if(!is.null(p)){
    fn <- file.path(OUTDIR, paste0("Fig_topGenes_byCluster_", r, ".png"))
    ggsave(fn, p, width=12, height=7, dpi=300)
    cat("[write] ", fn, "\n", sep="")
  }
}

# ----------------------------
# Quick console summaries
# ----------------------------
cat("\n=== Overall gene enrichment: how many genes at q<0.05? ===\n")
print(G[, .(
  n_genes=.N,
  n_q05=sum(q_overall < 0.05, na.rm=TRUE),
  n_q10=sum(q_overall < 0.10, na.rm=TRUE),
  best_q=min(q_overall, na.rm=TRUE)
), by=region][order(region)])

cat("\n=== Gene×cluster enrichment: how many tests at q<0.05 per cluster? ===\n")
print(GC[, .(
  n_tests=.N,
  n_q05=sum(q_fisher < 0.05, na.rm=TRUE),
  n_q10=sum(q_fisher < 0.10, na.rm=TRUE),
  best_q=min(q_fisher, na.rm=TRUE)
), by=.(region, driver_cluster)][order(region, driver_cluster)])

cat("\n[OK] done. Outputs in:\n  ", OUTDIR, "\n", sep="")


IN_LM=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0/LM_perSNP_mtCluster_manual_plus_treePC12.tsv.gz \
Rscript /mnt/spareHD_2/nu_287/q2_parallelism/10_gene_enrichment_from_cluster_LM.R




