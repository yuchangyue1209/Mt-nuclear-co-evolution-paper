#delta af～delta mt mt sync to af file
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# --------------------
# Paths (env override)
# --------------------
SYNC <- Sys.getenv("SYNC", unset="/work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.sync")
POPS <- Sys.getenv("POPS", unset="/work/cyu/poolseq/PPalign_output/mtDNA_bam/mt_pops_in_sync_order.short.txt")
OUTD <- Sys.getenv("OUTD", unset="/work/cyu/poolseq/PPalign_output/mtDNA_bam")
dir.create(OUTD, showWarnings=FALSE, recursive=TRUE)

OUT_GZ <- Sys.getenv("OUT", unset=file.path(OUTD, "mtAF_long.major.pipegz.tsv.gz"))

# --------------------
# Params
# --------------------
MIN_COV  <- as.integer(Sys.getenv("MIN_COV", unset="10"))
COV_MODE <- Sys.getenv("COV_MODE", unset="ATCG")  # "ATCG" or "ALL"
CHUNK    <- as.integer(Sys.getenv("CHUNK", unset="20000"))

# --------------------
# Helpers
# --------------------
if(!file.exists(POPS)) stop("POPS file not found: ", POPS)
if(!file.exists(SYNC)) stop("SYNC file not found: ", SYNC)

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

parse_sync_cell <- function(cell){
  if(is.na(cell) || cell=="" || cell==".") return(c(NA,NA,NA,NA,NA,NA))
  parts <- strsplit(cell, ":", fixed=TRUE)[[1]]
  if(length(parts) < 6) return(c(NA,NA,NA,NA,NA,NA))
  suppressWarnings(as.numeric(parts[1:6]))
}

major_minor <- function(A, T, C, G){
  v <- c(A=A, T=T, C=C, G=G)
  if(any(!is.finite(v))) return(list(majorBase=NA_character_, majorCount=NA_real_,
                                     minorBase=NA_character_, minorCount=NA_real_))
  o <- order(v, decreasing=TRUE, na.last=TRUE)
  list(
    majorBase  = names(v)[o[1]],
    majorCount = v[o[1]],
    minorBase  = names(v)[o[2]],
    minorCount = v[o[2]]
  )
}

open_sync_con <- function(path){
  if(grepl("\\.gz$", path, ignore.case=TRUE)) gzfile(path, open="rt") else file(path, open="rt")
}

# --------------------
# Load pop names (must match sync column order)
# --------------------
pops <- fread(POPS, header=FALSE)[[1]]
pops <- normalize_pop(pops)

cat("[info] #pops from POPS:", length(pops), "\n")
cat("[info] first pops:", paste(head(pops, 10), collapse=", "), "\n")

# --------------------
# Open gzip output pipe (guaranteed valid .gz)
# This overwrites OUT_GZ without needing rm
# --------------------
cmd <- sprintf("gzip -c > %s", shQuote(OUT_GZ))
out_con <- pipe(cmd, open="wb")
on.exit(close(out_con), add=TRUE)

# write header (tab-separated)
header <- paste(
  c("chr","pos","ref","pop","cov","A","T","C","G","N","del","majorBase","majorAF","minorBase","minorAF"),
  collapse="\t"
)
writeLines(header, con=out_con, sep="\n")

# --------------------
# Stream read sync file, write in chunks
# --------------------
sync_con <- open_sync_con(SYNC)
on.exit(close(sync_con), add=TRUE)

line_buf <- character(CHUNK)
total_lines <- 0L
written <- 0L

repeat {
  n <- 0L
  while(n < CHUNK){
    x <- readLines(sync_con, n=1, warn=FALSE)
    if(length(x)==0) break
    n <- n + 1L
    line_buf[n] <- x
  }
  if(n == 0L) break

  lines <- line_buf[seq_len(n)]
  total_lines <- total_lines + n

  spl <- strsplit(lines, "\t", fixed=TRUE)
  ncol_expected <- 3L + length(pops)

  keep_i <- which(vapply(spl, length, integer(1)) == ncol_expected)
  if(length(keep_i)==0) next
  spl <- spl[keep_i]

  chr <- vapply(spl, `[[`, character(1), 1)
  pos <- as.integer(vapply(spl, `[[`, character(1), 2))
  ref <- vapply(spl, `[[`, character(1), 3)

  sample_cells <- lapply(spl, function(v) v[4:length(v)])
  M <- do.call(rbind, sample_cells)
  colnames(M) <- pops

  DTm <- as.data.table(M)
  DTm[, `:=`(chr=chr, pos=pos, ref=ref)]

  long <- melt(DTm,
               id.vars=c("chr","pos","ref"),
               variable.name="pop",
               value.name="cell",
               variable.factor=FALSE)

  counts <- t(vapply(long$cell, parse_sync_cell, numeric(6)))
  colnames(counts) <- c("A","T","C","G","N","del")
  long[, c("A","T","C","G","N","del") := as.data.table(counts)]
  long[, cell := NULL]

  if(COV_MODE == "ALL"){
    long[, cov := A+T+C+G+N+del]
  } else {
    long[, cov := A+T+C+G]
  }

  long <- long[is.finite(cov) & cov >= MIN_COV]
  if(nrow(long)==0) next

  mm <- mapply(major_minor, long$A, long$T, long$C, long$G, SIMPLIFY=FALSE)
  long[, majorBase := vapply(mm, `[[`, character(1), "majorBase")]
  long[, majorCount := vapply(mm, `[[`, numeric(1), "majorCount")]
  long[, minorBase := vapply(mm, `[[`, character(1), "minorBase")]
  long[, minorCount := vapply(mm, `[[`, numeric(1), "minorCount")]

  long[, majorAF := majorCount / cov]
  long[, minorAF := minorCount / cov]
  long[, c("majorCount","minorCount") := NULL]

  setcolorder(long, c("chr","pos","ref","pop","cov","A","T","C","G","N","del",
                      "majorBase","majorAF","minorBase","minorAF"))

  # write chunk WITHOUT header; use write.table to connection
  write.table(long,
              file=out_con,
              sep="\t",
              row.names=FALSE,
              col.names=FALSE,
              quote=FALSE)

  written <- written + nrow(long)

  if(total_lines %% (CHUNK*5) == 0){
    cat("[info] processed sync lines:", total_lines, "| wrote long rows:", written, "\n")
  }
}

cat("\n[OK] done.\n")
cat("SYNC:", SYNC, "\n")
cat("POPS:", POPS, "\n")
cat("OUT :", OUT_GZ, "\n")
cat("Processed sync lines:", total_lines, "\n")
cat("Wrote long rows     :", written, "\n")



zcat /work/cyu/poolseq/PPalign_output/mtDNA_bam/mtAF_long.major.pipegz.tsv.gz \
| awk -F'\t' 'NR==1 || !($2>=15653 && $2<=16543)' \
| gzip > /work/cyu/poolseq/PPalign_output/mtDNA_bam/mtAF_long.major.noDloop.tsv.gz



#mtΔAF（noDloop）
#/work/cyu/poolseq/PPalign_output/mtDNA_bam/02_mtDeltaAF_long.noDloop.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 02_mtDeltaAF_long.noDloop.R
# Build mtDeltaAF long table from mtAF_long.major.noDloop.tsv.gz
#
# Region-specific marine baseline:
#   AK: RS
#   BC: SAY
#
# Input long AF must contain:
#   chr, pos, ref, pop, cov, and either:
#     - af   (generic)
#     - majorAF (recommended for your file)
#
# Outputs:
#   mtDeltaAF_long.noDloop.tsv.gz
#   mtDeltaAF_pop_summary.noDloop.tsv
# ============================================================

MT_AF <- Sys.getenv("MT_AF",
  unset="/work/cyu/poolseq/PPalign_output/mtDNA_bam/mtAF_long.major.noDloop.tsv.gz"
)
OUTDIR <- Sys.getenv("OUTDIR",
  unset="/work/cyu/poolseq/PPalign_output/mtDNA_bam"
)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_DELTA <- file.path(OUTDIR, "mtDeltaAF_long.noDloop.tsv.gz")
OUT_SUM   <- file.path(OUTDIR, "mtDeltaAF_pop_summary.noDloop.tsv")

MIN_COV <- as.numeric(Sys.getenv("MIN_COV", unset="30"))

AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")
AK_marine <- "RS"
BC_marine <- "SAY"

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}
infer_region <- function(pop){
  fifelse(pop %in% c(AK_fresh, AK_marine), "AK",
          fifelse(pop %in% c(BC_fresh, BC_marine), "BC", NA_character_))
}

stopifnot(file.exists(MT_AF))

cat("[read] ", MT_AF, "\n", sep="")
AF <- fread(cmd = paste("zcat", shQuote(MT_AF)))

# ---- choose AF column ----
af_col <- NA_character_
if("af" %in% names(AF)) {
  af_col <- "af"
} else if("majorAF" %in% names(AF)) {
  af_col <- "majorAF"
} else if("minorAF" %in% names(AF)) {
  af_col <- "minorAF"
}
if(!nzchar(af_col) || is.na(af_col)){
  cat("[ERROR] cannot find AF column. Need 'af' or 'majorAF' (or 'minorAF').\n")
  cat("[FOUND] columns:\n"); print(names(AF))
  stop("mtAF_long input mismatch.")
}
cat("[info] using AF column: ", af_col, "\n", sep="")

need <- c("chr","pos","ref","pop","cov")
miss <- setdiff(need, names(AF))
if(length(miss) > 0){
  cat("[ERROR] input missing required columns: ", paste(miss, collapse=", "), "\n", sep="")
  cat("[FOUND] columns:\n"); print(names(AF))
  stop("mtAF_long input mismatch.")
}

AF[, pop := normalize_pop(pop)]
AF[, region := infer_region(pop)]
AF <- AF[!is.na(region)]

# keep only pops of interest
keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]

# numeric AF
AF[, AF_use := suppressWarnings(as.numeric(get(af_col)))]

AF <- AF[is.finite(AF_use) & is.finite(cov) & cov >= MIN_COV]
AF[, site := paste(chr, pos, sep=":")]

cat("[info] rows after cov+AF filter: ", nrow(AF), "\n", sep="")
cat("[info] pops kept: ", uniqueN(AF$pop), "\n", sep="")

# -----------------------------
# Marine baseline per site per region
# -----------------------------
marine_ak <- AF[region=="AK" & pop==AK_marine, .(region, site, marine_af=AF_use, marine_cov=cov)]
marine_bc <- AF[region=="BC" & pop==BC_marine, .(region, site, marine_af=AF_use, marine_cov=cov)]
MAR <- rbindlist(list(marine_ak, marine_bc), use.names=TRUE, fill=TRUE)
MAR <- MAR[, .(
  marine_af  = mean(marine_af,  na.rm=TRUE),
  marine_cov = mean(marine_cov, na.rm=TRUE)
), by=.(region, site)]

# -----------------------------
# Freshwater rows + delta
# -----------------------------
FRESH <- AF[
  (region=="AK" & pop %in% AK_fresh) |
  (region=="BC" & pop %in% BC_fresh),
  .(region, site, chr, pos, ref, pop, af=AF_use, cov)
]

DT <- merge(FRESH, MAR, by=c("region","site"), all.x=TRUE)
DT <- DT[is.finite(marine_af)]
DT[, deltaAF_mt := af - marine_af]

setorder(DT, region, chr, pos, pop)

# -----------------------------
# Write outputs (gzip safely)
# -----------------------------
fwrite(
  DT[, .(region, chr, pos, ref, site, pop,
         af, cov, marine_af, marine_cov, deltaAF_mt,
         af_col_used = af_col)],
  file = OUT_DELTA,
  sep  = "\t",
  compress = "gzip"
)

SUM <- DT[is.finite(deltaAF_mt), .(
  n_sites = .N,
  mean_deltaAF = mean(deltaAF_mt, na.rm=TRUE),
  sd_deltaAF   = sd(deltaAF_mt,   na.rm=TRUE),
  med_abs_deltaAF = median(abs(deltaAF_mt), na.rm=TRUE),
  mean_cov = mean(cov, na.rm=TRUE),
  mean_marine_cov = mean(marine_cov, na.rm=TRUE)
), by=.(region, pop)][order(region, pop)]

fwrite(SUM, OUT_SUM, sep="\t")

cat("\n[OK] wrote:\n")
cat("  - ", OUT_DELTA, "\n", sep="")
cat("  - ", OUT_SUM, "\n", sep="")
cat("\n[check] head summary:\n")
print(head(SUM, 10))


#03
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =======================
# Inputs
# =======================
NU_DELTA <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
MT_DELTA <- "/work/cyu/poolseq/PPalign_output/mtDNA_bam/mtDeltaAF_long.noDloop.tsv.gz"
COV_FILE <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"  # pop treePC1 treePC2 ...

# =======================
# Outputs
# =======================
OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_popLevel.noAMO.noDloop"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_POP  <- file.path(OUTDIR, "popLevel_nuDeltaAF_vs_mtDeltaAF.tsv.gz")
OUT_LM   <- file.path(OUTDIR, "LM_popLevel_nuDeltaAF_vs_mtDeltaAF_plus_treePC12.tsv")
FIG_MEAN <- file.path(OUTDIR, "Fig_meanNuDeltaAF_vs_meanMtDeltaAF.png")
FIG_SD   <- file.path(OUTDIR, "Fig_sdNuDeltaAF_vs_sdMtDeltaAF.png")

# =======================
# Helpers
# =======================
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

zread <- function(f){
  if(!file.exists(f)) stop("File not found: ", f)
  fread(cmd = paste("zcat", shQuote(f)))
}

# =======================
# Read covariates
# =======================
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))
COV <- unique(COV[, .(pop, treePC1, treePC2)])

# =======================
# Read nu ΔAF long
# =======================
cat("[read] nu deltaAF: ", NU_DELTA, "\n", sep="")
NU <- zread(NU_DELTA)
need_nu <- c("region","pop","deltaAF","treePC1","treePC2")
miss <- setdiff(need_nu, names(NU))
if(length(miss)>0) stop("NU missing cols: ", paste(miss, collapse=", "))

NU[, pop := normalize_pop(pop)]
NU <- NU[is.finite(deltaAF)]

NU_POP <- NU[, .(
  n_sites_nu = .N,
  mean_nuDeltaAF = mean(deltaAF, na.rm=TRUE),
  sd_nuDeltaAF   = sd(deltaAF,   na.rm=TRUE),
  med_abs_nuDeltaAF = median(abs(deltaAF), na.rm=TRUE)
), by=.(region, pop)]

# =======================
# Read mt ΔAF long
# =======================
cat("[read] mt deltaAF: ", MT_DELTA, "\n", sep="")
MT <- zread(MT_DELTA)
# your mt file columns include: region pop deltaAF_mt af cov marine_af marine_cov ...
need_mt <- c("region","pop","deltaAF_mt","cov","marine_cov")
miss <- setdiff(need_mt, names(MT))
if(length(miss)>0) stop("MT missing cols: ", paste(miss, collapse=", "))

MT[, pop := normalize_pop(pop)]
MT <- MT[is.finite(deltaAF_mt)]

MT_POP <- MT[, .(
  n_sites_mt = .N,
  mean_mtDeltaAF = mean(deltaAF_mt, na.rm=TRUE),
  sd_mtDeltaAF   = sd(deltaAF_mt,   na.rm=TRUE),
  med_abs_mtDeltaAF = median(abs(deltaAF_mt), na.rm=TRUE),
  mean_cov_mt = mean(cov, na.rm=TRUE),
  mean_marine_cov = mean(marine_cov, na.rm=TRUE)
), by=.(region, pop)]

# =======================
# Merge pop-level tables
# =======================
POP <- merge(NU_POP, MT_POP, by=c("region","pop"), all=FALSE)
POP <- merge(POP, COV, by="pop", all=FALSE)

setorder(POP, region, pop)

fwrite(POP, OUT_POP, sep="\t", compress="gzip")
cat("[write] ", OUT_POP, "\n", sep="")

cat("\n[info] pops kept after merge:\n")
print(POP[, .N, by=region][order(region)])

# =======================
# Pop-level LM (within each region)
# =======================
lm_one <- function(dd, y, x){
  fml <- as.formula(paste0(y, " ~ ", x, " + treePC1 + treePC2"))
  m <- lm(fml, data=dd)
  sm <- summary(m)
  co <- coef(sm)
  if(!(x %in% rownames(co))) return(NULL)
  data.table(
    response=y, predictor=x,
    beta = unname(co[x,"Estimate"]),
    se   = unname(co[x,"Std. Error"]),
    t    = unname(co[x,"t value"]),
    p    = unname(co[x,"Pr(>|t|)"]),
    R2   = unname(sm$r.squared),
    adjR2= unname(sm$adj.r.squared),
    n    = nrow(dd)
  )
}

LM <- rbindlist(lapply(sort(unique(POP$region)), function(r){
  dd <- POP[region==r]
  out1 <- lm_one(dd, "mean_nuDeltaAF", "mean_mtDeltaAF")
  out2 <- lm_one(dd, "sd_nuDeltaAF",   "sd_mtDeltaAF")
  out  <- rbindlist(list(out1, out2), use.names=TRUE, fill=TRUE)
  if(nrow(out)==0) return(NULL)
  out[, region := r]
  out
}), use.names=TRUE, fill=TRUE)

fwrite(LM, OUT_LM, sep="\t")
cat("[write] ", OUT_LM, "\n", sep="")

# =======================
# Plots
# =======================
p_mean <- ggplot(POP, aes(x=mean_mtDeltaAF, y=mean_nuDeltaAF)) +
  geom_point(size=2, alpha=0.85) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~region, scales="free") +
  theme_classic(base_size=14) +
  labs(x="mean mtΔAF (pop vs marine baseline; noDloop)",
       y="mean nuΔAF (pop vs marine baseline; noAMO)",
       title="Pop-level coupling: mean nuΔAF vs mean mtΔAF (+treePC1-2)")

ggsave(FIG_MEAN, p_mean, width=7.6, height=4.6, dpi=300)
cat("[write] ", FIG_MEAN, "\n", sep="")

p_sd <- ggplot(POP, aes(x=sd_mtDeltaAF, y=sd_nuDeltaAF)) +
  geom_point(size=2, alpha=0.85) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~region, scales="free") +
  theme_classic(base_size=14) +
  labs(x="sd mtΔAF (noDloop)",
       y="sd nuΔAF (noAMO)",
       title="Pop-level coupling: sd nuΔAF vs sd mtΔAF (+treePC1-2)")

ggsave(FIG_SD, p_sd, width=7.6, height=4.6, dpi=300)
cat("[write] ", FIG_SD, "\n", sep="")

cat("\n[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")


#04_nuDeltaAF_SNPlevel_vs_mtDeltaAF_popSummary.R
#!/usr/bin/env Rscript
# ============================================================
# 04_nuDeltaAF_SNPlevel_vs_mtDeltaAF_popSummary.R
#
# Goal (per region AK/BC):
#   For each nuclear SNP, regress nu ΔAF across populations on a
#   population-level summary of mt ΔAF (computed from mtDeltaAF_long).
#
# You asked to run ALL of these (per SNP, within each region):
#   1) nu ΔAF            ~ mean_mtΔAF       + treePC1 + treePC2
#   2) nu ΔAF            ~ mean_abs_mtΔAF   + treePC1 + treePC2
#   3) nu ΔAF            ~ sd_mtΔAF         + treePC1 + treePC2
#   4) |nu ΔAF| (MAG)    ~ mean_mtΔAF       + treePC1 + treePC2
#   5) |nu ΔAF| (MAG)    ~ mean_abs_mtΔAF   + treePC1 + treePC2
#   6) |nu ΔAF| (MAG)    ~ sd_mtΔAF         + treePC1 + treePC2
#
# Inputs:
#   NU_DELTA: deltaAF_long.noAMO.tsv.gz  (nuclear SNP-level deltaAF long)
#   MT_DELTA: mtDeltaAF_long.noDloop.tsv.gz (mt site-level deltaAF long)
#   COV_FILE: covariates.treePC.tsv (fallback if NU_DELTA lacks treePC)
#
# Outputs (OUTDIR):
#   popLevel_mtDeltaAF_summary.tsv.gz
#   snpLevel_table.noAMO.noDloop.tsv.gz
#   LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.tsv.gz
#   LOG_summary.txt
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ----------------------------
# Inputs
# ----------------------------
NU_DELTA <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
MT_DELTA <- "/work/cyu/poolseq/PPalign_output/mtDNA_bam/mtDeltaAF_long.noDloop.tsv.gz"
COV_FILE <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

# ----------------------------
# Outputs
# ----------------------------
OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.noAMO.noDloop"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_MTPOP <- file.path(OUTDIR, "popLevel_mtDeltaAF_summary.tsv.gz")
OUT_DT    <- file.path(OUTDIR, "snpLevel_table.noAMO.noDloop.tsv.gz")
OUT_LM    <- file.path(OUTDIR, "LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.tsv.gz")
OUT_LOG   <- file.path(OUTDIR, "LOG_summary.txt")

# ----------------------------
# Params
# ----------------------------
MIN_POPS_PER_SNP_AK <- 5
MIN_POPS_PER_SNP_BC <- 8

# ----------------------------
# Helpers
# ----------------------------
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

zread <- function(f){
  if(!file.exists(f)) stop("File not found: ", f)
  fread(cmd = paste("zcat", shQuote(f)))
}

# fast lm for y ~ x + t1 + t2
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

# rescue treePC1/2 after merges (handle .x/.y)
rescue_treePC <- function(DT){
  rescue_one <- function(base, xcol, ycol){
    if(base %in% names(DT)) return(invisible(NULL))
    if(xcol %in% names(DT) && ycol %in% names(DT)){
      DT[, (base) := fifelse(is.finite(get(xcol)), get(xcol), get(ycol))]
      return(invisible(NULL))
    }
    if(xcol %in% names(DT)){ setnames(DT, xcol, base); return(invisible(NULL)) }
    if(ycol %in% names(DT)){ setnames(DT, ycol, base); return(invisible(NULL)) }
    stop("Cannot find ", base, " in DT. Columns: ", paste(names(DT), collapse=", "))
  }
  rescue_one("treePC1", "treePC1.x", "treePC1.y")
  rescue_one("treePC2", "treePC2.x", "treePC2.y")

  extra <- intersect(c("treePC1.x","treePC1.y","treePC2.x","treePC2.y"), names(DT))
  if(length(extra) > 0) DT[, (extra) := NULL]
  DT
}

# ----------------------------
# 1) mt pop-level summary
# ----------------------------
cat("[read] mt deltaAF: ", MT_DELTA, "\n", sep="")
MT <- zread(MT_DELTA)

need_mt <- c("region","pop","deltaAF_mt","cov","marine_cov")
miss <- setdiff(need_mt, names(MT))
if(length(miss) > 0){
  cat("[ERROR] MT missing cols: ", paste(miss, collapse=", "), "\n", sep="")
  cat("[FOUND] MT cols:\n"); print(names(MT))
  stop("MT input mismatch.")
}

MT[, pop := normalize_pop(pop)]
MT <- MT[is.finite(deltaAF_mt) & is.finite(cov)]

MT_POP <- MT[, .(
  n_sites_mt = .N,
  mean_mtDeltaAF      = mean(deltaAF_mt, na.rm=TRUE),
  mean_abs_mtDeltaAF  = mean(abs(deltaAF_mt), na.rm=TRUE),
  sd_mtDeltaAF        = sd(deltaAF_mt, na.rm=TRUE),
  med_abs_mtDeltaAF   = median(abs(deltaAF_mt), na.rm=TRUE),
  mean_cov_mt         = mean(cov, na.rm=TRUE),
  mean_marine_cov     = mean(marine_cov, na.rm=TRUE)
), by=.(region, pop)]

fwrite(MT_POP, OUT_MTPOP, sep="\t", compress="gzip")
cat("[write] ", OUT_MTPOP, "\n", sep="")

# ----------------------------
# 2) nuclear deltaAF long + covariates
# ----------------------------
cat("[read] nu deltaAF: ", NU_DELTA, "\n", sep="")
NU <- zread(NU_DELTA)

need_nu_core <- c("region","snp","chr","pos","gene","pop","deltaAF")
miss <- setdiff(need_nu_core, names(NU))
if(length(miss) > 0){
  cat("[ERROR] NU missing cols: ", paste(miss, collapse=", "), "\n", sep="")
  cat("[FOUND] NU cols:\n"); print(names(NU))
  stop("NU input mismatch.")
}

NU[, pop := normalize_pop(pop)]
NU <- NU[is.finite(deltaAF)]

# COV fallback (only used if NU lacks treePC columns)
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))
COV <- unique(COV[, .(pop, treePC1, treePC2)])

# Merge: NU (per SNP x pop) + MT_POP (per pop) + COV (treePC)
DT <- merge(NU, MT_POP, by=c("region","pop"), all=FALSE)
DT <- merge(DT, COV, by="pop", all.x=TRUE)  # all.x=TRUE to keep NU rows; rescue will fill

DT <- rescue_treePC(DT)
stopifnot(all(c("treePC1","treePC2") %in% names(DT)))

# drop rows missing predictors
DT <- DT[
  is.finite(deltaAF) &
  is.finite(mean_mtDeltaAF) &
  is.finite(mean_abs_mtDeltaAF) &
  is.finite(sd_mtDeltaAF) &
  is.finite(treePC1) & is.finite(treePC2)
]

# enforce min pops per SNP within region
DT[, n_in_snp := .N, by=.(region, snp)]
DT <- DT[(region=="AK" & n_in_snp >= MIN_POPS_PER_SNP_AK) | (region=="BC" & n_in_snp >= MIN_POPS_PER_SNP_BC)]
DT[, n_in_snp := NULL]

# Save merged SNP-level table (one row per region×snp×pop)
fwrite(DT, OUT_DT, sep="\t", compress="gzip")
cat("[write] ", OUT_DT, "\n", sep="")

cat("\n[info] pops kept after merge:\n")
print(unique(DT[, .N, by=.(region, pop)])[order(region, -N)])

# ----------------------------
# 3) per-SNP regressions
# ----------------------------
cat("[run] per-SNP LMs (within region)\n")

# Model grid you requested
model_grid <- data.table(
  response_name = c("nuDeltaAF", "nuDeltaAF", "nuDeltaAF",
                    "abs_nuDeltaAF", "abs_nuDeltaAF", "abs_nuDeltaAF"),
  y_col         = c("deltaAF",   "deltaAF",   "deltaAF",
                    "abs_deltaAF","abs_deltaAF","abs_deltaAF"),
  x_name        = c("mean_mtDeltaAF", "mean_abs_mtDeltaAF", "sd_mtDeltaAF",
                    "mean_mtDeltaAF", "mean_abs_mtDeltaAF", "sd_mtDeltaAF"),
  x_col         = c("mean_mtDeltaAF", "mean_abs_mtDeltaAF", "sd_mtDeltaAF",
                    "mean_mtDeltaAF", "mean_abs_mtDeltaAF", "sd_mtDeltaAF")
)

DT[, abs_deltaAF := abs(deltaAF)]

# run
regions <- sort(unique(DT$region))
out_list <- vector("list", 1000L)
k <- 0L

setkey(DT, region, snp)

for(r in regions){
  snps <- unique(DT[region==r, snp])
  for(s in snps){
    dd <- DT[list(r, s)]
    if(nrow(dd) < 5) next

    for(i in seq_len(nrow(model_grid))){
      y <- dd[[ model_grid$y_col[i] ]]
      x <- dd[[ model_grid$x_col[i] ]]
      t1 <- dd$treePC1
      t2 <- dd$treePC2

      fit <- fast_lm_1x2t(y, x, t1, t2)
      if(is.null(fit)) next

      k <- k + 1L
      out_list[[k]] <- data.table(
        region = r,
        snp    = s,
        chr    = dd$chr[1],
        pos    = dd$pos[1],
        gene   = dd$gene[1],
        response = model_grid$response_name[i],
        predictor = model_grid$x_name[i],
        beta = fit$beta,
        se   = fit$se,
        t    = fit$t,
        p    = fit$p,
        R2   = fit$R2,
        adjR2= fit$adjR2,
        n    = fit$n
      )
    }
  }
}

LM <- rbindlist(out_list[seq_len(k)], use.names=TRUE, fill=TRUE)
setorder(LM, region, response, predictor, p)

fwrite(LM, OUT_LM, sep="\t", compress="gzip")
cat("[write] ", OUT_LM, "\n", sep="")

# ----------------------------
# 4) Log summary
# ----------------------------
log_txt <- c(
  sprintf("[INFO] NU_DELTA: %s", NU_DELTA),
  sprintf("[INFO] MT_DELTA: %s", MT_DELTA),
  sprintf("[INFO] COV_FILE: %s", COV_FILE),
  "",
  sprintf("[INFO] merged DT rows: %d", nrow(DT)),
  sprintf("[INFO] unique pops: %d", uniqueN(DT$pop)),
  sprintf("[INFO] unique SNPs: %d", uniqueN(DT$snp)),
  sprintf("[INFO] LM rows: %d", nrow(LM)),
  "",
  "=== pops per region (after filters) ===",
  capture.output(DT[, .N, by=.(region, pop)][order(region, pop)]),
  "",
  "=== quick LM summary: prop p<0.05 by region×response×predictor ===",
  capture.output(
    LM[, .(
      n_tests=.N,
      prop_p05=mean(p < 0.05, na.rm=TRUE),
      min_p=min(p, na.rm=TRUE),
      med_p=median(p, na.rm=TRUE)
    ), by=.(region, response, predictor)][order(region, response, predictor)]
  )
)
writeLines(log_txt, OUT_LOG)
cat("[write] ", OUT_LOG, "\n", sep="")

cat("\n[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")


