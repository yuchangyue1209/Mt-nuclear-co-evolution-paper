#MT NUCLEAR E
#M124 mtpc1+ all tree pc 
cyu@stickleback:/mnt/spareHD_2/nu_287$ cat run_hill_assoc72_noPheno_withMETA_csv.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# INPUTS (edit paths)
# =========================
AF72   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
COV    <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
META   <- "/mnt/spareHD_2/nu_287/meta_pop_env.csv"   # <-- your CSV here

# optional Q2 for bridge
Q2_PAR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_perSNP.tsv.gz"

OUTPFX <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/HILL72_noPheno"

N_TREEPC <- 4
DROP_RECENT <- FALSE   # TRUE: keep only Marine/Freshwater for Habitat-based models
SEED <- 1

pick_first <- function(x, candidates){
  hit <- candidates[candidates %in% x]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}

# =========================
# load
# =========================
message("[1/5] read AF72: ", AF72)
af   <- fread(AF72)

message("[2/5] read COV: ", COV)
cov  <- fread(COV)

message("[2/5] read META (CSV): ", META)
meta <- fread(META)   # fread auto-detects CSV/TSV; works for .csv

# ---- AF columns ----
col_pop  <- pick_first(names(af), c("pop","Population","population","Pop","POP"))
col_gene <- pick_first(names(af), c("gene","Gene"))
col_chr  <- pick_first(names(af), c("chr","CHR","chrom","Chrom","chromosome"))
col_pos  <- pick_first(names(af), c("pos","POS","position"))
col_af   <- pick_first(names(af), c("af","AF","nuAF","allele_freq","freq","p"))

if(any(is.na(c(col_pop,col_gene,col_pos,col_af))))
  stop("AF72 missing required cols. Need pop/gene/pos/AF.")

setnames(af, col_pop, "pop")
setnames(af, col_gene, "gene")
setnames(af, col_pos, "pos")
setnames(af, col_af, "AF")
if(!is.na(col_chr)) setnames(af, col_chr, "chr")
if(!("chr" %in% names(af))) af[, chr := "NA"]

# ---- COV columns ----
col_pop_cov <- pick_first(names(cov), c("pop","Population","population","Pop","POP"))
col_mito1   <- pick_first(names(cov), c("mitoPC1","mtPC1"))
col_mito2   <- pick_first(names(cov), c("mitoPC2","mtPC2"))

if(any(is.na(c(col_pop_cov,col_mito1))))
  stop("COV missing pop/mitoPC1. Expect columns like pop, mitoPC1, treePC1..")

setnames(cov, col_pop_cov, "pop")
setnames(cov, col_mito1, "mitoPC1")
if(!is.na(col_mito2)) setnames(cov, col_mito2, "mitoPC2")

treepcs <- paste0("treePC", 1:N_TREEPC)
treepcs <- treepcs[treepcs %in% names(cov)]
if(length(treepcs) < 2) stop("Too few treePC columns in COV.")

# ---- META columns (env) ----
m_pop <- pick_first(names(meta), c("pop","Population","population","Pop","POP"))
m_reg <- pick_first(names(meta), c("Region","region"))
m_lat <- pick_first(names(meta), c("Latitude","lat","latitude"))
m_lon <- pick_first(names(meta), c("Longitude","lon","longitude"))
m_hab <- pick_first(names(meta), c("Habitat","habitat"))

if(any(is.na(c(m_pop,m_reg,m_lat,m_lon,m_hab))))
  stop("META CSV missing required env cols. Need Population/Region/Latitude/Longitude/Habitat (or pop/...).")

setnames(meta, m_pop, "pop")
setnames(meta, m_reg, "Region")
setnames(meta, m_lat, "Latitude")
setnames(meta, m_lon, "Longitude")
setnames(meta, m_hab, "Habitat")

meta <- unique(meta[, .(pop, Region, Latitude, Longitude, Habitat)])

# =========================
# merge AF + COV + META
# =========================
message("[3/5] merge AF + COV + META")
cov2 <- cov[, c("pop","mitoPC1","mitoPC2",treepcs), with=FALSE]
dt <- merge(af, cov2,  by="pop", all.x=TRUE)
dt <- merge(dt, meta, by="pop", all.x=TRUE)

dt <- dt[!is.na(AF) & !is.na(mitoPC1)]
dt <- dt[!is.na(Habitat) & !is.na(Region) & !is.na(Latitude) & !is.na(Longitude)]

if(DROP_RECENT){
  dt <- dt[Habitat %in% c("Marine","Freshwater")]
}

dt[, Habitat := factor(Habitat)]
dt[, Region  := factor(Region)]

# SNP id
dt[, snp_id := paste(chr, pos, gene, sep=":")]

# =========================
# formulas
# =========================
rhs_tree <- paste(treepcs, collapse=" + ")

# Hill pieces (no phenotype):
# M1: Gn~Gmt | demography
# M2: Gn~Gmt + E | demography   (variance partition / "分账")
# M3: Gn~Gmt*E + ...           (Gmt×E interaction)
f_M1 <- as.formula(paste0("AF ~ mitoPC1 + ", rhs_tree))
f_M2 <- as.formula(paste0("AF ~ mitoPC1 + Habitat + Region + Latitude + Longitude + ", rhs_tree))
f_M3 <- as.formula(paste0("AF ~ mitoPC1 * Habitat + Region + Latitude + Longitude + ", rhs_tree))

get_coef_p <- function(fit, term){
  s <- summary(fit)$coefficients
  if(!(term %in% rownames(s))) return(c(NA_real_, NA_real_))
  c(s[term,"Estimate"], s[term,"Pr(>|t|)"])
}

# =========================
# per-SNP regression
# =========================
message("[4/5] per-SNP regression (M1/M2/M3)")
set.seed(SEED)

snps <- unique(dt$snp_id)
res_list <- vector("list", length(snps))

for(i in seq_along(snps)){
  sid <- snps[i]
  d <- dt[snp_id == sid]

  if(nrow(d) < (2 + length(treepcs))) next

  fit1 <- try(lm(f_M1, data=d), silent=TRUE)
  if(inherits(fit1,"try-error")) next
  b1 <- get_coef_p(fit1, "mitoPC1")

  fit2 <- try(lm(f_M2, data=d), silent=TRUE)
  b2 <- if(!inherits(fit2,"try-error")) get_coef_p(fit2, "mitoPC1") else c(NA,NA)

  fit3 <- try(lm(f_M3, data=d), silent=TRUE)
  if(!inherits(fit3,"try-error")){
    b3 <- get_coef_p(fit3, "mitoPC1")
    sm <- summary(fit3)$coefficients
    int_rows <- grep("^mitoPC1:Habitat", rownames(sm), value=TRUE)
    if(length(int_rows)>0){
      trm <- int_rows[1]
      bi  <- c(sm[trm,"Estimate"], sm[trm,"Pr(>|t|)"])
    } else {
      trm <- NA_character_; bi <- c(NA,NA)
    }
  } else {
    b3 <- c(NA,NA); trm <- NA_character_; bi <- c(NA,NA)
  }

  res_list[[i]] <- data.table(
    snp_id=sid,
    chr=d$chr[1], pos=d$pos[1], gene=d$gene[1],
    n_pop=nrow(d),
    beta_mitoPC1_M1=b1[1], p_mitoPC1_M1=b1[2],
    beta_mitoPC1_M2=b2[1], p_mitoPC1_M2=b2[2],
    beta_mitoPC1_M3=b3[1], p_mitoPC1_M3=b3[2],
    int_term=trm,
    beta_mitoPC1xHab=bi[1], p_mitoPC1xHab=bi[2]
  )

  if(i %% 2000 == 0) message("  ... ", i, "/", length(snps))
}

res <- rbindlist(res_list, fill=TRUE)
if(nrow(res)==0) stop("No results produced. Check merges / pop labels.")

# BH
res[, BH_mitoPC1_M1 := p.adjust(p_mitoPC1_M1, method="BH")]
res[, BH_mitoPC1_M2 := p.adjust(p_mitoPC1_M2, method="BH")]
res[, BH_mitoPC1_M3 := p.adjust(p_mitoPC1_M3, method="BH")]
res[, BH_mitoPC1xHab:= p.adjust(p_mitoPC1xHab, method="BH")]

out_snp  <- paste0(OUTPFX, "_perSNP_models_withEnv.tsv.gz")
fwrite(res, out_snp, sep="\t")
message("[ok] wrote: ", out_snp)

# gene summary
gene_sum <- res[, .(
  n_snps=.N,
  n_hit_M1=sum(BH_mitoPC1_M1<0.05, na.rm=TRUE),
  prop_hit_M1=mean(BH_mitoPC1_M1<0.05, na.rm=TRUE),
  n_hit_M2=sum(BH_mitoPC1_M2<0.05, na.rm=TRUE),
  prop_hit_M2=mean(BH_mitoPC1_M2<0.05, na.rm=TRUE),
  n_hit_int=sum(BH_mitoPC1xHab<0.05, na.rm=TRUE),
  prop_hit_int=mean(BH_mitoPC1xHab<0.05, na.rm=TRUE)
), by=gene][order(-prop_hit_M1, -n_hit_M1, gene)]

out_gene <- paste0(OUTPFX, "_gene_summary_withEnv.tsv")
fwrite(gene_sum, out_gene, sep="\t")
message("[ok] wrote: ", out_gene)

# optional bridge
if(file.exists(Q2_PAR)){
  q2 <- fread(Q2_PAR)
  q2_chr  <- pick_first(names(q2), c("chr","CHR","chrom","Chrom"))
  q2_pos  <- pick_first(names(q2), c("pos","POS","position"))
  q2_gene <- pick_first(names(q2), c("gene","Gene"))
  q2_y    <- pick_first(names(q2), c("is_parallel","isParallel","parallel","parallel_174","parallel_gold15"))
  if(!any(is.na(c(q2_chr,q2_pos,q2_gene,q2_y)))){
    setnames(q2, q2_chr, "chr"); setnames(q2, q2_pos, "pos"); setnames(q2, q2_gene, "gene"); setnames(q2, q2_y, "is_parallel")
    q2[, snp_id := paste(chr, pos, gene, sep=":")]

    br <- merge(q2[, .(snp_id, is_parallel)],
                res[, .(snp_id, gene, p_mitoPC1_M1)],
                by="snp_id", all.x=TRUE)
    br <- br[!is.na(is_parallel) & !is.na(p_mitoPC1_M1)]
    br[, mt_strength := -log10(p_mitoPC1_M1)]
    fit <- glm(is_parallel ~ mt_strength + gene, data=br, family=binomial())
    sm  <- summary(fit)$coefficients
    out_bridge <- data.table(term=rownames(sm),
                             beta=sm[,"Estimate"],
                             se=sm[,"Std. Error"],
                             z=sm[,"z value"],
                             p=sm[,"Pr(>|z|)"])
    out_b <- paste0(OUTPFX, "_bridge_parallelism_glm.tsv")
    fwrite(out_bridge, out_b, sep="\t")
    message("[ok] wrote: ", out_b)
  } else {
    message("[warn] Q2_PAR present but required columns not found; skipped bridge.")
  }
}

message("[DONE] Hill (no phenotype) models with META CSV completed.")








等价raw q1   
#mt12 tree12 m123
#/mnt/spareHD_2/nu_287/run_hill_assoc72_main_PC12_tree12_M123.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ==========================================================
# Main analysis (Hill): mtPC1+mtPC2 + treePC1+treePC2
#   M1: AF ~ mitoPC1 + mitoPC2 + treePC1 + treePC2
#   M2: + Habitat + Region + Latitude + Longitude
#   M3: + mitoPC1:Habitat interaction (PC2 as covariate)
# Output:
#   per-SNP table (with BH/FDR)
#   per-gene summary table
# ==========================================================

AF72 <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
COV  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
META <- "/mnt/spareHD_2/nu_287/meta_pop_env.csv"

OUTD <- "/mnt/spareHD_2/nu_287/_assoc72_subunit"
OUTPFX <- file.path(OUTD, "HILL72_MAIN_mt12_tree12_M123")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

SEED <- 1
DROP_RECENT <- FALSE  # set TRUE if you want only Marine+Freshwater

# -------------------------
# helpers
# -------------------------
pick_first <- function(x, candidates){
  hit <- candidates[candidates %in% x]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}

get_coef_p <- function(fit, term){
  s <- summary(fit)$coefficients
  if(!(term %in% rownames(s))) return(c(NA_real_, NA_real_))
  c(unname(s[term,"Estimate"]), unname(s[term,"Pr(>|t|)"]))
}

# ==========================================================
# [1/6] read inputs
# ==========================================================
message("[1/6] read AF72: ", AF72)
af <- fread(AF72)

message("[2/6] read COV: ", COV)
cov <- fread(COV)

message("[3/6] read META CSV: ", META)
meta <- fread(META)

# ==========================================================
# normalize AF columns
# ==========================================================
af_pop  <- pick_first(names(af), c("pop","Population","Pop","population"))
af_gene <- pick_first(names(af), c("gene","Gene"))
af_chr  <- pick_first(names(af), c("chr","CHR","chrom","Chrom","chromosome"))
af_pos  <- pick_first(names(af), c("pos","POS","position"))
af_af   <- pick_first(names(af), c("AF","af","nuAF","allele_freq","freq","p"))

if(any(is.na(c(af_pop, af_gene, af_pos, af_af)))){
  stop("AF72 missing required columns: need pop/gene/pos/AF")
}

setnames(af, af_pop,  "pop")
setnames(af, af_gene, "gene")
setnames(af, af_pos,  "pos")
setnames(af, af_af,   "AF")
if(!is.na(af_chr)) setnames(af, af_chr, "chr")
if(!("chr" %in% names(af))) af[, chr := "NA"]

# normalize pop labels: "10_THE" -> "THE"
af[, pop := trimws(as.character(pop))]
af[, pop := sub("^[0-9]+_", "", pop)]

# ==========================================================
# normalize COV columns
# ==========================================================
cov_pop <- pick_first(names(cov), c("pop","Population","Pop","population"))
if(is.na(cov_pop)) stop("COV missing pop column")
setnames(cov, cov_pop, "pop")
cov[, pop := trimws(as.character(pop))]

# mtPC1/2
m1 <- pick_first(names(cov), c("mitoPC1","mtPC1"))
m2 <- pick_first(names(cov), c("mitoPC2","mtPC2"))
if(is.na(m1) || is.na(m2)) stop("COV missing mitoPC1/mitoPC2 columns")
setnames(cov, m1, "mitoPC1")
setnames(cov, m2, "mitoPC2")

# treePC1/2 (force)
need_tree <- c("treePC1","treePC2")
if(!all(need_tree %in% names(cov))){
  stop("COV missing treePC1/treePC2")
}

# ==========================================================
# normalize META columns
# ==========================================================
m_pop <- pick_first(names(meta), c("Population","pop","Pop","population"))
m_reg <- pick_first(names(meta), c("Region","region"))
m_lat <- pick_first(names(meta), c("Latitude","lat","latitude"))
m_lon <- pick_first(names(meta), c("Longitude","lon","longitude"))
m_hab <- pick_first(names(meta), c("Habitat","habitat"))

if(any(is.na(c(m_pop,m_reg,m_lat,m_lon,m_hab)))){
  stop("META missing required cols: Population/Region/Latitude/Longitude/Habitat")
}

setnames(meta, m_pop, "Population")
setnames(meta, m_reg, "Region")
setnames(meta, m_lat, "Latitude")
setnames(meta, m_lon, "Longitude")
setnames(meta, m_hab, "Habitat")

meta[, Population := trimws(as.character(Population))]
meta[, Region     := trimws(as.character(Region))]
meta[, Habitat    := trimws(as.character(Habitat))]

meta_env <- unique(meta[, .(Population, Region, Latitude, Longitude, Habitat)])

# ==========================================================
# [4/6] merge AF + COV + META
# ==========================================================
message("[4/6] merge AF + COV + META")
cov_keep <- cov[, .(pop, mitoPC1, mitoPC2, treePC1, treePC2)]

dt <- merge(af, cov_keep, by="pop", all.x=TRUE)
dt <- merge(dt, meta_env, by.x="pop", by.y="Population", all.x=TRUE)

message(sprintf("[diag] nrow(dt)=%d  n_pop=%d  n_snp=%d",
                nrow(dt), length(unique(dt$pop)),
                length(unique(paste(dt$chr, dt$pos, dt$gene, sep=":")))))
message(sprintf("[diag] NA mitoPC1=%d  NA mitoPC2=%d  NA Habitat=%d  NA Region=%d",
                sum(is.na(dt$mitoPC1)), sum(is.na(dt$mitoPC2)),
                sum(is.na(dt$Habitat)), sum(is.na(dt$Region))))

# hard filters: must have AF + PCs + env
dt <- dt[!is.na(AF) & !is.na(mitoPC1) & !is.na(mitoPC2)]
dt <- dt[!is.na(treePC1) & !is.na(treePC2)]
dt <- dt[!is.na(Habitat) & !is.na(Region) & !is.na(Latitude) & !is.na(Longitude)]

if(DROP_RECENT){
  dt <- dt[Habitat %in% c("Marine","Freshwater")]
}

dt[, Habitat := factor(Habitat)]
dt[, Region  := factor(Region)]
dt[, snp_id := paste(chr, pos, gene, sep=":")]

# ==========================================================
# define models
# ==========================================================
f_M1 <- AF ~ mitoPC1 + mitoPC2 + treePC1 + treePC2
f_M2 <- AF ~ mitoPC1 + mitoPC2 + Habitat + Region + Latitude + Longitude + treePC1 + treePC2
# interaction: only mitoPC1 * Habitat, mitoPC2 as covariate
f_M3 <- AF ~ mitoPC1*Habitat + mitoPC2 + Region + Latitude + Longitude + treePC1 + treePC2

# ==========================================================
# [5/6] per-SNP regression
# ==========================================================
message("[5/6] per-SNP regression (M1/M2/M3; mtPC1+mtPC2; treePC1+treePC2)")
set.seed(SEED)

snps <- unique(dt$snp_id)
res_list <- vector("list", length(snps))

for(i in seq_along(snps)){
  sid <- snps[i]
  d <- dt[snp_id == sid]

  # quick sanity: need enough rows
  if(nrow(d) < 10) next

  fit1 <- try(lm(f_M1, data=d), silent=TRUE)
  if(inherits(fit1, "try-error")) next
  b1_pc1 <- get_coef_p(fit1, "mitoPC1")
  b1_pc2 <- get_coef_p(fit1, "mitoPC2")

  fit2 <- try(lm(f_M2, data=d), silent=TRUE)
  if(!inherits(fit2, "try-error")){
    b2_pc1 <- get_coef_p(fit2, "mitoPC1")
    b2_pc2 <- get_coef_p(fit2, "mitoPC2")
  } else {
    b2_pc1 <- c(NA,NA); b2_pc2 <- c(NA,NA)
  }

  fit3 <- try(lm(f_M3, data=d), silent=TRUE)
  if(!inherits(fit3, "try-error")){
    b3_pc1 <- get_coef_p(fit3, "mitoPC1")
    b3_pc2 <- get_coef_p(fit3, "mitoPC2")

    sm <- summary(fit3)$coefficients
    int_rows <- grep("^mitoPC1:Habitat", rownames(sm), value=TRUE)
    if(length(int_rows)>0){
      trm <- int_rows[1]
      bi  <- c(unname(sm[trm,"Estimate"]), unname(sm[trm,"Pr(>|t|)"]))
    } else {
      trm <- NA_character_; bi <- c(NA,NA)
    }
  } else {
    b3_pc1 <- c(NA,NA); b3_pc2 <- c(NA,NA); trm <- NA_character_; bi <- c(NA,NA)
  }

  res_list[[i]] <- data.table(
    snp_id=sid,
    chr=d$chr[1], pos=d$pos[1], gene=d$gene[1],
    n_pop=nrow(d),

    beta_pc1_M1=b1_pc1[1], p_pc1_M1=b1_pc1[2],
    beta_pc2_M1=b1_pc2[1], p_pc2_M1=b1_pc2[2],

    beta_pc1_M2=b2_pc1[1], p_pc1_M2=b2_pc1[2],
    beta_pc2_M2=b2_pc2[1], p_pc2_M2=b2_pc2[2],

    beta_pc1_M3=b3_pc1[1], p_pc1_M3=b3_pc1[2],
    beta_pc2_M3=b3_pc2[1], p_pc2_M3=b3_pc2[2],

    int_term=trm,
    beta_pc1xHab=bi[1], p_pc1xHab=bi[2]
  )

  if(i %% 2000 == 0) message("  ... ", i, "/", length(snps))
}

res <- rbindlist(res_list, fill=TRUE)
if(nrow(res)==0) stop("No results produced. Check merges/pop labels/filters.")

# BH/FDR
res[, BH_pc1_M1 := p.adjust(p_pc1_M1, method="BH")]
res[, BH_pc2_M1 := p.adjust(p_pc2_M1, method="BH")]
res[, BH_pc1_M2 := p.adjust(p_pc1_M2, method="BH")]
res[, BH_pc2_M2 := p.adjust(p_pc2_M2, method="BH")]
res[, BH_pc1_M3 := p.adjust(p_pc1_M3, method="BH")]
res[, BH_pc2_M3 := p.adjust(p_pc2_M3, method="BH")]
res[, BH_pc1xHab:= p.adjust(p_pc1xHab, method="BH")]

out_snp <- paste0(OUTPFX, "_perSNP.tsv.gz")
fwrite(res, out_snp, sep="\t")
message("[ok] wrote: ", out_snp)

# ==========================================================
# [6/6] gene-level summary
# ==========================================================
gene_sum <- res[, .(
  n_snps=.N,

  n_hit_pc1_M1=sum(BH_pc1_M1 < 0.05, na.rm=TRUE),
  prop_hit_pc1_M1=mean(BH_pc1_M1 < 0.05, na.rm=TRUE),
  n_hit_pc2_M1=sum(BH_pc2_M1 < 0.05, na.rm=TRUE),
  prop_hit_pc2_M1=mean(BH_pc2_M1 < 0.05, na.rm=TRUE),

  n_hit_pc1_M2=sum(BH_pc1_M2 < 0.05, na.rm=TRUE),
  prop_hit_pc1_M2=mean(BH_pc1_M2 < 0.05, na.rm=TRUE),
  n_hit_pc2_M2=sum(BH_pc2_M2 < 0.05, na.rm=TRUE),
  prop_hit_pc2_M2=mean(BH_pc2_M2 < 0.05, na.rm=TRUE),

  n_hit_pc1xHab=sum(BH_pc1xHab < 0.05, na.rm=TRUE),
  prop_hit_pc1xHab=mean(BH_pc1xHab < 0.05, na.rm=TRUE)
), by=gene][order(-prop_hit_pc1_M1, -n_hit_pc1_M1, gene)]

out_gene <- paste0(OUTPFX, "_gene_summary.tsv")
fwrite(gene_sum, out_gene, sep="\t")
message("[ok] wrote: ", out_gene)

message("[DONE] Main Hill analysis finished.")














#nu~ e meta
/mnt/spareHD_2/nu_287/_assoc72_subunit/run_hill72_env_E123.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

AF72 <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
COV  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
META <- "/mnt/spareHD_2/nu_287/meta_pop_env.csv"

OUTD   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit"
OUTPFX <- file.path(OUTD, "HILL72_ENVONLY_E123")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

pick_first <- function(x, candidates){
  hit <- candidates[candidates %in% x]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}
get_coef_p <- function(fit, term){
  s <- summary(fit)$coefficients
  if(!(term %in% rownames(s))) return(c(NA_real_, NA_real_))
  c(unname(s[term,"Estimate"]), unname(s[term,"Pr(>|t|)"]))
}

af   <- fread(AF72)
cov  <- fread(COV)
meta <- fread(META)

# --- normalize AF columns ---
af_pop  <- pick_first(names(af), c("pop","Population","Pop","population"))
af_gene <- pick_first(names(af), c("gene","Gene"))
af_chr  <- pick_first(names(af), c("chr","CHR","chrom","Chrom","chromosome"))
af_pos  <- pick_first(names(af), c("pos","POS","position"))
af_af   <- pick_first(names(af), c("AF","af","nuAF","allele_freq","freq","p"))
if(any(is.na(c(af_pop, af_gene, af_pos, af_af)))) stop("AF72 missing pop/gene/pos/AF")
setnames(af, af_pop,  "pop")
setnames(af, af_gene, "gene")
setnames(af, af_pos,  "pos")
setnames(af, af_af,   "AF")
if(!is.na(af_chr)) setnames(af, af_chr, "chr")
if(!("chr" %in% names(af))) af[, chr := "NA"]

af[, pop := toupper(sub("^[0-9]+_", "", trimws(pop)))]

# --- normalize cov ---
setnames(cov, pick_first(names(cov), c("pop","Population","Pop","population")), "pop")
cov[, pop := toupper(trimws(pop))]
# require mitoPC1/2 & treePC1/2
m1 <- pick_first(names(cov), c("mitoPC1","mtPC1"))
m2 <- pick_first(names(cov), c("mitoPC2","mtPC2"))
if(is.na(m1) || is.na(m2)) stop("COV missing mitoPC1/2")
setnames(cov, m1, "mitoPC1")
setnames(cov, m2, "mitoPC2")
if(!all(c("treePC1","treePC2") %in% names(cov))) stop("COV missing treePC1/2")

# --- normalize meta ---
setnames(meta,
         old=c("Population","Region","Latitude","Longitude","Habitat"),
         new=c("Population","Region","Latitude","Longitude","Habitat"),
         skip_absent=TRUE)
meta[, Population := toupper(trimws(Population))]
meta[, Region := factor(trimws(Region))]
meta[, Habitat := factor(trimws(Habitat))]
meta_env <- unique(meta[, .(Population, Region, Latitude, Longitude, Habitat)])

# --- merge ---
dt <- merge(af, cov[, .(pop, mitoPC1, mitoPC2, treePC1, treePC2)], by="pop", all.x=TRUE)
dt <- merge(dt, meta_env, by.x="pop", by.y="Population", all.x=TRUE)

dt <- dt[!is.na(AF) & !is.na(Habitat) & !is.na(Region) & !is.na(Latitude) & !is.na(Longitude)]
dt <- dt[!is.na(treePC1) & !is.na(treePC2) & !is.na(mitoPC1) & !is.na(mitoPC2)]
dt[, snp_id := paste(chr, pos, gene, sep=":")]

cat(sprintf("[diag] nrow=%d  n_pop=%d  n_snp=%d\n",
            nrow(dt), uniqueN(dt$pop), uniqueN(dt$snp_id)))

# --- models ---
# E1: environment (control tree)
f_E1 <- AF ~ Habitat + Region + Latitude + Longitude + treePC1 + treePC2

# E2: environment (control tree + mito)  -> “纯环境”
f_E2 <- AF ~ Habitat + Region + Latitude + Longitude + mitoPC1 + mitoPC2 + treePC1 + treePC2

# E3: env effect depends on mito background
f_E3 <- AF ~ Habitat*mitoPC1 + Region + Latitude + Longitude + mitoPC2 + treePC1 + treePC2

# --- per-SNP ---
snps <- unique(dt$snp_id)
res_list <- vector("list", length(snps))

for(i in seq_along(snps)){
  sid <- snps[i]
  d <- dt[snp_id == sid]
  if(nrow(d) < 10) next

  fit1 <- try(lm(f_E1, data=d), silent=TRUE)
  fit2 <- try(lm(f_E2, data=d), silent=TRUE)
  fit3 <- try(lm(f_E3, data=d), silent=TRUE)

  # 你最关心：Habitat（整体环境类变量）是否显著
  # 这里用 drop1 做 Habitat 的整体 F-test（比单个水平更合理）
  get_hab_p <- function(fit){
    if(inherits(fit,"try-error")) return(NA_real_)
    a <- try(drop1(fit, test="F"), silent=TRUE)
    if(inherits(a,"try-error")) return(NA_real_)
    if(!("Habitat" %in% rownames(a))) return(NA_real_)
    a["Habitat","Pr(>F)"]
  }

  pHab_E1 <- get_hab_p(fit1)
  pHab_E2 <- get_hab_p(fit2)

  # interaction：Habitat:mitoPC1 的整体检验
  get_int_p <- function(fit){
    if(inherits(fit,"try-error")) return(NA_real_)
    a <- try(drop1(fit, test="F"), silent=TRUE)
    if(inherits(a,"try-error")) return(NA_real_)
    rn <- rownames(a)
    hit <- rn[grepl("^Habitat:mitoPC1$|^mitoPC1:Habitat$", rn)]
    if(length(hit)==0) return(NA_real_)
    a[hit[1],"Pr(>F)"]
  }
  pInt_E3 <- get_int_p(fit3)

  res_list[[i]] <- data.table(
    snp_id=sid, chr=d$chr[1], pos=d$pos[1], gene=d$gene[1], n_pop=nrow(d),
    pHab_E1=pHab_E1,
    pHab_E2=pHab_E2,
    pInt_E3=pInt_E3
  )
}

res <- rbindlist(res_list, fill=TRUE)
res[, BH_Hab_E1 := p.adjust(pHab_E1, method="BH")]
res[, BH_Hab_E2 := p.adjust(pHab_E2, method="BH")]
res[, BH_Int_E3 := p.adjust(pInt_E3, method="BH")]

out_snp <- paste0(OUTPFX, "_perSNP.tsv.gz")
fwrite(res, out_snp, sep="\t")
cat("[ok] wrote:", out_snp, "\n")

# --- gene summary ---
gene_sum <- res[, .(
  n_snps=.N,

  n_hit_Hab_E1 = sum(BH_Hab_E1 < 0.05, na.rm=TRUE),
  prop_hit_Hab_E1 = mean(BH_Hab_E1 < 0.05, na.rm=TRUE),

  n_hit_Hab_E2 = sum(BH_Hab_E2 < 0.05, na.rm=TRUE),
  prop_hit_Hab_E2 = mean(BH_Hab_E2 < 0.05, na.rm=TRUE),

  n_hit_Int_E3 = sum(BH_Int_E3 < 0.05, na.rm=TRUE),
  prop_hit_Int_E3 = mean(BH_Int_E3 < 0.05, na.rm=TRUE)
), by=gene][order(-prop_hit_Hab_E1, -n_hit_Hab_E1, gene)]

out_gene <- paste0(OUTPFX, "_gene_summary.tsv")
fwrite(gene_sum, out_gene, sep="\t")
cat("[ok] wrote:", out_gene, "\n")


#E123 EBLOCK
cat /mnt/spareHD_2/nu_287/_assoc72_subunit/run_hill72_E2_full_terms.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

AF72 <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
COV  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
META <- "/mnt/spareHD_2/nu_287/meta_pop_env.csv"

OUTD   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit"
OUTPFX <- file.path(OUTD, "HILL72_E2_FULLTERMS")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

SEED <- 1
MIN_NPOP <- 10

pick_first <- function(x, candidates){
  hit <- candidates[candidates %in% x]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}

# drop1整体F检验取某一行的p
get_drop1_p <- function(fit, term){
  if(inherits(fit,"try-error") || is.null(fit)) return(NA_real_)
  a <- try(drop1(fit, test="F"), silent=TRUE)
  if(inherits(a,"try-error")) return(NA_real_)
  if(!(term %in% rownames(a))) return(NA_real_)
  a[term, "Pr(>F)"]
}

# 用drop1做一个“块检验”：比较 full vs reduced
get_block_p <- function(d, f_full, f_reduced){
  fit0 <- try(lm(f_reduced, data=d), silent=TRUE)
  fit1 <- try(lm(f_full, data=d), silent=TRUE)
  if(inherits(fit0,"try-error") || inherits(fit1,"try-error")) return(NA_real_)
  a <- try(anova(fit0, fit1), silent=TRUE)
  if(inherits(a,"try-error")) return(NA_real_)
  # 第二行是 added terms 的p
  as.numeric(a$`Pr(>F)`[2])
}

message("[1/5] read inputs")
af   <- fread(AF72)
cov  <- fread(COV)
meta <- fread(META)

# --- normalize AF ---
af_pop  <- pick_first(names(af), c("pop","Population","Pop","population"))
af_gene <- pick_first(names(af), c("gene","Gene"))
af_chr  <- pick_first(names(af), c("chr","CHR","chrom","Chrom","chromosome"))
af_pos  <- pick_first(names(af), c("pos","POS","position"))
af_af   <- pick_first(names(af), c("AF","af","nuAF","allele_freq","freq","p"))
if(any(is.na(c(af_pop, af_gene, af_pos, af_af)))) stop("AF72 missing pop/gene/pos/AF")

setnames(af, af_pop,  "pop")
setnames(af, af_gene, "gene")
setnames(af, af_pos,  "pos")
setnames(af, af_af,   "AF")
if(!is.na(af_chr)) setnames(af, af_chr, "chr")
if(!("chr" %in% names(af))) af[, chr := "NA"]

af[, pop := toupper(sub("^[0-9]+_", "", trimws(pop)))]

# --- normalize COV ---
cov_pop <- pick_first(names(cov), c("pop","Population","Pop","population"))
if(is.na(cov_pop)) stop("COV missing pop")
setnames(cov, cov_pop, "pop")
cov[, pop := toupper(trimws(pop))]

m1 <- pick_first(names(cov), c("mitoPC1","mtPC1"))
m2 <- pick_first(names(cov), c("mitoPC2","mtPC2"))
if(is.na(m1) || is.na(m2)) stop("COV missing mitoPC1/2")
setnames(cov, m1, "mitoPC1")
setnames(cov, m2, "mitoPC2")
if(!all(c("treePC1","treePC2") %in% names(cov))) stop("COV missing treePC1/2")

# --- normalize META ---
setnames(meta,
         old=c("Population","Region","Latitude","Longitude","Habitat"),
         new=c("Population","Region","Latitude","Longitude","Habitat"),
         skip_absent=TRUE)

meta[, Population := toupper(trimws(Population))]
meta[, Region  := factor(trimws(Region))]
meta[, Habitat := factor(trimws(Habitat))]
meta_env <- unique(meta[, .(Population, Region, Latitude, Longitude, Habitat)])

message("[2/5] merge")
dt <- merge(af, cov[, .(pop, mitoPC1, mitoPC2, treePC1, treePC2)], by="pop", all.x=TRUE)
dt <- merge(dt, meta_env, by.x="pop", by.y="Population", all.x=TRUE)

dt <- dt[!is.na(AF) & !is.na(Habitat) & !is.na(Region) & !is.na(Latitude) & !is.na(Longitude)]
dt <- dt[!is.na(treePC1) & !is.na(treePC2) & !is.na(mitoPC1) & !is.na(mitoPC2)]
dt[, snp_id := paste(chr, pos, gene, sep=":")]

message(sprintf("[diag] nrow=%d  n_pop=%d  n_snp=%d",
                nrow(dt), uniqueN(dt$pop), uniqueN(dt$snp_id)))

# ==========================================================
# E2 full model and reduced models for block tests
# ==========================================================
f_E2_full <- AF ~ Habitat + Region + Latitude + Longitude + mitoPC1 + mitoPC2 + treePC1 + treePC2

# 1) Geo block: Latitude+Longitude
f_E2_noGeo <- AF ~ Habitat + Region + mitoPC1 + mitoPC2 + treePC1 + treePC2

# 2) Entire E block: Habitat+Region+Latitude+Longitude
f_E2_noEblock <- AF ~ mitoPC1 + mitoPC2 + treePC1 + treePC2

message("[3/5] per-SNP E2 term tests")
set.seed(SEED)

snps <- unique(dt$snp_id)
res_list <- vector("list", length(snps))

for(i in seq_along(snps)){
  sid <- snps[i]
  d <- dt[snp_id == sid]
  if(nrow(d) < MIN_NPOP) next

  fit <- try(lm(f_E2_full, data=d), silent=TRUE)
  if(inherits(fit,"try-error")) next

  # 单项/块：用 drop1（对factor是整体F检验）
  p_hab <- get_drop1_p(fit, "Habitat")
  p_reg <- get_drop1_p(fit, "Region")
  # 连续变量也可以drop1整体检验（更一致）
  p_lat <- get_drop1_p(fit, "Latitude")
  p_lon <- get_drop1_p(fit, "Longitude")

  # Geo块检验（Lat+Lon一起）
  p_geo <- get_block_p(d, f_E2_full, f_E2_noGeo)

  # 整个E块检验（Hab+Reg+Lat+Lon一起）
  p_Eblk <- get_block_p(d, f_E2_full, f_E2_noEblock)

  res_list[[i]] <- data.table(
    snp_id=sid, chr=d$chr[1], pos=d$pos[1], gene=d$gene[1], n_pop=nrow(d),
    p_Hab_E2=p_hab,
    p_Region_E2=p_reg,
    p_Lat_E2=p_lat,
    p_Lon_E2=p_lon,
    p_GeoBlock_E2=p_geo,
    p_EBlock_E2=p_Eblk
  )

  if(i %% 2000 == 0) message("  ... ", i, "/", length(snps))
}

res <- rbindlist(res_list, fill=TRUE)
if(nrow(res)==0) stop("No results produced. Check merges/labels/filters.")

message("[4/5] BH correction")
res[, BH_Hab_E2      := p.adjust(p_Hab_E2, method="BH")]
res[, BH_Region_E2   := p.adjust(p_Region_E2, method="BH")]
res[, BH_Lat_E2      := p.adjust(p_Lat_E2, method="BH")]
res[, BH_Lon_E2      := p.adjust(p_Lon_E2, method="BH")]
res[, BH_GeoBlock_E2 := p.adjust(p_GeoBlock_E2, method="BH")]
res[, BH_EBlock_E2   := p.adjust(p_EBlock_E2, method="BH")]

out_snp <- paste0(OUTPFX, "_perSNP.tsv.gz")
fwrite(res, out_snp, sep="\t")
message("[ok] wrote: ", out_snp)

# gene summary：每个term各有多少显著SNP
gene_sum <- res[, .(
  n_snps=.N,

  n_hit_Hab_E2      = sum(BH_Hab_E2 < 0.05, na.rm=TRUE),
  prop_hit_Hab_E2   = mean(BH_Hab_E2 < 0.05, na.rm=TRUE),

  n_hit_Region_E2   = sum(BH_Region_E2 < 0.05, na.rm=TRUE),
  prop_hit_Region_E2= mean(BH_Region_E2 < 0.05, na.rm=TRUE),

  n_hit_Lat_E2      = sum(BH_Lat_E2 < 0.05, na.rm=TRUE),
  prop_hit_Lat_E2   = mean(BH_Lat_E2 < 0.05, na.rm=TRUE),

  n_hit_Lon_E2      = sum(BH_Lon_E2 < 0.05, na.rm=TRUE),
  prop_hit_Lon_E2   = mean(BH_Lon_E2 < 0.05, na.rm=TRUE),

  n_hit_GeoBlock_E2 = sum(BH_GeoBlock_E2 < 0.05, na.rm=TRUE),
  prop_hit_GeoBlock_E2 = mean(BH_GeoBlock_E2 < 0.05, na.rm=TRUE),

  n_hit_EBlock_E2   = sum(BH_EBlock_E2 < 0.05, na.rm=TRUE),
  prop_hit_EBlock_E2= mean(BH_EBlock_E2 < 0.05, na.rm=TRUE)
), by=gene][order(-prop_hit_EBlock_E2, -n_hit_EBlock_E2, gene)]

out_gene <- paste0(OUTPFX, "_gene_summary.tsv")
fwrite(gene_sum, out_gene, sep="\t")
message("[ok] wrote: ", out_gene)

message("[5/5] DONE")
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/_assoc72_subunit$ 



#drop 
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ==========================================================
# Main analysis (Hill): mtPC1+mtPC2 + treePC1+treePC2
#   M1: AF ~ mitoPC1 + mitoPC2 + treePC1 + treePC2
#   M2: + Habitat + Region + Latitude + Longitude
#   M3: + mitoPC1 * Habitat interaction (PC2 as covariate)
# Output:
#   per-SNP table (with BH/FDR)
#   per-gene summary table
#
# CHANGE vs original:
#   - Use drop1(F-test) to test the *overall* interaction mitoPC1:Habitat
#   - Still optionally keep one coefficient beta for reference (not used for p/FDR)
# ==========================================================

AF72 <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
COV  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
META <- "/mnt/spareHD_2/nu_287/meta_pop_env.csv"

OUTD <- "/mnt/spareHD_2/nu_287/_assoc72_subunit"
OUTPFX <- file.path(OUTD, "HILL72_MAIN_mt12_tree12_M123_drop1Int")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

SEED <- 1
DROP_RECENT <- FALSE  # TRUE -> only Marine+Freshwater

# -------------------------
# helpers
# -------------------------
pick_first <- function(x, candidates){
  hit <- candidates[candidates %in% x]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}

get_coef_p <- function(fit, term){
  s <- summary(fit)$coefficients
  if(!(term %in% rownames(s))) return(c(NA_real_, NA_real_))
  c(unname(s[term,"Estimate"]), unname(s[term,"Pr(>|t|)"]))
}

# NEW: overall p for interaction term using drop1 (F-test)
get_drop1_p <- function(fit, term_regex){
  if(inherits(fit,"try-error")) return(NA_real_)
  a <- try(drop1(fit, test="F"), silent=TRUE)
  if(inherits(a,"try-error")) return(NA_real_)
  rn <- rownames(a)
  hit <- rn[grepl(term_regex, rn)]
  if(length(hit)==0) return(NA_real_)
  a[hit[1], "Pr(>F)"]
}

# Optional: pick one interaction coefficient just to record a beta (NOT the overall test)
get_any_int_beta_p <- function(fit, term_regex){
  if(inherits(fit,"try-error")) return(c(NA_real_, NA_real_, NA_character_))
  sm <- summary(fit)$coefficients
  rn <- rownames(sm)
  hit <- rn[grepl(term_regex, rn)]
  if(length(hit)==0) return(c(NA_real_, NA_real_, NA_character_))
  trm <- hit[1]
  c(unname(sm[trm,"Estimate"]), unname(sm[trm,"Pr(>|t|)"]), trm)
}

# ==========================================================
# [1/6] read inputs
# ==========================================================
message("[1/6] read AF72: ", AF72)
af <- fread(AF72)

message("[2/6] read COV: ", COV)
cov <- fread(COV)

message("[3/6] read META CSV: ", META)
meta <- fread(META)

# ==========================================================
# normalize AF columns
# ==========================================================
af_pop  <- pick_first(names(af), c("pop","Population","Pop","population"))
af_gene <- pick_first(names(af), c("gene","Gene"))
af_chr  <- pick_first(names(af), c("chr","CHR","chrom","Chrom","chromosome"))
af_pos  <- pick_first(names(af), c("pos","POS","position"))
af_af   <- pick_first(names(af), c("AF","af","nuAF","allele_freq","freq","p"))

if(any(is.na(c(af_pop, af_gene, af_pos, af_af)))){
  stop("AF72 missing required columns: need pop/gene/pos/AF")
}

setnames(af, af_pop,  "pop")
setnames(af, af_gene, "gene")
setnames(af, af_pos,  "pos")
setnames(af, af_af,   "AF")
if(!is.na(af_chr)) setnames(af, af_chr, "chr")
if(!("chr" %in% names(af))) af[, chr := "NA"]

# normalize pop labels: "10_THE" -> "THE"
af[, pop := trimws(as.character(pop))]
af[, pop := sub("^[0-9]+_", "", pop)]

# ==========================================================
# normalize COV columns
# ==========================================================
cov_pop <- pick_first(names(cov), c("pop","Population","Pop","population"))
if(is.na(cov_pop)) stop("COV missing pop column")
setnames(cov, cov_pop, "pop")
cov[, pop := trimws(as.character(pop))]

m1 <- pick_first(names(cov), c("mitoPC1","mtPC1"))
m2 <- pick_first(names(cov), c("mitoPC2","mtPC2"))
if(is.na(m1) || is.na(m2)) stop("COV missing mitoPC1/mitoPC2 columns")
setnames(cov, m1, "mitoPC1")
setnames(cov, m2, "mitoPC2")

need_tree <- c("treePC1","treePC2")
if(!all(need_tree %in% names(cov))){
  stop("COV missing treePC1/treePC2")
}

# ==========================================================
# normalize META columns
# ==========================================================
m_pop <- pick_first(names(meta), c("Population","pop","Pop","population"))
m_reg <- pick_first(names(meta), c("Region","region"))
m_lat <- pick_first(names(meta), c("Latitude","lat","latitude"))
m_lon <- pick_first(names(meta), c("Longitude","lon","longitude"))
m_hab <- pick_first(names(meta), c("Habitat","habitat"))

if(any(is.na(c(m_pop,m_reg,m_lat,m_lon,m_hab)))){
  stop("META missing required cols: Population/Region/Latitude/Longitude/Habitat")
}

setnames(meta, m_pop, "Population")
setnames(meta, m_reg, "Region")
setnames(meta, m_lat, "Latitude")
setnames(meta, m_lon, "Longitude")
setnames(meta, m_hab, "Habitat")

meta[, Population := trimws(as.character(Population))]
meta[, Region     := trimws(as.character(Region))]
meta[, Habitat    := trimws(as.character(Habitat))]

meta_env <- unique(meta[, .(Population, Region, Latitude, Longitude, Habitat)])

# ==========================================================
# [4/6] merge AF + COV + META
# ==========================================================
message("[4/6] merge AF + COV + META")
cov_keep <- cov[, .(pop, mitoPC1, mitoPC2, treePC1, treePC2)]

dt <- merge(af, cov_keep, by="pop", all.x=TRUE)
dt <- merge(dt, meta_env, by.x="pop", by.y="Population", all.x=TRUE)

message(sprintf("[diag] nrow(dt)=%d  n_pop=%d  n_snp=%d",
                nrow(dt), length(unique(dt$pop)),
                length(unique(paste(dt$chr, dt$pos, dt$gene, sep=":")))))
message(sprintf("[diag] NA mitoPC1=%d  NA mitoPC2=%d  NA Habitat=%d  NA Region=%d",
                sum(is.na(dt$mitoPC1)), sum(is.na(dt$mitoPC2)),
                sum(is.na(dt$Habitat)), sum(is.na(dt$Region))))

dt <- dt[!is.na(AF) & !is.na(mitoPC1) & !is.na(mitoPC2)]
dt <- dt[!is.na(treePC1) & !is.na(treePC2)]
dt <- dt[!is.na(Habitat) & !is.na(Region) & !is.na(Latitude) & !is.na(Longitude)]

if(DROP_RECENT){
  dt <- dt[Habitat %in% c("Marine","Freshwater")]
}

dt[, Habitat := factor(Habitat)]
dt[, Region  := factor(Region)]
dt[, snp_id := paste(chr, pos, gene, sep=":")]

# ==========================================================
# define models
# ==========================================================
f_M1 <- AF ~ mitoPC1 + mitoPC2 + treePC1 + treePC2
f_M2 <- AF ~ mitoPC1 + mitoPC2 + Habitat + Region + Latitude + Longitude + treePC1 + treePC2
f_M3 <- AF ~ mitoPC1*Habitat + mitoPC2 + Region + Latitude + Longitude + treePC1 + treePC2

# ==========================================================
# [5/6] per-SNP regression
# ==========================================================
message("[5/6] per-SNP regression (M1/M2/M3; mtPC1+mtPC2; treePC1+treePC2)")
set.seed(SEED)

snps <- unique(dt$snp_id)
res_list <- vector("list", length(snps))

for(i in seq_along(snps)){
  sid <- snps[i]
  d <- dt[snp_id == sid]
  if(nrow(d) < 10) next

  fit1 <- try(lm(f_M1, data=d), silent=TRUE)
  if(inherits(fit1, "try-error")) next
  b1_pc1 <- get_coef_p(fit1, "mitoPC1")
  b1_pc2 <- get_coef_p(fit1, "mitoPC2")

  fit2 <- try(lm(f_M2, data=d), silent=TRUE)
  if(!inherits(fit2, "try-error")){
    b2_pc1 <- get_coef_p(fit2, "mitoPC1")
    b2_pc2 <- get_coef_p(fit2, "mitoPC2")
  } else {
    b2_pc1 <- c(NA,NA); b2_pc2 <- c(NA,NA)
  }

  fit3 <- try(lm(f_M3, data=d), silent=TRUE)
  if(!inherits(fit3, "try-error")){
    b3_pc1 <- get_coef_p(fit3, "mitoPC1")
    b3_pc2 <- get_coef_p(fit3, "mitoPC2")

    # NEW: overall interaction p-value (F-test)
    pInt_M3 <- get_drop1_p(fit3, "^(mitoPC1:Habitat|Habitat:mitoPC1)$")

    # OPTIONAL: record one specific interaction coefficient (beta + p) for debugging/plotting
    any_int <- get_any_int_beta_p(fit3, "(mitoPC1:Habitat|Habitat:mitoPC1)")
    beta_anyInt <- as.numeric(any_int[1])
    p_anyInt    <- as.numeric(any_int[2])
    term_anyInt <- as.character(any_int[3])
  } else {
    b3_pc1 <- c(NA,NA); b3_pc2 <- c(NA,NA)
    pInt_M3 <- NA_real_
    beta_anyInt <- NA_real_; p_anyInt <- NA_real_; term_anyInt <- NA_character_
  }

  res_list[[i]] <- data.table(
    snp_id=sid,
    chr=d$chr[1], pos=d$pos[1], gene=d$gene[1],
    n_pop=nrow(d),

    beta_pc1_M1=b1_pc1[1], p_pc1_M1=b1_pc1[2],
    beta_pc2_M1=b1_pc2[1], p_pc2_M1=b1_pc2[2],

    beta_pc1_M2=b2_pc1[1], p_pc1_M2=b2_pc1[2],
    beta_pc2_M2=b2_pc2[1], p_pc2_M2=b2_pc2[2],

    beta_pc1_M3=b3_pc1[1], p_pc1_M3=b3_pc1[2],
    beta_pc2_M3=b3_pc2[1], p_pc2_M3=b3_pc2[2],

    # NEW outputs:
    pInt_M3 = pInt_M3,              # overall interaction test (recommended)
    any_int_term = term_anyInt,     # one coefficient name (optional)
    beta_anyInt  = beta_anyInt,     # one coefficient beta (optional)
    p_anyInt     = p_anyInt         # one coefficient p (optional)
  )

  if(i %% 2000 == 0) message("  ... ", i, "/", length(snps))
}

res <- rbindlist(res_list, fill=TRUE)
if(nrow(res)==0) stop("No results produced. Check merges/pop labels/filters.")

# BH/FDR
res[, BH_pc1_M1 := p.adjust(p_pc1_M1, method="BH")]
res[, BH_pc2_M1 := p.adjust(p_pc2_M1, method="BH")]
res[, BH_pc1_M2 := p.adjust(p_pc1_M2, method="BH")]
res[, BH_pc2_M2 := p.adjust(p_pc2_M2, method="BH")]
res[, BH_pc1_M3 := p.adjust(p_pc1_M3, method="BH")]
res[, BH_pc2_M3 := p.adjust(p_pc2_M3, method="BH")]

# NEW: BH for overall interaction
res[, BH_Int_M3 := p.adjust(pInt_M3, method="BH")]

out_snp <- paste0(OUTPFX, "_perSNP.tsv.gz")
fwrite(res, out_snp, sep="\t")
message("[ok] wrote: ", out_snp)

# ==========================================================
# [6/6] gene-level summary
# ==========================================================
gene_sum <- res[, .(
  n_snps=.N,

  n_hit_pc1_M1=sum(BH_pc1_M1 < 0.05, na.rm=TRUE),
  prop_hit_pc1_M1=mean(BH_pc1_M1 < 0.05, na.rm=TRUE),
  n_hit_pc2_M1=sum(BH_pc2_M1 < 0.05, na.rm=TRUE),
  prop_hit_pc2_M1=mean(BH_pc2_M1 < 0.05, na.rm=TRUE),

  n_hit_pc1_M2=sum(BH_pc1_M2 < 0.05, na.rm=TRUE),
  prop_hit_pc1_M2=mean(BH_pc1_M2 < 0.05, na.rm=TRUE),
  n_hit_pc2_M2=sum(BH_pc2_M2 < 0.05, na.rm=TRUE),
  prop_hit_pc2_M2=mean(BH_pc2_M2 < 0.05, na.rm=TRUE),

  # NEW: overall interaction hits
  n_hit_Int_M3=sum(BH_Int_M3 < 0.05, na.rm=TRUE),
  prop_hit_Int_M3=mean(BH_Int_M3 < 0.05, na.rm=TRUE)
), by=gene][order(-prop_hit_pc1_M1, -n_hit_pc1_M1, gene)]

out_gene <- paste0(OUTPFX, "_gene_summary.tsv")
fwrite(gene_sum, out_gene, sep="\t")
message("[ok] wrote: ", out_gene)

message("[DONE] Main Hill analysis finished.")




#e12 table
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

E1 <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/HILL72_E1_FULLTERMS_perSNP.tsv.gz"
E2 <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/HILL72_E2_FULLTERMS_perSNP.tsv.gz"

count_one <- function(f, label){
  x <- fread(f)

  # 必备列
  need <- c("snp_id","BH_Hab_E1","BH_GeoBlock_E1","BH_EBlock_E1")
  if(grepl("_E2_", f) || label=="E2"){
    need <- c("snp_id","BH_Hab_E2","BH_GeoBlock_E2","BH_EBlock_E2")
  }
  miss <- setdiff(need, names(x))
  if(length(miss)>0) stop(label, " missing: ", paste(miss, collapse=", "))

  if(label=="E1"){
    n_total <- nrow(x)
    n_hab   <- sum(x$BH_Hab_E1      < 0.05, na.rm=TRUE)
    n_geo   <- sum(x$BH_GeoBlock_E1 < 0.05, na.rm=TRUE)
    n_eblk  <- sum(x$BH_EBlock_E1   < 0.05, na.rm=TRUE)
  } else {
    n_total <- nrow(x)
    n_hab   <- sum(x$BH_Hab_E2      < 0.05, na.rm=TRUE)
    n_geo   <- sum(x$BH_GeoBlock_E2 < 0.05, na.rm=TRUE)
    n_eblk  <- sum(x$BH_EBlock_E2   < 0.05, na.rm=TRUE)
  }

  data.table(
    model = label,
    n_total = n_total,
    n_sig_Hab = n_hab,   prop_sig_Hab = n_hab/n_total,
    n_sig_GeoBlock = n_geo, prop_sig_GeoBlock = n_geo/n_total,
    n_sig_EBlock = n_eblk,  prop_sig_EBlock = n_eblk/n_total
  )
}

res <- rbind(
  count_one(E1, "E1"),
  count_one(E2, "E2")
)

print(res)

cat("\nKey readout (EBlock):\n")
cat(sprintf("E1: %d / %d = %.4f\n", res[model=="E1"]$n_sig_EBlock, res[model=="E1"]$n_total, res[model=="E1"]$prop_sig_EBlock))
cat(sprintf("E2: %d / %d = %.4f\n", res[model=="E2"]$n_sig_EBlock, res[model=="E2"]$n_total, res[model=="E2"]$prop_sig_EBlock))



#统计
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ========== inputs ==========
MT <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/HILL72_MAIN_mt12_tree12_M123_drop1Int_perSNP.tsv.gz"
E  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/HILL72_E2_FULLTERMS_perSNP.tsv.gz"

OUTD <- "/mnt/spareHD_2/nu_287/_assoc72_subunit"
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

OUT_TSV <- file.path(OUTD, "HILL72_SNPLEVEL_VENN_mtPC1_vs_EBlock_counts.tsv")
OUT_PNG <- file.path(OUTD, "HILL72_SNPLEVEL_VENN_mtPC1_vs_EBlock.png")
OUT_PDF <- file.path(OUTD, "HILL72_SNPLEVEL_VENN_mtPC1_vs_EBlock.pdf")

# ========== read ==========
m <- fread(MT)
e <- fread(E)

# sanity check
stopifnot(all(c("snp_id","BH_pc1_M2") %in% names(m)))
stopifnot(all(c("snp_id","BH_EBlock_E2") %in% names(e)))

# merge: only SNPs present in both result tables
dt <- merge(
  m[, .(snp_id, BH_pc1_M2)],
  e[, .(snp_id, BH_EBlock_E2)],
  by="snp_id",
  all=FALSE
)

# flags
dt[, mt_sig  := BH_pc1_M2 < 0.05]
dt[, env_sig := BH_EBlock_E2 < 0.05]

# counts
mt_only <- dt[ mt_sig & !env_sig, .N]
e_only  <- dt[!mt_sig &  env_sig, .N]
overlap <- dt[ mt_sig &  env_sig, .N]
none    <- dt[!mt_sig & !env_sig, .N]
total   <- nrow(dt)

cat(sprintf("[diag] merged SNPs = %d\n", total))
cat(sprintf("mt_only=%d  E_only=%d  overlap=%d  none=%d\n", mt_only, e_only, overlap, none))

# write counts table
res <- data.table(
  category=c("mt_only","E_only","overlap(mt&E)","none","TOTAL"),
  n=c(mt_only, e_only, overlap, none, total)
)
fwrite(res, OUT_TSV, sep="\t")
cat("[ok] wrote:", OUT_TSV, "\n")

# ========== plot Venn (proportional) ==========
# Prefer eulerr (proportional circles)
if(!requireNamespace("eulerr", quietly=TRUE)) {
  message("Installing eulerr ...")
  install.packages("eulerr", repos="https://cloud.r-project.org")
}
library(eulerr)

fit <- eulerr::euler(c(
  "mtPC1"        = mt_only + overlap,
  "Eblock"       = e_only + overlap,
  "mtPC1&Eblock" = overlap
))

# PNG
png(OUT_PNG, width=1400, height=1100, res=200)
plot(
  fit,
  quantities = list(type = "counts"),
  labels = list(font = 2),
  main = "SNP-level overlap: mtPC1 (BH<0.05) vs E-block (BH<0.05)\n(E-block = Habitat + Region + Lat + Lon; model E2 controls mitoPC1/2 + treePC1/2)"
)
dev.off()
cat("[ok] wrote:", OUT_PNG, "\n")

# PDF
pdf(OUT_PDF, width=7, height=6)
plot(
  fit,
  quantities = list(type = "counts"),
  labels = list(font = 2),
  main = "SNP-level overlap: mtPC1 vs E-block"
)
dev.off()
cat("[ok] wrote:", OUT_PDF, "\n")









#table mixed dominant
library(data.table)

# -------- input --------
M_file    <- "HILL72_MAIN_mt12_tree12_M123_gene_summary.tsv"
E_file    <- "HILL72_ENVONLY_E123_gene_summary.tsv"
Eblk_file <- "HILL72_E2_FULLTERMS_gene_summary.tsv"  # 你有的话就保留

M    <- fread(M_file)
E    <- fread(E_file)
Eblk <- fread(Eblk_file)

# -------- sanity: ensure key column names --------
stopifnot("gene" %in% names(M))
stopifnot("gene" %in% names(E))
stopifnot("gene" %in% names(Eblk))

# -------- pick columns (match YOUR M table) --------
M_keep <- M[, .(
  gene,
  n_snps,
  prop_hit_pc1_M1,
  prop_hit_pc1_M2,
  prop_hit_pc1xHab   # <- 你的 M3 interaction 命中比例
)]

# -------- E table: try common column names; print names if missing --------
need_E <- c("prop_hit_Hab_E1","prop_hit_Hab_E2","prop_hit_Int_E3")
miss_E <- setdiff(need_E, names(E))
if(length(miss_E) > 0){
  cat("[E] missing:", paste(miss_E, collapse=", "), "\n")
  cat("[E] available columns:\n")
  print(names(E))
  stop("E gene summary column names do not match. See printed names(E).")
}
E_keep <- E[, .(
  gene,
  prop_hit_Hab_E1,
  prop_hit_Hab_E2,
  prop_hit_Int_E3
)]

# -------- E2 full terms (E-block) --------
need_B <- c("prop_hit_EBlock_E2")
miss_B <- setdiff(need_B, names(Eblk))
if(length(miss_B) > 0){
  cat("[Eblk] missing:", paste(miss_B, collapse=", "), "\n")
  cat("[Eblk] available columns:\n")
  print(names(Eblk))
  stop("E2_FULLTERMS gene summary column names do not match. See printed names(Eblk).")
}
Eblk_keep <- Eblk[, .(
  gene,
  prop_hit_EBlock_E2
)]

# -------- merge --------
res <- Reduce(function(x,y) merge(x,y, by="gene", all=TRUE),
              list(M_keep, E_keep, Eblk_keep))

# NA -> 0 for prop columns
prop_cols <- setdiff(names(res), c("gene","n_snps"))
for (j in prop_cols) set(res, which(is.na(res[[j]])), j, 0)

# classify (optional)
res[, signal_class := fifelse(
  prop_hit_pc1_M1 > 0 & prop_hit_Hab_E2 == 0, "Mitonuclear-dominant",
  fifelse(prop_hit_pc1_M1 > 0 & prop_hit_Hab_E2 > 0, "Mixed (mt + env)",
          fifelse(prop_hit_pc1_M1 == 0 & prop_hit_Hab_E2 > 0, "Environment-only", "No signal")))
]

# sort for Results
setorder(res, -prop_hit_pc1_M1, -prop_hit_pc1xHab, -prop_hit_Hab_E2, gene)

fwrite(res, "HILL72_gene_summary_M_E_merged.tsv", sep="\t")
cat("[ok] wrote HILL72_gene_summary_M_E_merged.tsv\n")



library(data.table)

# 读入三张 perSNP 表
M <- fread("HILL72_MAIN_mt12_tree12_M123_perSNP.tsv.gz")
E <- fread("HILL72_ENVONLY_E123_perSNP.tsv.gz")
B <- fread("HILL72_E2_FULLTERMS_perSNP.tsv.gz")

# 统一 key
setkey(M, snp_id)
setkey(E, snp_id)
setkey(B, snp_id)

# merge
X <- merge(M, E, by="snp_id", all=TRUE)
X <- merge(X, B, by="snp_id", all=TRUE)

# 给每个 SNP 打标签（你可以按需要调整优先级）
X[, label := fifelse(!is.na(BH_pc1_M1) & BH_pc1_M1 < 0.05, "mitoPC1",
              fifelse(!is.na(BH_pc2_M1) & BH_pc2_M1 < 0.05, "mitoPC2",
              fifelse(!is.na(BH_Hab_E2) & BH_Hab_E2 < 0.05, "Habitat_independent_E2",
              fifelse(!is.na(BH_EBlock_E2) & BH_EBlock_E2 < 0.05, "EnvBlock_independent_E2",
              fifelse(!is.na(BH_Int_E3) & BH_Int_E3 < 0.05, "Habitat×mitoPC1_E3",
              "none")))))]

# 输出一个可追溯表：每个 SNP 到底显著在哪
fwrite(X, "HILL72_perSNP_traceable_labels.tsv.gz", sep="\t")


library(data.table)

# ---- read ----
M <- fread("HILL72_MAIN_mt12_tree12_M123_perSNP.tsv.gz")
E <- fread("HILL72_ENVONLY_E123_perSNP.tsv.gz")
B <- fread("HILL72_E2_FULLTERMS_perSNP.tsv.gz")

# ---- helper: pick first existing column ----
pick <- function(dt, cand){
  hit <- cand[cand %in% names(dt)]
  if(length(hit)==0) NA_character_ else hit[1]
}

# ---- choose BH columns robustly ----
# mt hits: prefer M1 (main mitonuclear), fall back to any BH pc1/pc2
bh_pc1 <- pick(M, c("BH_pc1_M1","BH_mitoPC1_M1","BH_PC1_M1","BH_pc1"))
bh_pc2 <- pick(M, c("BH_pc2_M1","BH_mitoPC2_M1","BH_PC2_M1","BH_pc2"))

if(is.na(bh_pc1) || is.na(bh_pc2)){
  stop("Cannot find BH columns for mitoPC1/2 in M. Check names(M).")
}

# env independent hits: prefer E2 fullterms block (strongest definition)
# B has BH_Hab_E2 / BH_EBlock_E2 typically
bh_hab_E2   <- pick(B, c("BH_Hab_E2","BH_Habitat_E2"))
bh_Eblk_E2  <- pick(B, c("BH_EBlock_E2","BH_Eblk_E2","BH_EBlock"))

if(is.na(bh_hab_E2) && is.na(bh_Eblk_E2)){
  # fallback: use E table (envonly) if fullterms not available
  bh_hab_E2  <- pick(E, c("BH_Hab_E2","BH_Habitat_E2"))
  if(is.na(bh_hab_E2)){
    stop("Cannot find BH env columns (Hab/EBlock) in B or E. Check names(B), names(E).")
  }
}

# ---- merge (keep gene/chr/pos if present) ----
setkey(M, snp_id)
setkey(E, snp_id)
setkey(B, snp_id)

X <- merge(M, E, by="snp_id", all=TRUE, suffixes=c(".M",".E"))
X <- merge(X, B, by="snp_id", all=TRUE, suffixes=c("",".B"))

# ---- ensure gene column exists ----
gene_col <- pick(X, c("gene","gene.M","gene.E","gene.B"))
if(is.na(gene_col)) stop("No gene column found after merge.")
setnames(X, gene_col, "gene")

# ---- define hits ----
X[, mt_hit  := (get(bh_pc1) < 0.05) | (get(bh_pc2) < 0.05)]

if(!is.na(bh_Eblk_E2)){
  X[, env_hit := (get(bh_hab_E2) < 0.05) | (get(bh_Eblk_E2) < 0.05)]
} else {
  # only have habitat BH
  X[, env_hit := (get(bh_hab_E2) < 0.05)]
}

X[, both_hit := mt_hit & env_hit]

# optional: label
X[, label := fifelse(both_hit, "both",
              fifelse(mt_hit, "mt_only",
              fifelse(env_hit, "env_only", "none")))]

# ---- global summary ----
global <- X[, .(
  n_snps = .N,
  n_mt   = sum(mt_hit, na.rm=TRUE),
  n_env  = sum(env_hit, na.rm=TRUE),
  n_both = sum(both_hit, na.rm=TRUE),
  n_none = sum(label=="none", na.rm=TRUE),
  prop_mt   = mean(mt_hit, na.rm=TRUE),
  prop_env  = mean(env_hit, na.rm=TRUE),
  prop_both = mean(both_hit, na.rm=TRUE)
)]

print(global)

# ---- per-gene summary ----
by_gene <- X[, .(
  n_snps = .N,
  n_mt   = sum(mt_hit, na.rm=TRUE),
  n_env  = sum(env_hit, na.rm=TRUE),
  n_both = sum(both_hit, na.rm=TRUE),
  n_none = sum(label=="none", na.rm=TRUE),
  prop_mt   = mean(mt_hit, na.rm=TRUE),
  prop_env  = mean(env_hit, na.rm=TRUE),
  prop_both = mean(both_hit, na.rm=TRUE)
), by=gene][order(-prop_mt, -n_mt, gene)]

fwrite(by_gene, "HILL72_perGene_mt_env_overlap_summary.tsv", sep="\t")
fwrite(X[, .(snp_id, gene, label, mt_hit, env_hit, both_hit)], 
       "HILL72_perSNP_mt_env_labels.tsv.gz", sep="\t")

cat("\n[ok] wrote:\n  HILL72_perGene_mt_env_overlap_summary.tsv\n  HILL72_perSNP_mt_env_labels.tsv.gz\n")
cat("\n[diag] Using columns:\n")
cat("  mt BH:", bh_pc1, ",", bh_pc2, "\n")
cat("  env BH:", bh_hab_E2, if(!is.na(bh_Eblk_E2)) paste0(", ", bh_Eblk_E2) else "", "\n")






