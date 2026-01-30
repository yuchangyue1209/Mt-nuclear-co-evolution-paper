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




















#improved  mt~E
#!/usr/bin/env Rscript
# ==========================================================
# mt ~ Environment (meta level), joint mtPC1+mtPC2
# - Fits:
#   (E0) mtPC1/2 ~ Habitat + Region + Latitude + Longitude
#   (E1) mtPC1/2 ~ Habitat + Region + Latitude + Longitude + treePC1 + treePC2
# - Outputs:
#   1) OUT_LM:   per-PC linear models (PC1 and PC2 separately)
#   2) OUT_MAN:  joint MANOVA tests for (mtPC1, mtPC2) together
#      (recommend reading Pillai first)
# ==========================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- inputs ----
COV  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
META <- "/mnt/spareHD_2/nu_287/meta_pop_env.csv"

# ---- outputs ----
OUT_LM  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mt_env_lm_PC12.tsv"
OUT_MAN <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mt_env_manova_PC12.tsv"

# ==========================================================
# helpers
# ==========================================================
pick_first <- function(x, candidates){
  hit <- candidates[candidates %in% x]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}

norm_pop <- function(x){
  x <- toupper(trimws(as.character(x)))
  # tolerate "10_THE" / "10_THE_S10" style -> THE
  x <- sub("^[0-9]+_", "", x)
  x <- sub("_S[0-9]+$", "", x)
  x
}

grab_lm <- function(fit, tag){
  s <- summary(fit)
  co <- as.data.table(s$coefficients, keep.rownames="term")
  setnames(co,
           old=c("Estimate","Std. Error","t value","Pr(>|t|)"),
           new=c("beta","se","t","p"))
  co[, model := tag]
  co[, adjR2 := s$adj.r.squared]
  co[, R2    := s$r.squared]
  co[, n     := nobs(fit)]
  co[]
}

# Robust MANOVA extraction: keep only common columns across tests
grab_manova <- function(man, tag){
  tests <- c("Pillai","Wilks","Hotelling-Lawley","Roy")
  out <- rbindlist(lapply(tests, function(tt){
    a <- summary(man, test=tt)$stats
    # keep shared columns only (avoid data.table rbind column mismatch warnings)
    a <- a[, c("approx F","num Df","den Df","Pr(>F)"), drop=FALSE]
    x <- as.data.table(a, keep.rownames="term")
    setnames(x,
             old=c("approx F","num Df","den Df","Pr(>F)"),
             new=c("F","numDf","denDf","p"))
    x[, test  := tt]
    x[, model := tag]
    x[, n     := nrow(model.frame(man))]
    x[]
  }), fill=TRUE, use.names=TRUE)
  out
}

# ==========================================================
# read + normalize
# ==========================================================
cov  <- fread(COV)
meta <- fread(META)

# cov pop
cov_pop <- pick_first(names(cov), c("pop","Population","Pop","population"))
if(is.na(cov_pop)) stop("COV missing pop column")
setnames(cov, cov_pop, "pop")
cov[, pop := norm_pop(pop)]

# required PCs
m1 <- pick_first(names(cov), c("mitoPC1","mtPC1"))
m2 <- pick_first(names(cov), c("mitoPC2","mtPC2"))
if(is.na(m1) || is.na(m2)) stop("COV missing mitoPC1/mitoPC2")
setnames(cov, m1, "mitoPC1")
setnames(cov, m2, "mitoPC2")

need_tree <- c("treePC1","treePC2")
if(!all(need_tree %in% names(cov))) stop("COV missing treePC1/treePC2")

# meta columns
m_pop <- pick_first(names(meta), c("Population","pop","Pop","population"))
m_reg <- pick_first(names(meta), c("Region","region"))
m_lat <- pick_first(names(meta), c("Latitude","lat","latitude"))
m_lon <- pick_first(names(meta), c("Longitude","lon","longitude"))
m_hab <- pick_first(names(meta), c("Habitat","habitat"))

if(any(is.na(c(m_pop,m_reg,m_lat,m_lon,m_hab)))){
  stop("META missing required cols: Population/Region/Latitude/Longitude/Habitat")
}

setnames(meta, m_pop, "pop")
setnames(meta, m_reg, "Region")
setnames(meta, m_lat, "Latitude")
setnames(meta, m_lon, "Longitude")
setnames(meta, m_hab, "Habitat")

meta[, pop := norm_pop(pop)]
meta[, Region  := factor(trimws(as.character(Region)))]
meta[, Habitat := factor(trimws(as.character(Habitat)))]
meta[, Latitude  := as.numeric(Latitude)]
meta[, Longitude := as.numeric(Longitude)]

meta_env <- unique(meta[, .(pop, Region, Latitude, Longitude, Habitat)])

# ==========================================================
# merge
# ==========================================================
dt <- merge(
  cov[, .(pop, mitoPC1, mitoPC2, treePC1, treePC2)],
  meta_env,
  by="pop", all=FALSE
)

# hard filters
dt <- dt[is.finite(mitoPC1) & is.finite(mitoPC2)]
dt <- dt[is.finite(treePC1) & is.finite(treePC2)]
dt <- dt[is.finite(Latitude) & is.finite(Longitude)]
dt <- dt[!is.na(Habitat) & !is.na(Region)]

if(nrow(dt) < 10) stop("Too few rows after merge/filtering.")

cat(sprintf("[diag] n=%d pops=%d\n", nrow(dt), uniqueN(dt$pop)))
cat("[diag] Habitat levels:", paste(levels(dt$Habitat), collapse=", "), "\n")
cat("[diag] Region  levels:", paste(levels(dt$Region), collapse=", "), "\n")

# ==========================================================
# models
# ==========================================================
# E0: env only
f_pc1_E0 <- mitoPC1 ~ Habitat + Region + Latitude + Longitude
f_pc2_E0 <- mitoPC2 ~ Habitat + Region + Latitude + Longitude

# E1: env + tree background
f_pc1_E1 <- mitoPC1 ~ Habitat + Region + Latitude + Longitude + treePC1 + treePC2
f_pc2_E1 <- mitoPC2 ~ Habitat + Region + Latitude + Longitude + treePC1 + treePC2

# joint MANOVA formulas
f_joint_E0 <- cbind(mitoPC1, mitoPC2) ~ Habitat + Region + Latitude + Longitude
f_joint_E1 <- cbind(mitoPC1, mitoPC2) ~ Habitat + Region + Latitude + Longitude + treePC1 + treePC2

# ==========================================================
# fit (LM)
# ==========================================================
m1_pc1 <- lm(f_pc1_E0, data=dt)
m1_pc2 <- lm(f_pc2_E0, data=dt)
m2_pc1 <- lm(f_pc1_E1, data=dt)
m2_pc2 <- lm(f_pc2_E1, data=dt)

lm_res <- rbindlist(list(
  grab_lm(m1_pc1, "PC1_E0_env"),
  grab_lm(m2_pc1, "PC1_E1_env+tree"),
  grab_lm(m1_pc2, "PC2_E0_env"),
  grab_lm(m2_pc2, "PC2_E1_env+tree")
), fill=TRUE, use.names=TRUE)

fwrite(lm_res, OUT_LM, sep="\t")
cat("[ok] wrote lm:", OUT_LM, "\n")

# ==========================================================
# fit (MANOVA)
# ==========================================================
man0 <- manova(f_joint_E0, data=dt)
man1 <- manova(f_joint_E1, data=dt)

man_res <- rbindlist(list(
  grab_manova(man0, "PC12_E0_env"),
  grab_manova(man1, "PC12_E1_env+tree")
), fill=TRUE, use.names=TRUE)

fwrite(man_res, OUT_MAN, sep="\t")
cat("[ok] wrote manova:", OUT_MAN, "\n\n")

cat("How to read:\n")
cat("1) 看 OUT_MAN 里 Pillai test：Habitat/Region/Lat/Lon 在 PC12 的 joint p 是否显著。\n")
cat("2) 如果 joint 显著，再看 OUT_LM：是 PC1 还是 PC2 在驱动（以及方向/系数）。\n")
cat("3) 对比 E0 vs E1：如果 E1 显著性大幅消失，说明 mt-env 主要是树/空间结构造成。\n")





#mt～ e   distance to mean of 2 marine

cd /work/cyu/poolseq/PPalign_output/ann

cat > mt_sync_gene_dist.py <<'PY'
#!/usr/bin/env python3
import sys, glob, math
from pathlib import Path

def parse_counts(s):
    # PoPoolation sync field: "A:T:C:G:N:del"
    a,t,c,g,n,d = s.split(":")
    return list(map(int, (a,t,c,g,n,d)))

def freqs_4(ct):
    a,t,c,g,n,d = ct
    cov = a+t+c+g
    if cov == 0:
        return None
    return [a/cov, c/cov, g/cov, t/cov]  # A,C,G,T

def euclid(x, y):
    return math.sqrt(sum((xi-yi)**2 for xi, yi in zip(x,y)))

def mean_vec(v1, v2):
    return [(a+b)/2.0 for a,b in zip(v1,v2)]

if len(sys.argv) < 5:
    print("Usage: mt_sync_gene_dist.py '<glob_sync>' mt_pops.txt RS SAY > mt_gene_dist.tsv", file=sys.stderr)
    sys.exit(1)

glob_pat = sys.argv[1]   # "*_fixed.sync"
pop_file = sys.argv[2]   # mt_pops.txt
rs_name  = sys.argv[3]   # RS
say_name = sys.argv[4]   # SAY

pops = [l.strip() for l in open(pop_file) if l.strip()]
pop_idx = {p:i for i,p in enumerate(pops)}
if rs_name not in pop_idx or say_name not in pop_idx:
    raise SystemExit(f"[ERROR] RS/SAY not found in mt_pops.txt. Example pops: {pops[:10]}")

sync_files = sorted(glob.glob(glob_pat))
if not sync_files:
    raise SystemExit(f"[ERROR] No files matched: {glob_pat}")

print("gene\tpop\tdist_to_marineMean\tn_sites_used")

for fp in sync_files:
    gene = Path(fp).name
    gene = gene.replace("_fixed.sync","").replace(".sync","")

    vecs = {p: [] for p in pops}
    used = 0

    with open(fp) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")

            # Most common: chr pos ref pop1 pop2 ...
            counts = parts[3:3+len(pops)]
            # Fallback: chr pos ref extra pop1 pop2 ...
            if len(counts) < len(pops):
                counts = parts[4:4+len(pops)]
            if len(counts) < len(pops):
                continue

            site_freqs = []
            ok = True
            for cf in counts:
                fr = freqs_4(parse_counts(cf))
                if fr is None:
                    ok = False
                    break
                site_freqs.append(fr)
            if not ok:
                continue

            used += 1
            for p, fr in zip(pops, site_freqs):
                vecs[p].extend(fr)

    if used == 0:
        for p in pops:
            print(f"{gene}\t{p}\tNA\t0")
        continue

    marine = mean_vec(vecs[rs_name], vecs[say_name])

    for p in pops:
        d = euclid(vecs[p], marine) if len(vecs[p])==len(marine) else float("nan")
        print(f"{gene}\t{p}\t{d}\t{used}")
PY

chmod +x mt_sync_gene_dist.py


./mt_sync_gene_dist.py "*_fixed.sync" mt_pops.txt RS SAY > mt_gene_dist.tsv
head mt_gene_dist.tsv


awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{print $0,"dist_norm"; next}
{
  n=$4;
  if($3=="NA" || n==0){print $0,"NA"; next}
  norm = $3 / sqrt(4*n);
  print $0, norm
}' mt_gene_dist.tsv > mt_gene_dist.norm.tsv


library(data.table)

d <- fread("mt_gene_dist.norm.tsv")
meta <- fread("/mnt/spareHD_2/nu_287/meta_pop_env.csv")  # 你已有

x <- merge(d, meta, by="pop")
x <- x[is.finite(dist_norm)]

fit_one <- function(df){
  m0 <- lm(dist_norm ~ Habitat + Region + Latitude + Longitude, data=df)
  m1 <- lm(dist_norm ~ Habitat + Region + Latitude + Longitude + treePC1 + treePC2, data=df)

  p0 <- drop1(m0, test="F")
  p1 <- drop1(m1, test="F")

  data.table(
    gene=df$gene[1],
    n=nrow(df),
    adjR2_0=summary(m0)$adj.r.squared,
    adjR2_1=summary(m1)$adj.r.squared,
    pHabitat_0=p0["Habitat","Pr(>F)"],
    pRegion_0 =p0["Region","Pr(>F)"],
    pHabitat_1=p1["Habitat","Pr(>F)"],
    pRegion_1 =p1["Region","Pr(>F)"]
  )
}

res <- x[, fit_one(.SD), by=gene]
res[, qHabitat_0 := p.adjust(pHabitat_0, "BH")]
res[, qHabitat_1 := p.adjust(pHabitat_1, "BH")]
fwrite(res, "mt_gene_env_assoc.tsv", sep="\t")







#mt～e    pergene pi



#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

IN  <- "/mnt/spareHD_2/nu_287/Poolinfo_with_13gene_details.csv"
OUT <- "/mnt/spareHD_2/nu_287/mt_pi13_env_assoc.tsv"

dt <- fread(IN)

# ====== 基本清理 ======
# 确保 Habitat/Region 是因子（你的数据里 Region = Alaska/BC）
dt[, Habitat := factor(Habitat)]
dt[, Region  := factor(Region)]

# 可选：去掉 Marine 只看 FW/Recent（你想怎么做都行）
# dt <- dt[Habitat != "Marine"]

# ====== 把 13 个基因的 pi 列拉成长表 ======
pi_cols <- grep("^pi_(ND1|ND2|ND3|ND4L|ND4|ND5|ND6|COX1|COX2|COX3|CYTB|ATP6|ATP8)$",
                names(dt), value = TRUE)

long <- melt(
  dt,
  id.vars = c("Population","Region","Latitude","Longitude","Habitat"),
  measure.vars = pi_cols,
  variable.name = "gene",
  value.name = "pi"
)

# gene 名从 "pi_ND5" -> "ND5"
long[, gene := sub("^pi_", "", gene)]

# 去掉缺失
long <- long[is.finite(pi)]

# ====== 每个 gene 跑一个 lm：pi ~ Habitat + Region + Lat + Lon ======
fit_one <- function(x) {
  m <- lm(pi ~ Habitat + Region + Latitude + Longitude, data = x)
  a <- anova(m)

  # 给每个 term 一个 p（和你之前 dist~E 的风格一致）
  data.table(
    gene = unique(x$gene),
    n = nrow(x),
    adjR2 = summary(m)$adj.r.squared,
    pHabitat = a["Habitat", "Pr(>F)"],
    pRegion  = a["Region",  "Pr(>F)"],
    pLat     = a["Latitude","Pr(>F)"],
    pLon     = a["Longitude","Pr(>F)"]
  )
}

res <- long[, fit_one(.SD), by = gene]

# ====== 多重校正（在 13 个基因内）=====
res[, qHabitat := p.adjust(pHabitat, "BH")]
res[, qRegion  := p.adjust(pRegion,  "BH")]
res[, qLat     := p.adjust(pLat,     "BH")]
res[, qLon     := p.adjust(pLon,     "BH")]

setorder(res, qHabitat, qRegion)

fwrite(res, OUT, sep = "\t")
cat("Wrote:", OUT, "\n")
print(res)
