
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