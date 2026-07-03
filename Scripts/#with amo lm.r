#with amo lm


#/mnt/spareHD_2/nu_287/
#MT_TREE=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile \
##OUT_PC=/mnt/spareHD_2/nu_287/meta/mtPC_from_mtTree_noDloop_withAMO.tsv.gz \
# 00_make_mtPC_from_tree.R
Rscript - <<'RS'
library(ape)

tr <- read.tree("/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile")

D <- cophenetic(tr)
pco <- cmdscale(D, k=5, eig=TRUE)

expl <- 100 * pco$eig / sum(pco$eig)
cat(sprintf("Explained variance (with AMO): PC1=%.2f%% PC2=%.2f%% PC3=%.2f%% PC4=%.2f%% PC5=%.2f%%\n",
            expl[1], expl[2], expl[3], expl[4], expl[5]))
RS
Warning message:
package ‘ape’ was built under R version 4.4.2 
Explained variance (with AMO): PC1=52.52% PC2=26.84% PC3=20.28% PC4=5.21% PC5=2.41%

#mtpc plot
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

mt_tree <- "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile"
out_dir <- "/mnt/spareHD_2/nu_287/_assoc72_subunit"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

tr <- read.tree(mt_tree)

D <- cophenetic(tr)
pco <- cmdscale(D, k=5, eig=TRUE)

expl <- 100 * pco$eig / sum(pco$eig)
cat(sprintf("[mitoPC explained withAMO] PC1=%.2f%% PC2=%.2f%% PC3=%.2f%% PC4=%.2f%% PC5=%.2f%%\n",
            expl[1], expl[2], expl[3], expl[4], expl[5]))

PC <- data.frame(
  pop    = rownames(pco$points),
  mitoPC1 = pco$points[,1],
  mitoPC2 = pco$points[,2],
  mitoPC3 = pco$points[,3],
  mitoPC4 = pco$points[,4],
  mitoPC5 = pco$points[,5]
)

PC <- PC %>%
  mutate(region = case_when(
    pop %in% c("FG","LG","SR","SL","TL","WB","WT","WK","LB") ~ "Alaska",
    pop %in% c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB") ~ "BC",
    pop %in% c("RS","SAY") ~ "Marine",
    pop %in% c("PACH","FRED","SC","CH") ~ "Recent",
    TRUE ~ "Other"
  ))

p12 <- ggplot(PC, aes(mitoPC1, mitoPC2, color=region)) +
  geom_point(size=3, alpha=0.95) +
  ggrepel::geom_text_repel(aes(label=pop), size=3, max.overlaps=Inf, seed=1) +
  labs(
    x=sprintf("mitoPC1 (%.1f%%)", expl[1]),
    y=sprintf("mitoPC2 (%.1f%%)", expl[2]),
    title="Mitochondrial PCoA (with AMO)"
  ) +
  theme_bw(base_size=12) +
  theme(
    panel.grid.major = element_line(linewidth=0.3),
    panel.grid.minor = element_line(linewidth=0.2),
    legend.position="right"
  )

png(file.path(out_dir, "mitoPC_withAMO_PC1_PC2.theme_bw.png"), width=2200, height=1600, res=300)
print(p12)
dev.off()

pdf(file.path(out_dir, "mitoPC_withAMO_PC1_PC2.theme_bw.pdf"), width=7.2, height=5.2)
print(p12)
dev.off()

# explained variance 表（补充材料/记录）
exp_out <- file.path(out_dir, "mitoPC_withAMO_explained.tsv")
fwrite(data.table(PC=paste0("PC",1:5), explained_percent=expl[1:5]), exp_out, sep="\t")
cat("[write] ", exp_out, "\n")

# mitoPC 坐标表（后续 merge covariates 用）
pc_out <- file.path(out_dir, "mitoPC_withAMO.tsv")
fwrite(as.data.table(PC), pc_out, sep="\t")
cat("[write] ", pc_out, "\n")



#af～ mtpc12345
cat > /mnt/spareHD_2/nu_287/_assoc72_subunit/11_AF_vs_mitoPC1to5_plus_treePC12_withAMO.R <<'RS'
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

AF_FILE  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
PC_FILE  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_withAMO.tsv"
COV_FILE <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_withAMO"
MIN_DEPTH <- 10L

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

OUT_LM  <- file.path(OUTDIR, "AF_LM_perSNP_mitoPC1to5_plus_treePC12.withAMO.tsv.gz")
OUT_SUM <- file.path(OUTDIR, "AF_LM_summary_by_mitoPC.withAMO.tsv")

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

fast_lm_1x2t <- function(y, x, t1, t2){
  ok <- is.finite(y) & is.finite(x) & is.finite(t1) & is.finite(t2)
  y <- y[ok]; x <- x[ok]; t1 <- t1[ok]; t2 <- t2[ok]
  n <- length(y)
  if(n < 6) return(NULL)

  X <- cbind(1, x, t1, t2)
  fit <- lm.fit(X, y)

  p <- ncol(X)
  df_res <- n - p
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

# ----------------------------
# Read AF
# ----------------------------
cat("[read] AF: ", AF_FILE, "\n", sep="")
AF <- fread(AF_FILE)
need_af <- c("chr","pos","gene","pop","af","depth")
miss_af <- setdiff(need_af, names(AF))
if(length(miss_af)>0) stop("AF missing columns: ", paste(miss_af, collapse=", "))

AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

# ★关键：同一 pop 多个 sample 汇总为每 snp×pop 一行（depth 加权）
AF_POP <- AF[, .(
  af = ifelse(sum(depth, na.rm=TRUE) > 0,
              sum(af * depth, na.rm=TRUE) / sum(depth, na.rm=TRUE),
              mean(af, na.rm=TRUE)),
  depth = sum(depth, na.rm=TRUE)
), by=.(snp, chr, pos, gene, pop)]

cat("[info] AF_POP rows: ", nrow(AF_POP), "\n", sep="")
cat("[info] unique pops (AF_POP): ", uniqueN(AF_POP$pop), "\n", sep="")

# ----------------------------
# Read mitoPC
# ----------------------------
cat("[read] PC: ", PC_FILE, "\n", sep="")
PC <- fread(PC_FILE)
pc_pop_col <- names(PC)[grep("^pop$|population|pop_id|site|id", tolower(names(PC)))[1]]
if(is.na(pc_pop_col)) stop("Cannot find pop column in PC file.")
setnames(PC, pc_pop_col, "pop")
PC[, pop := normalize_pop(pop)]

pc_cols <- paste0("mitoPC", 1:5)
miss_pc <- setdiff(c("pop", pc_cols), names(PC))
if(length(miss_pc)>0) stop("PC file missing columns: ", paste(miss_pc, collapse=", "))

PC <- unique(PC[, c("pop", pc_cols), with=FALSE], by="pop")

# ----------------------------
# Read treePC
# ----------------------------
cat("[read] COV(treePC): ", COV_FILE, "\n", sep="")
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]

need_cov <- c("pop","treePC1","treePC2")
miss_cov <- setdiff(need_cov, names(COV))
if(length(miss_cov)>0) stop("COV missing columns: ", paste(miss_cov, collapse=", "))

COV <- unique(COV[, ..need_cov], by="pop")

# ----------------------------
# Merge predictors
# ----------------------------
PRED <- merge(PC, COV, by="pop", all=FALSE)
DT <- merge(
  AF_POP[, .(snp, chr, pos, gene, pop, af)],
  PRED,
  by="pop",
  all=FALSE
)

cat("[info] rows after merge: ", nrow(DT), "\n", sep="")
cat("[info] unique pops in merged: ", uniqueN(DT$pop), "\n", sep="")
cat("[info] unique SNPs: ", uniqueN(DT$snp), "\n", sep="")

ncheck <- DT[, .(n_pop=uniqueN(pop), n_row=.N), by=.(snp)]
cat("[check] n_pop summary:\n")
print(summary(ncheck$n_pop))

setkey(DT, snp)

# ----------------------------
# Per-SNP regressions
# ----------------------------
res_list <- vector("list", length(pc_cols))
names(res_list) <- pc_cols

for(pc in pc_cols){
  cat("[lm] ", pc, "\n", sep="")
  tmp <- DT[, {
    fit <- fast_lm_1x2t(af, get(pc), treePC1, treePC2)
    if(is.null(fit)) NULL else .(
      chr = chr[1],
      pos = pos[1],
      gene= gene[1],
      mitoPC = pc,
      beta_mitoPC = fit$beta,
      se_mitoPC   = fit$se,
      t_mitoPC    = fit$t,
      p_mitoPC    = fit$p,
      R2          = fit$R2,
      adjR2       = fit$adjR2,
      n           = fit$n
    )
  }, by=.(snp)]
  res_list[[pc]] <- tmp
  cat("  -> rows:", nrow(tmp), "\n")
}

RES <- rbindlist(res_list, use.names=TRUE, fill=TRUE)
setorder(RES, mitoPC, p_mitoPC)
fwrite(RES, OUT_LM, sep="\t", compress="gzip")
cat("[write] ", OUT_LM, "\n", sep="")

# ----------------------------
# Summary per mitoPC
# ----------------------------
RES2 <- RES[is.finite(p_mitoPC)]
SUM <- RES2[, .(
  n_snps = .N,
  min_p  = min(p_mitoPC, na.rm=TRUE),
  med_p  = median(p_mitoPC, na.rm=TRUE),
  prop_p05 = mean(p_mitoPC < 0.05, na.rm=TRUE),
  prop_p10 = mean(p_mitoPC < 0.10, na.rm=TRUE),
  mean_beta = mean(beta_mitoPC, na.rm=TRUE),
  prop_beta_pos = mean(beta_mitoPC > 0, na.rm=TRUE),
  median_n = median(n, na.rm=TRUE)
), by=.(mitoPC)]

setorder(SUM, -prop_p05, med_p)
fwrite(SUM, OUT_SUM, sep="\t")
cat("[write] ", OUT_SUM, "\n", sep="")

cat("\n=== summary by mitoPC ===\n")
print(SUM)

cat("\n[OK] done. Output dir:\n  ", OUTDIR, "\n", sep="")
RS

Rscript /mnt/spareHD_2/nu_287/_assoc72_subunit/11_AF_vs_mitoPC1to5_plus_treePC12_withAMO.R


#lamda
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

in_file <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_withAMO/AF_LM_perSNP_mitoPC1to5_plus_treePC12.withAMO.tsv.gz"
out_dir <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_withAMO/QQ_lambda"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

DT <- fread(cmd=paste("zcat", shQuote(in_file)))
DT <- DT[is.finite(p_mitoPC) & p_mitoPC>0 & p_mitoPC<=1]

# lambda_GC = median(chi^2)/median_chi^2_null
lambda_gc <- function(p){
  chisq <- qchisq(1 - p, df=1)
  median(chisq, na.rm=TRUE) / qchisq(0.5, df=1)  # qchisq(0.5,1)=0.4549...
}

SUM <- DT[, .(
  n = .N,
  lambda = lambda_gc(p_mitoPC),
  min_p = min(p_mitoPC),
  med_p = median(p_mitoPC),
  prop_p05 = mean(p_mitoPC < 0.05)
), by=.(mitoPC)]

setorder(SUM, mitoPC)
fwrite(SUM, file.path(out_dir, "lambda_summary.withAMO.tsv"), sep="\t")
cat("[write] ", file.path(out_dir, "lambda_summary.withAMO.tsv"), "\n", sep="")
print(SUM)

# QQ plot per mitoPC
qqplot_dt <- DT[, {
  p <- sort(p_mitoPC)
  n <- length(p)
  .(exp = -log10(ppoints(n)),
    obs = -log10(p))
}, by=.(mitoPC)]

p <- ggplot(qqplot_dt, aes(exp, obs)) +
  geom_point(alpha=0.35, size=0.6) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~mitoPC, ncol=3) +
  theme_classic(base_size=13) +
  labs(x="Expected -log10(p)", y="Observed -log10(p)",
       title="QQ plots: AF ~ mitoPCk + treePC1 + treePC2 (with AMO)")

ggsave(file.path(out_dir, "QQ_AF_mitoPC1to5.withAMO.png"), p, width=10.5, height=6.5, dpi=300)
cat("[write] ", file.path(out_dir, "QQ_AF_mitoPC1to5.withAMO.png"), "\n", sep="")
summary(m2)
[1] 0.01481018

Call:
lm(formula = mitoPC1 ~ treePC1 + treePC2 + treePC3, data = df)

Residuals:
       Min         1Q     Median         3Q        Max 
-0.0048864 -0.0001780  0.0002680  0.0004720  0.0009656 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  5.248e-21  2.150e-04   0.000    1.000
treePC1      5.004e-04  1.243e-03   0.403    0.691
treePC2     -6.900e-04  1.946e-03  -0.355    0.726
treePC3     -5.684e-04  2.362e-03  -0.241    0.812

Residual standard error: 0.001117 on 23 degrees of freedom
Multiple R-squared:  0.01481,	Adjusted R-squared:  -0.1137 
F-statistic: 0.1153 on 3 and 23 DF,  p-value: 0.9503

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

in_file  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_withAMO/AF_LM_perSNP_mitoPC1to5_plus_treePC12.withAMO.tsv.gz"
out_dir  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_withAMO/BH"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

out_file <- file.path(out_dir, "AF_LM_perSNP_mitoPC1to5_plus_treePC12.withAMO.BH.tsv.gz")
out_sum  <- file.path(out_dir, "BH_summary_by_mitoPC.withAMO.tsv")

DT <- fread(cmd=paste("zcat", shQuote(in_file)))
DT <- DT[is.finite(p_mitoPC) & p_mitoPC>0 & p_mitoPC<=1]

# BH within each mitoPC group (each is its own family of tests)
DT[, BH_mitoPC := p.adjust(p_mitoPC, method="BH"), by=.(mitoPC)]

fwrite(DT, out_file, sep="\t", compress="gzip")
cat("[write] ", out_file, "\n", sep="")

SUM <- DT[, .(
  n = .N,
  n_p05 = sum(p_mitoPC < 0.05),
  prop_p05 = mean(p_mitoPC < 0.05),
  n_BH05 = sum(BH_mitoPC < 0.05),
  prop_BH05 = mean(BH_mitoPC < 0.05),
  min_p = min(p_mitoPC),
  min_BH = min(BH_mitoPC),
  med_p = median(p_mitoPC),
  med_BH = median(BH_mitoPC)
), by=.(mitoPC)]

setorder(SUM, mitoPC)
fwrite(SUM, out_sum, sep="\t")
cat("[write] ", out_sum, "\n", sep="")
print(SUM)
Rscript /mnt/spareHD_2/nu_287/_assoc72_subunit/13_add_BH_and_summary_AF_mitoPC1to5.withAMO.R
[write] /mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_withAMO/BH/AF_LM_perSNP_mitoPC1to5_plus_treePC12.withAMO.BH.tsv.gz
[write] /mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_withAMO/BH/BH_summary_by_mitoPC.withAMO.tsv
    mitoPC     n n_p05   prop_p05 n_BH05   prop_BH05        min_p       min_BH
    <char> <int> <int>      <num>  <int>       <num>        <num>        <num>
1: mitoPC1 42759 21528 0.50347295  21510 0.503051989 3.692522e-12 2.214020e-11
2: mitoPC2 42759   996 0.02329334      0 0.000000000 1.448947e-03 7.661555e-01
3: mitoPC3 42759   902 0.02109497      0 0.000000000 4.853841e-04 7.079290e-01
4: mitoPC4 42759   663 0.01550551    487 0.011389415 1.871532e-15 3.177762e-12
5: mitoPC5 42759  1093 0.02556187    426 0.009962815 2.384926e-09 1.019771e-04
          med_p       med_BH
          <num>        <num>
1: 0.0004580765 0.0009161315
2: 0.5289571281 0.7661555449
3: 0.5586081046 0.7079289847
4: 0.7383347058 0.9981330849
5: 0.9443495696 0.9858696757

#af～ mt cluster
cat /mnt/spareHD_2/nu_287/_assoc72_subunit/13_AF_vs_mtCluster_plus_treePC12_withAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# withAMO | per-SNP regression:
#   af ~ mtCluster + treePC1 + treePC2
# Test mtCluster term by nested-model F-test:
#   reduced: af ~ treePC1 + treePC2
#   full:    af ~ mtCluster + treePC1 + treePC2
# Outputs include BH(q) across SNPs.
# ============================================================

AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR    <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mtCluster_withAMO"
MIN_DEPTH <- 10L

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

OUT_LM  <- file.path(OUTDIR, "AF_LM_perSNP_mtCluster_plus_treePC12.withAMO.tsv.gz")
OUT_SUM <- file.path(OUTDIR, "AF_LM_summary_mtCluster.withAMO.tsv")

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

lambda_gc <- function(p){
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if(length(p) < 10) return(NA_real_)
  chisq <- qchisq(1 - p, df=1)
  median(chisq, na.rm=TRUE) / qchisq(0.5, df=1)
}

lm_cluster_ftest <- function(dt){
  dt <- dt[is.finite(af) & is.finite(treePC1) & is.finite(treePC2) & !is.na(mtCluster)]
  if(nrow(dt) < 6) return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))
  if(length(unique(dt$mtCluster)) < 2) return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))

  fit_full <- try(lm(af ~ mtCluster + treePC1 + treePC2, data=dt), silent=TRUE)
  fit_red  <- try(lm(af ~            treePC1 + treePC2, data=dt), silent=TRUE)
  if(inherits(fit_full, "try-error") || inherits(fit_red, "try-error")){
    return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))
  }

  a <- try(anova(fit_red, fit_full), silent=TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2){
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

# ----------------------------
# Read AF (withAMO)
# ----------------------------
cat("[read] AF: ", AF_FILE, "\n", sep="")
AF <- fread(AF_FILE)
need_af <- c("chr","pos","gene","pop","af","depth")
miss_af <- setdiff(need_af, names(AF))
if(length(miss_af)>0) stop("AF missing columns: ", paste(miss_af, collapse=", "))

AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[info] AF rows after filters: ", nrow(AF), "\n", sep="")
cat("[info] unique pops: ", uniqueN(AF$pop), "\n", sep="")

# ----------------------------
# Read mtCluster manual (must include AMO)
# ----------------------------
cat("[read] mtCluster: ", CL_FILE, "\n", sep="")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
CL <- CL[!(is.na(pop) | pop=="")]
need_cl <- c("pop","mtCluster")
miss_cl <- setdiff(need_cl, names(CL))
if(length(miss_cl)>0) stop("Cluster file missing columns: ", paste(miss_cl, collapse=", "))

CL[, pop := normalize_pop(pop)]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

cat("[info] pops per cluster:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# ----------------------------
# Read treePC (withAMO)
# ----------------------------
cat("[read] treePC: ", COV_FILE, "\n", sep="")
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
need_cov <- c("pop","treePC1","treePC2")
miss_cov <- setdiff(need_cov, names(COV))
if(length(miss_cov)>0) stop("COV missing columns: ", paste(miss_cov, collapse=", "))

COV <- unique(COV[, ..need_cov], by="pop")

# ----------------------------
# Merge predictors into AF
# ----------------------------
PRED <- merge(COV, CL, by="pop", all=FALSE)

DT <- merge(
  AF[, .(snp, chr, pos, gene, pop, af)],
  PRED,
  by="pop",
  all=FALSE
)

cat("[info] rows after merge: ", nrow(DT), "\n", sep="")
cat("[info] unique SNPs: ", uniqueN(DT$snp), "\n", sep="")
cat("[info] unique pops in merged: ", uniqueN(DT$pop), "\n", sep="")

setkey(DT, snp)

# ----------------------------
# Per-SNP LM: F-test for mtCluster term
# ----------------------------
cat("[lm] per SNP F-test for mtCluster...\n")
LM <- DT[, {
  out <- lm_cluster_ftest(.SD)
  if(!is.finite(out$p)) NULL else .(
    chr=chr[1], pos=pos[1], gene=gene[1],
    n=out$n, df1=out$df1, df2=out$df2, F_cluster=out$F, p_cluster=out$p
  )
}, by=.(snp)]

# BH across SNPs (one family)
LM[, q_cluster := p.adjust(p_cluster, method="BH")]

fwrite(LM, OUT_LM, sep="\t", compress="gzip")
cat("[write] ", OUT_LM, "\n", sep="")

# ----------------------------
# Summary + lambda
# ----------------------------
p <- LM$p_cluster
q <- LM$q_cluster

SUM <- data.table(
  n_snps = nrow(LM),
  min_p = suppressWarnings(min(p, na.rm=TRUE)),
  med_p = suppressWarnings(median(p, na.rm=TRUE)),
  prop_p05 = mean(p < 0.05, na.rm=TRUE),
  prop_p10 = mean(p < 0.10, na.rm=TRUE),
  n_q05 = sum(q < 0.05, na.rm=TRUE),
  prop_q05 = mean(q < 0.05, na.rm=TRUE),
  lambda_gc = lambda_gc(p)
)

fwrite(SUM, OUT_SUM, sep="\t")
cat("[write] ", OUT_SUM, "\n", sep="")
print(SUM)

cat("\n[OK] done. Output dir:\n  ", OUTDIR, "\n", sep="")



cat /mnt/spareHD_2/nu_287/_assoc72_subunit/16_geneLevel_AF_vs_mtCluster_plus_treePC12_withAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/geneAF_vs_mtCluster_withAMO"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

MIN_DEPTH <- 10L
MIN_SNP_GENE <- 20L

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

lm_cluster_ftest_gene <- function(d){
  d <- d[is.finite(geneMeanAF) & is.finite(treePC1) & is.finite(treePC2) & !is.na(mtCluster)]
  if(uniqueN(d$pop) < 10) return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=uniqueN(d$pop)))
  if(length(unique(d$mtCluster)) < 2) return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=uniqueN(d$pop)))

  fit_full <- try(lm(geneMeanAF ~ mtCluster + treePC1 + treePC2, data=d), silent=TRUE)
  fit_red  <- try(lm(geneMeanAF ~            treePC1 + treePC2, data=d), silent=TRUE)
  if(inherits(fit_full, "try-error") || inherits(fit_red, "try-error")){
    return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=uniqueN(d$pop)))
  }

  a <- try(anova(fit_red, fit_full), silent=TRUE)
  if(inherits(a, "try-error") || nrow(a) < 2){
    return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=uniqueN(d$pop)))
  }

  list(
    p   = as.numeric(a$`Pr(>F)`[2]),
    F   = as.numeric(a$F[2]),
    df1 = as.numeric(a$Df[2]),
    df2 = as.numeric(a$Res.Df[2]),
    n   = uniqueN(d$pop)
  )
}

# ----------------------------
# Read AF and compute geneMeanAF per pop (withAMO)
# ----------------------------
cat("[read] AF\n")
AF <- fread(AF_FILE)
stopifnot(all(c("gene","pop","af","depth") %in% names(AF)))

AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= MIN_DEPTH]

cat("[summarize] gene mean AF per pop\n")
GAF <- AF[, .(
  n_snps = .N,
  geneMeanAF = mean(af, na.rm=TRUE)
), by=.(gene, pop)]

G_tot <- GAF[, .(n_snps_total = sum(n_snps)), by=gene]
GAF <- merge(GAF, G_tot, by="gene", all.x=TRUE)
GAF <- GAF[n_snps_total >= MIN_SNP_GENE]

cat("[info] genes kept: ", uniqueN(GAF$gene), "\n", sep="")

# ----------------------------
# Read mtCluster manual (must include AMO)
# ----------------------------
cat("[read] mtCluster\n")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
CL <- CL[!(is.na(pop) | pop=="")]
stopifnot(all(c("pop","mtCluster") %in% names(CL)))
CL[, pop := normalize_pop(pop)]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

cat("[info] pops per cluster:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# ----------------------------
# Read treePC (withAMO)
# ----------------------------
cat("[read] treePC\n")
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

# ----------------------------
# Merge predictors into gene-level table
# ----------------------------
PRED <- merge(COV, CL, by="pop", all=FALSE)
DT <- merge(GAF, PRED, by="pop", all=FALSE)

cat("[info] pops used: ", uniqueN(DT$pop), "\n", sep="")
cat("[info] genes used: ", uniqueN(DT$gene), "\n", sep="")

# ----------------------------
# Per-gene F-test for mtCluster term + BH
# ----------------------------
cat("[lm] per gene F-test for mtCluster\n")
RES <- DT[, {
  out <- lm_cluster_ftest_gene(.SD)
  if(!is.finite(out$p)) NULL else .(
    n_pop = out$n,
    df1 = out$df1, df2 = out$df2,
    F_cluster = out$F,
    p_cluster = out$p
  )
}, by=.(gene)]

RES[, q_cluster := p.adjust(p_cluster, method="BH")]
RES <- merge(RES, unique(GAF[, .(gene, n_snps_total)], by="gene"), by="gene", all.x=TRUE)

out_full <- file.path(OUTDIR, "geneLevel_AF_LM_byGene_mtCluster_plus_treePC12.withAMO.tsv.gz")
out_sum  <- file.path(OUTDIR, "geneLevel_AF_summary_mtCluster.withAMO.tsv")

fwrite(RES[order(q_cluster, p_cluster)], out_full, sep="\t", compress="gzip")
cat("[write] ", out_full, "\n", sep="")

SUM <- RES[, .(
  n_genes = .N,
  min_p = min(p_cluster, na.rm=TRUE),
  med_p = median(p_cluster, na.rm=TRUE),
  n_q05 = sum(q_cluster < 0.05, na.rm=TRUE),
  n_q10 = sum(q_cluster < 0.10, na.rm=TRUE)
)]

fwrite(SUM, out_sum, sep="\t")
cat("[write] ", out_sum, "\n", sep="")
print(SUM)

cat("\n=== top 10 genes by p_cluster ===\n")
print(head(RES[order(p_cluster)], 10))

cat("\n[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")
















#/mnt/spareHD_2/nu_287/_assoc72_subunit/20_partialRDA_OXPHOS_AF_withAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

# =========================
# withAMO | partial RDA
#   Y_hel ~ mitoPC1..5 + Condition(treePC1 + treePC2)
#   where Y = pop × SNP AF matrix (complete cases)
#   Y_hel = Hellinger(Y) via decostand(., "hellinger")
# =========================

AF_FILE  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
PC_FILE  <- "/mnt/spareHD_2/nu_287/meta/mtPC_from_mtTree_noDloop_withAMO.tsv.gz"   # 你现在有的 withAMO mtPC
COV_FILE <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/rda_mitoPC_withAMO"
INDIR    <- file.path(OUTDIR, "rda_inputs")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
dir.create(INDIR,  showWarnings=FALSE, recursive=TRUE)

MIN_DEPTH <- 10L

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

cat("[read] AF\n")
AF <- fread(AF_FILE)
stopifnot(all(c("chr","pos","gene","pop","af","depth") %in% names(AF)))
AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[read] mitoPC\n")
PC <- fread(PC_FILE)
# 兼容 pop 列名
pc_pop_col <- names(PC)[grep("^pop$|population|pop_id|site|id", tolower(names(PC)))[1]]
if(is.na(pc_pop_col)) stop("Cannot find pop column in PC file.")
setnames(PC, pc_pop_col, "pop")
PC[, pop := normalize_pop(pop)]
pc_cols <- paste0("mtPC", 1:5)
# 兼容你文件列名是 mtPC1..5 而不是 mitoPC1..5
if(!all(pc_cols %in% names(PC))){
  # 尝试 mitoPC1..5
  alt <- paste0("mitoPC", 1:5)
  if(all(alt %in% names(PC))){
    setnames(PC, alt, pc_cols)
  } else {
    stop("PC file missing mtPC1-5 (or mitoPC1-5).")
  }
}
PC <- unique(PC[, c("pop", pc_cols), with=FALSE], by="pop")

cat("[read] treePC\n")
COV <- fread(COV_FILE)
stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))
COV[, pop := normalize_pop(pop)]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

# 预测变量合并（withAMO：不删 AMO）
META <- merge(PC, COV, by="pop", all=FALSE)
setorder(META, pop)

cat("[info] pops: ", nrow(META), "\n", sep="")
cat("[info] pops: ", paste(META$pop, collapse=","), "\n", sep="")

# -------------------------
# Build Y matrix: pop × snp
# -------------------------
cat("[build] AF matrix\n")
DT <- merge(AF[, .(pop, snp, af)], META[, .(pop)], by="pop", all=FALSE)
# dcast: rows=pop, cols=snp
Ydt <- dcast(DT, pop ~ snp, value.var="af")
Ypop <- Ydt$pop
Y <- as.matrix(Ydt[, -"pop"])

cat("[info] SNPs total: ", ncol(Y), "\n", sep="")
# complete SNPs across all pops
ok <- colSums(!is.na(Y)) == nrow(Y)
Y <- Y[, ok, drop=FALSE]
cat("[info] SNPs complete: ", ncol(Y), "\n", sep="")

# 保存 inputs（可复用，不用每次重建）
saveRDS(list(Y=Y, pops=Ypop), file.path(INDIR, "Y_AF_matrix.withAMO.rds"))
saveRDS(META, file.path(INDIR, "meta_predictors.withAMO.rds"))

# Hellinger
Y_hel <- decostand(Y, method="hellinger")

# 组装 meta2（按 Y 的 pop 顺序）
meta2 <- META[match(Ypop, pop)]
stopifnot(all(meta2$pop == Ypop))

cat("[rda] partial RDA: Y_hel ~ mtPC1-5 + Condition(treePC1-2)\n")
form <- as.formula("Y_hel ~ mtPC1 + mtPC2 + mtPC3 + mtPC4 + mtPC5 + Condition(treePC1 + treePC2)")
fit <- rda(form, data=meta2)

# permutation tests
a_global <- anova.cca(fit, permutations=999)
a_terms  <- anova.cca(fit, by="term", permutations=999)
r2 <- RsquareAdj(fit)

# write results
OUT_TXT <- file.path(OUTDIR, "partialRDA_mitoPC_withAMO.results.txt")
OUT_RDS <- file.path(OUTDIR, "partialRDA_mitoPC_withAMO.fit.rds")

txt <- c(
  "partial RDA: Y_hel ~ mtPC1..5 + Condition(treePC1 + treePC2)",
  "",
  paste0("AF_FILE: ", AF_FILE),
  paste0("PC_FILE: ", PC_FILE),
  paste0("COV_FILE:", COV_FILE),
  paste0("MIN_DEPTH: ", MIN_DEPTH),
  "",
  paste0("Pops: ", nrow(meta2)),
  paste0("SNPs complete: ", ncol(Y)),
  "",
  "=== anova(global) ===",
  capture.output(a_global),
  "",
  "=== anova(by term) ===",
  capture.output(a_terms),
  "",
  "=== RsquareAdj ===",
  capture.output(r2)
)

writeLines(txt, OUT_TXT)
saveRDS(fit, OUT_RDS)

cat("[write] ", OUT_TXT, "\n", sep="")
cat("[write] ", OUT_RDS, "\n", sep="")
cat("[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")

#21_partialRDA_mtCluster_withAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

# =========================
# withAMO | partial RDA (factor)
#   Y_hel ~ mtCluster + Condition(treePC1 + treePC2)
# =========================

AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/rda_mtCluster_withAMO"
INDIR    <- file.path(OUTDIR, "rda_inputs")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
dir.create(INDIR,  showWarnings=FALSE, recursive=TRUE)

MIN_DEPTH <- 10L

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

cat("[read] AF\n")
AF <- fread(AF_FILE)
stopifnot(all(c("chr","pos","gene","pop","af","depth") %in% names(AF)))
AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[read] mtCluster\n")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
stopifnot(all(c("pop","mtCluster") %in% names(CL)))
CL[, pop := normalize_pop(pop)]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

cat("[read] treePC\n")
COV <- fread(COV_FILE)
stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))
COV[, pop := normalize_pop(pop)]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

META <- merge(COV, CL, by="pop", all=FALSE)
setorder(META, pop)

cat("[info] pops: ", nrow(META), "\n", sep="")
cat("[info] pops per cluster:\n")
print(as.data.table(table(META$mtCluster))[order(-N)])

# Build Y
cat("[build] Y AF matrix\n")
DT <- merge(AF[, .(pop, snp, af)], META[, .(pop)], by="pop", all=FALSE)
Ydt <- dcast(DT, pop ~ snp, value.var="af")
Ypop <- Ydt$pop
Y <- as.matrix(Ydt[, -"pop"])
cat("[info] SNPs total: ", ncol(Y), "\n", sep="")

ok <- colSums(!is.na(Y)) == nrow(Y)
Y <- Y[, ok, drop=FALSE]
cat("[info] SNPs complete: ", ncol(Y), "\n", sep="")

saveRDS(list(Y=Y, pops=Ypop), file.path(INDIR, "Y_AF_matrix.withAMO.rds"))
saveRDS(META, file.path(INDIR, "meta_predictors.withAMO.rds"))

Y_hel <- decostand(Y, method="hellinger")
meta2 <- META[match(Ypop, pop)]
stopifnot(all(meta2$pop == Ypop))

cat("[rda] partial RDA: Y_hel ~ mtCluster + Condition(treePC1-2)\n")
fit <- rda(Y_hel ~ mtCluster + Condition(treePC1 + treePC2), data=meta2)

a_global <- anova.cca(fit, permutations=999)
a_terms  <- anova.cca(fit, by="term", permutations=999)
r2 <- RsquareAdj(fit)

OUT_TXT <- file.path(OUTDIR, "partialRDA_mtCluster_withAMO.results.txt")
OUT_RDS <- file.path(OUTDIR, "partialRDA_mtCluster_withAMO.fit.rds")

txt <- c(
  "partial RDA: Y_hel ~ mtCluster + Condition(treePC1 + treePC2)",
  "",
  paste0("AF_FILE: ", AF_FILE),
  paste0("CL_FILE: ", CL_FILE),
  paste0("COV_FILE:", COV_FILE),
  paste0("MIN_DEPTH: ", MIN_DEPTH),
  "",
  paste0("Pops: ", nrow(meta2)),
  paste0("SNPs complete: ", ncol(Y)),
  "",
  "=== anova(global) ===",
  capture.output(a_global),
  "",
  "=== anova(by term) ===",
  capture.output(a_terms),
  "",
  "=== RsquareAdj ===",
  capture.output(r2)
)

writeLines(txt, OUT_TXT)
saveRDS(fit, OUT_RDS)

cat("[write] ", OUT_TXT, "\n", sep="")
cat("[write] ", OUT_RDS, "\n", sep="")
cat("[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")

#22_partialdbRDA_mtCluster_withAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

# =========================
# withAMO | partial dbRDA
#   D = vegdist(Hellinger(Y), method="bray")
#   capscale: D ~ mtCluster + Condition(treePC1 + treePC2)
# =========================

AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/dbrda_mtCluster_withAMO"
INDIR    <- file.path(OUTDIR, "rda_inputs")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
dir.create(INDIR,  showWarnings=FALSE, recursive=TRUE)

MIN_DEPTH <- 10L
DIST_METHOD <- "bray"

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

cat("[read] AF\n")
AF <- fread(AF_FILE)
stopifnot(all(c("chr","pos","gene","pop","af","depth") %in% names(AF)))
AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[read] mtCluster\n")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
stopifnot(all(c("pop","mtCluster") %in% names(CL)))
CL[, pop := normalize_pop(pop)]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

cat("[read] treePC\n")
COV <- fread(COV_FILE)
stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))
COV[, pop := normalize_pop(pop)]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

META <- merge(COV, CL, by="pop", all=FALSE)
setorder(META, pop)

cat("[info] pops: ", nrow(META), "\n", sep="")
cat("[info] pops per cluster:\n")
print(as.data.table(table(META$mtCluster))[order(-N)])

# Build Y
cat("[build] Y AF matrix\n")
DT <- merge(AF[, .(pop, snp, af)], META[, .(pop)], by="pop", all=FALSE)
Ydt <- dcast(DT, pop ~ snp, value.var="af")
Ypop <- Ydt$pop
Y <- as.matrix(Ydt[, -"pop"])

cat("[info] SNPs total: ", ncol(Y), "\n", sep="")
ok <- colSums(!is.na(Y)) == nrow(Y)
Y <- Y[, ok, drop=FALSE]
cat("[info] SNPs complete: ", ncol(Y), "\n", sep="")

saveRDS(list(Y=Y, pops=Ypop), file.path(INDIR, "Y_AF_matrix.withAMO.rds"))
saveRDS(META, file.path(INDIR, "meta_predictors.withAMO.rds"))

# Hellinger + distance
Y_hel <- decostand(Y, method="hellinger")
D <- vegdist(Y_hel, method=DIST_METHOD)

meta2 <- META[match(Ypop, pop)]
stopifnot(all(meta2$pop == Ypop))

cat("[dbrda] capscale: D(Y) ~ mtCluster + Condition(treePC1-2)\n")
fit <- capscale(D ~ mtCluster + Condition(treePC1 + treePC2), data=meta2)

a_global <- anova.cca(fit, permutations=999)
a_terms  <- anova.cca(fit, by="term", permutations=999)
r2 <- RsquareAdj(fit)

OUT_TXT <- file.path(OUTDIR, "partialdbRDA_mtCluster_withAMO.results.txt")
OUT_RDS <- file.path(OUTDIR, "partialdbRDA_mtCluster_withAMO.fit.rds")

txt <- c(
  "partial dbRDA (capscale): vegdist(Hellinger(Y)) ~ mtCluster + Condition(treePC1 + treePC2)",
  "",
  paste0("DIST_METHOD: ", DIST_METHOD),
  paste0("AF_FILE: ", AF_FILE),
  paste0("CL_FILE: ", CL_FILE),
  paste0("COV_FILE:", COV_FILE),
  paste0("MIN_DEPTH: ", MIN_DEPTH),
  "",
  paste0("Pops: ", nrow(meta2)),
  paste0("SNPs complete: ", ncol(Y)),
  "",
  "=== anova(global) ===",
  capture.output(a_global),
  "",
  "=== anova(by term) ===",
  capture.output(a_terms),
  "",
  "=== RsquareAdj ===",
  capture.output(r2)
)

writeLines(txt, OUT_TXT)
saveRDS(fit, OUT_RDS)

cat("[write] ", OUT_TXT, "\n", sep="")
cat("[write] ", OUT_RDS, "\n", sep="")
cat("[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")





#delta af～ mtpc12345
#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })

af_file  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
pc_file  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_withAMO.tsv"      # <-- CHANGED
cov_file <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

out_dir  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

out_delta <- file.path(out_dir, "deltaAF_long.withAMO.tsv.gz")
out_mito  <- file.path(out_dir, "LM_perSNP_mitoPCk_plus_treePC12.withAMO.tsv.gz")
out_tree  <- file.path(out_dir, "LM_perSNP_treeOnly.withAMO.tsv.gz")
out_log   <- file.path(out_dir, "LOG_summary.txt")

min_depth <- 10
AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")  # <-- INCLUDE AMO
AK_marine <- "RS"
BC_marine <- "SAY"
min_n_AK <- 5
min_n_BC <- 8  # BC now has 13 pops; threshold 8 still fine

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}
infer_region <- function(pop){
  fifelse(pop %in% c(AK_fresh, AK_marine), "AK",
          fifelse(pop %in% c(BC_fresh, BC_marine), "BC", NA_character_))
}

# ===== fast LM helpers (unchanged) =====
fast_lm_1x2t <- function(y, x, t1, t2){
  ok <- is.finite(y) & is.finite(x) & is.finite(t1) & is.finite(t2)
  y <- y[ok]; x <- x[ok]; t1 <- t1[ok]; t2 <- t2[ok]
  n <- length(y); if(n < 5) return(NULL)
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
  n <- length(y); if(n < 4) return(NULL)
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
AF <- fread(cmd=paste("zcat", shQuote(af_file)), sep=",", header=TRUE)
stopifnot(all(c("chr","pos","gene","pop","af","depth") %in% names(AF)))
AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= min_depth]
AF[, snp := paste(chr, pos, gene, sep=":")]

keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]
AF[, region := infer_region(pop)]
AF <- AF[!is.na(region)]

# ========= read mitoPC + treePC =========
PC <- fread(pc_file)
PC[, pop := normalize_pop(pop)]                       # <-- CHANGED
COV <- fread(cov_file)
COV[, pop := normalize_pop(pop)]                      # <-- CHANGED
stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))

pc_cols <- grep("^mitoPC[1-5]$", names(PC), value=TRUE)
stopifnot(length(pc_cols) >= 1)

PC  <- unique(PC[, c("pop", pc_cols), with=FALSE], by="pop")
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

# DEBUG: missing pops in predictors
pred_pops <- Reduce(intersect, list(PC$pop, COV$pop, unique(AF$pop)))
miss_in_PC  <- setdiff(unique(AF$pop), PC$pop)
miss_in_COV <- setdiff(unique(AF$pop), COV$pop)

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

SNP_ANN <- unique(AF[, .(snp, chr, pos, gene)])
DT <- merge(DT, SNP_ANN, by="snp", all.x=TRUE)

DT <- merge(DT, PRED, by="pop", all=FALSE)

DT[, n_in_snp := .N, by=.(region, snp)]
DT <- DT[(region=="AK" & n_in_snp >= min_n_AK) | (region=="BC" & n_in_snp >= min_n_BC)]
DT[, n_in_snp := NULL]

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

mito_res_list <- vector("list", length(pc_cols)); names(mito_res_list) <- pc_cols
for(pc in pc_cols){
  tmp <- DT[, {
    fit <- fast_lm_1x2t(deltaAF, get(pc), treePC1, treePC2)
    if(is.null(fit)) return(NULL)
    .(chr=chr[1], pos=pos[1], gene=gene[1],
      mitoPC=pc, beta_mitoPC=fit$beta, se_mitoPC=fit$se, t_mitoPC=fit$t, p_mitoPC=fit$p,
      R2=fit$R2, adjR2=fit$adjR2, n=fit$n)
  }, by=.(region, snp)]
  mito_res_list[[pc]] <- tmp
}
mito_res <- rbindlist(mito_res_list, use.names=TRUE, fill=TRUE)
setorder(mito_res, region, mitoPC, p_mitoPC)
fwrite(mito_res, out_mito, sep="\t", compress="gzip")

# ========= log =========
log_txt <- c(
  sprintf("[INFO] unique AF pops: %d (%s)", uniqueN(AF$pop), paste(sort(unique(AF$pop)), collapse=",")),
  sprintf("[INFO] PC pops: %d | COV pops: %d", nrow(PC), nrow(COV)),
  sprintf("[WARN] missing in PC: %s", ifelse(length(miss_in_PC)==0, "NONE", paste(miss_in_PC, collapse=","))),
  sprintf("[WARN] missing in COV: %s", ifelse(length(miss_in_COV)==0, "NONE", paste(miss_in_COV, collapse=","))),
  sprintf("[INFO] deltaAF observations: %d", nrow(DT)),
  sprintf("[INFO] mitoPC cols used: %s", paste(pc_cols, collapse=",")),
  sprintf("[INFO] tree-only rows: %d", nrow(tree_res)),
  sprintf("[INFO] mito rows: %d", nrow(mito_res)),
  "",
  "=== pops used after merge (DT) ===",
  capture.output(DT[, .N, by=.(region, pop)][order(region, -N)])
)
writeLines(log_txt, out_log)

cat("[OK] wrote:\n")
cat(" - ", out_delta, "\n", sep="")
cat(" - ", out_tree,  "\n", sep="")
cat(" - ", out_mito,  "\n", sep="")
cat(" - ", out_log,   "\n", sep="")



#delta af～ mtcluster
#/mnt/spareHD_2/nu_287/q2_parallelism/09_mtCluster_deltaAF_SNPlevel_LM_perm.withAMO.requireTree.R
MT_TREE=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile \
DELTA=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz \
OUTDIR=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_withAMO \
K=4 N_PERM=200 SEED=1 \
Rscript /mnt/spareHD_2/nu_287/q2_parallelism/09_mtCluster_deltaAF_SNPlevel_LM_perm.withAMO.requireTree.R
#!/usr/bin/env Rscript
# ============================================================
# 09_mtCluster_deltaAF_SNPlevel_LM_perm.withAMO.requireTree.R
#
# REQUIRE explicit MT_TREE (no auto-search).
# Hard-check that AMO exists in tree tips (normalized).
#
# Model per SNP within region:
#   deltaAF ~ mtCluster + treePC1 + treePC2
# Test mtCluster via nested-model F-test (anova(reduced, full)).
#
# Permutation:
#   permute mtCluster labels among POPs within each region,
#   recompute F for each SNP, empirical p_perm.
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
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz"
)

OUT_DIR <- Sys.getenv("OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_withAMO"
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

stop_msg_tree <- function(){
  cat("\n[FATAL] MT_TREE is required and must exist.\n")
  cat("Example:\n")
  cat("  MT_TREE=/work/.../aligned_mt_noDloop.withAMO.iqtree.treefile \\\n")
  cat("  DELTA=/.../deltaAF_long.withAMO.tsv.gz OUTDIR=/... N_PERM=200 K=4 \\\n")
  cat("  Rscript 09_mtCluster_deltaAF_SNPlevel_LM_perm.withAMO.requireTree.R\n\n")
  stop("MT_TREE not set / not found.")
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
DEL <- DEL[is.finite(deltaAF) & is.finite(treePC1) & is.finite(treePC2)]
cat("[info] rows after NA filter: ", nrow(DEL), "\n", sep="")
cat("[info] unique pops in DELTA: ", uniqueN(DEL$pop), "\n", sep="")

# ----------------------------
# Require MT_TREE and build clusters
# ----------------------------
TREE_FILE <- Sys.getenv("MT_TREE", unset = "")
if(!nzchar(TREE_FILE) || !file.exists(TREE_FILE)) stop_msg_tree()

cat("[read tree] ", TREE_FILE, "\n", sep="")
tr <- read.tree(TREE_FILE)
cat("[tree] ntip = ", length(tr$tip.label), "\n", sep="")
cat("[tree] first tips: ", paste(head(tr$tip.label, 10), collapse=", "), "\n", sep="")

# Hard check AMO exists in tree (normalized)
tip_norm <- normalize_pop(tr$tip.label)
if(!("AMO" %in% tip_norm)){
  cat("\n[FATAL] AMO not found in MT_TREE tips (normalized).\n")
  cat("TREE_FILE = ", TREE_FILE, "\n", sep="")
  cat("First 30 tips (normalized):\n")
  print(head(tip_norm, 30))
  stop("Not a withAMO treefile. Please provide a tree that includes AMO.")
}

# Distances -> clustering
D  <- cophenetic(tr)
hc <- hclust(as.dist(D), method="average")
cl <- cutree(hc, k=K_CLUST)

mt_cluster <- data.table(
  pop = normalize_pop(names(cl)),
  mt_cluster = as.integer(cl)
)
mt_cluster <- unique(mt_cluster, by="pop")

# Keep only pops that appear in DELTA
pops_in_DEL <- unique(DEL$pop)
mt_cluster  <- mt_cluster[pop %in% pops_in_DEL]

# Merge cluster into DEL
DEL <- merge(DEL, mt_cluster, by="pop", all.x=TRUE)

cat("\n[info] rows deltaAF: ", nrow(DEL), "\n", sep="")
cat("[info] pops missing mt_cluster rows: ", sum(is.na(DEL$mt_cluster)), "\n", sep="")
if(sum(is.na(DEL$mt_cluster)) > 0){
  cat("[WARN] pops without mt_cluster (unique):\n")
  print(unique(DEL[is.na(mt_cluster), .(pop)])[order(pop)])
}

# Drop missing cluster rows
DEL <- DEL[!is.na(mt_cluster)]
DEL[, cluster := factor(paste0("C", mt_cluster))]

# Save pop->cluster table
pop_region <- unique(DEL[, .(pop, region)], by="pop")
pop_cluster_out <- merge(pop_region, unique(DEL[, .(pop, mt_cluster)], by="pop"), by="pop", all.x=TRUE)
out_popcl <- file.path(OUT_DIR, "mtCluster_by_pop.withAMO.tsv")
fwrite(pop_cluster_out[order(region, mt_cluster, pop)], out_popcl, sep="\t")
cat("[write] ", out_popcl, "\n", sep="")

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
    out <- lm_cluster_Ftest(dt)

    cm <- dt[, .(mu=mean(deltaAF, na.rm=TRUE), n=.N), by=cluster][order(cluster)]
    cm_str <- paste0(cm$cluster, ":", sprintf("%.6g", cm$mu), "(n=", cm$n, ")", collapse=";")

    idx <- idx + 1L
    lm_list[[idx]] <- data.table(
      region=r, snp=s,
      chr=dt$chr[1], pos=dt$pos[1], gene=dt$gene[1],
      n=out$n, n_clusters=length(unique(dt$cluster)),
      F_cluster=out$F, df1=out$df1, df2=out$df2, p_cluster=out$p,
      cluster_means=cm_str
    )
  }
}

LM <- rbindlist(lm_list[seq_len(idx)], use.names=TRUE, fill=TRUE)

lm_out <- file.path(OUT_DIR, "LM_perSNP_mtCluster_plus_treePC12.withAMO.tsv.gz")
fwrite(LM, lm_out, sep="\t")
cat("\n[write] ", lm_out, "\n", sep="")

SUM <- LM[, .(
  n_snps = .N,
  prop_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
  prop_p10 = mean(p_cluster < 0.10, na.rm=TRUE),
  min_p = suppressWarnings(min(p_cluster, na.rm=TRUE)),
  med_p = suppressWarnings(median(p_cluster, na.rm=TRUE))
), by=region][order(region)]

out_sum <- file.path(OUT_DIR, "LM_cluster_summary_by_region.withAMO.tsv")
fwrite(SUM, out_sum, sep="\t")
cat("[write] ", out_sum, "\n\n", sep="")
cat("=== summary (cluster term) ===\n")
print(SUM)

# ----------------------------
# Permutation test
# ----------------------------
perm_out_file <- file.path(OUT_DIR, "Perm_perSNP_mtCluster_plus_treePC12.withAMO.tsv.gz")
perm_global_file <- file.path(OUT_DIR, "Perm_global_summary.withAMO.tsv")

if(N_PERM <= 0){
  cat("\n[perm] N_PERM=0 -> skip permutation.\n")
} else {

  cat("\n[perm] running permutations: N_PERM=", N_PERM, "\n", sep="")

  pop_cluster_dt <- unique(DEL[, .(pop, region, cluster)], by=c("pop","region"))

  OBS <- LM[, .(region, snp, F_obs=F_cluster)]
  setkey(OBS, region, snp)
  OBS[, ge_count := 0L]
  OBS[, n_perm_used := 0L]

  BASE <- DEL[, .(region, snp, pop, deltaAF, treePC1, treePC2)]
  setkey(BASE, region, snp)

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

        dt2 <- copy(dt)
        dt2[, cluster := cluster_perm]
        out <- lm_cluster_Ftest(dt2)
        if(!is.finite(out$F)) next

        f_obs <- OBS[list(r, s), F_obs]
        if(!is.finite(f_obs)) next

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

  G <- PERM[, .(
    n_snps=.N,
    prop_LM_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
    prop_perm_p05 = mean(p_perm < 0.05, na.rm=TRUE),
    med_perm_p = suppressWarnings(median(p_perm, na.rm=TRUE)),
    n_perm_used_med = suppressWarnings(median(n_perm_used, na.rm=TRUE))
  ), by=region][order(region)]
  fwrite(G, perm_global_file, sep="\t")
  cat("[write] ", perm_global_file, "\n", sep="")

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
       title="SNP-level association: ΔAF ~ mtCluster + treePC1 + treePC2 (withAMO; tree-required)") +
  theme_classic(base_size=14)
ggsave(file.path(OUT_DIR, "Fig_pHist_mtCluster_LM.withAMO.png"), p_hist, width=7.2, height=6.5, dpi=300)

LMq <- LM[is.finite(p_cluster) & p_cluster > 0 & p_cluster <= 1]
LMq[, obs := -log10(sort(p_cluster))]
LMq[, exp := -log10(ppoints(.N))]
p_qq <- ggplot(LMq, aes(exp, obs)) +
  geom_point(alpha=0.4, size=0.8) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~region) +
  labs(x="Expected -log10(p)", y="Observed -log10(p)",
       title="QQ plot: mtCluster term p-values (withAMO; tree-required)") +
  theme_classic(base_size=14)
ggsave(file.path(OUT_DIR, "Fig_QQ_mtCluster_LM.withAMO.png"), p_qq, width=7.2, height=4.2, dpi=300)

if(N_PERM > 0 && file.exists(perm_out_file)){
  PERM <- safe_zread(perm_out_file)
  if("p_perm" %in% names(PERM)){
    p_hist2 <- ggplot(PERM[is.finite(p_perm)], aes(p_perm)) +
      geom_histogram(bins=50) +
      facet_wrap(~region, ncol=1, scales="free_y") +
      labs(x="Permutation p-value for mtCluster term", y="count",
           title=paste0("Permutation test (N_PERM=", N_PERM, "; withAMO; tree-required)")) +
      theme_classic(base_size=14)
    ggsave(file.path(OUT_DIR, "Fig_pHist_mtCluster_perm.withAMO.png"), p_hist2, width=7.2, height=6.5, dpi=300)
  }
}

cat("\n[OK] done. Outputs in:\n  ", OUT_DIR, "\n", sep="")






#perm0
#!/usr/bin/env Rscript
# ============================================================
# 09_mtClusterManual_deltaAF_SNPlevel_perm0.withAMO.v2.R
# deltaAF ~ mtCluster(manual) + treePC1 + treePC2
# per SNP, within region
# NO permutation (perm=0)
#
# Inputs:
#   DELTA_FILE: deltaAF_long.withAMO.tsv.gz
#   CLUSTER_FILE: mtCluster_manual.tsv  (must include AMO row)
#
# Outputs:
#   OUTDIR/
#     LM_perSNP_mtCluster_manual_plus_treePC12.withAMO.tsv.gz
#     LM_cluster_summary_by_region.withAMO.tsv
#     LM_driverCluster_share_sigSNP.withAMO.tsv
#     Fig_pHist_mtCluster_manual.withAMO.png
#     Fig_driverCluster_share_sigSNP.withAMO.png
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

DELTA_FILE <- Sys.getenv("DELTA",
  unset="/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz"
)

# IMPORTANT: default to the withAMO version (has AMO)
CLUSTER_FILE <- Sys.getenv("CLUSTER",
  unset="/mnt/spareHD_2/nu_287/_assoc72_subunit/mtCluster_manual.tsv"
)

OUTDIR <- Sys.getenv("OUTDIR",
  unset="/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_withAMO_perm0"
)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

safe_zread <- function(f){
  if(!file.exists(f)) stop("File not found: ", f)
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
cat("[read] deltaAF_long: ", DELTA_FILE, "\n", sep="")
DEL <- safe_zread(DELTA_FILE)
DEL[, pop := normalize_pop(pop)]

cat("[read] mtCluster_manual: ", CLUSTER_FILE, "\n", sep="")
if(!file.exists(CLUSTER_FILE)) stop("Cluster file not found: ", CLUSTER_FILE)

CL <- fread(CLUSTER_FILE, sep="\t", header=TRUE, fill=TRUE)
CL <- CL[!(is.na(pop) | pop=="")]

need <- c("pop","mtCluster")
miss <- setdiff(need, names(CL))
if(length(miss) > 0){
  stop("Cluster file missing columns: ", paste(miss, collapse=", "))
}

CL[, pop := normalize_pop(pop)]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

# Hard check AMO exists in cluster file
if(!("AMO" %in% CL$pop)){
  cat("\n[FATAL] AMO not found in CLUSTER_FILE.\n")
  cat("CLUSTER_FILE = ", CLUSTER_FILE, "\n", sep="")
  cat("First 30 pops in CLUSTER_FILE:\n")
  print(head(CL$pop, 30))
  stop("Use the withAMO mtCluster_manual.tsv (must contain AMO).")
}

cat("[info] pops per cluster:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# Merge
DEL <- merge(DEL, CL, by="pop", all.x=TRUE)

cat("[info] rows with cluster: ", sum(!is.na(DEL$mtCluster)), " / ", nrow(DEL), "\n", sep="")
if(any(is.na(DEL$mtCluster))){
  cat("[warn] pops missing mtCluster (unique):\n")
  print(unique(DEL[is.na(mtCluster), .(pop)])[order(pop)])
}

# Drop missing cluster rows (should NOT include AMO now)
DEL <- DEL[!is.na(mtCluster)]

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
      mean_deltaAF = mean(deltaAF, na.rm=TRUE),
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

out_lm <- file.path(OUTDIR, "LM_perSNP_mtCluster_manual_plus_treePC12.withAMO.tsv.gz")
fwrite(LM, out_lm, sep="\t")
cat("[write] ", out_lm, "\n", sep="")

# ----------------------------
# Summary
# ----------------------------
SUM <- LM[, .(
  n_snps = .N,
  prop_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
  prop_p10 = mean(p_cluster < 0.10, na.rm=TRUE),
  min_p = suppressWarnings(min(p_cluster, na.rm=TRUE)),
  med_p = suppressWarnings(median(p_cluster, na.rm=TRUE))
), by=region][order(region)]

out_sum <- file.path(OUTDIR, "LM_cluster_summary_by_region.withAMO.tsv")
fwrite(SUM, out_sum, sep="\t")
cat("[write] ", out_sum, "\n", sep="")
cat("\n=== summary ===\n")
print(SUM)

# ----------------------------
# Driver cluster share among p<0.05
# ----------------------------
SIG <- LM[is.finite(p_cluster) & p_cluster < 0.05]
DRIVER <- SIG[, .N, by=.(region, driver_cluster)]
DRIVER[, prop := N / sum(N), by=region]
setorder(DRIVER, region, -prop)

out_drv <- file.path(OUTDIR, "LM_driverCluster_share_sigSNP.withAMO.tsv")
fwrite(DRIVER, out_drv, sep="\t")
cat("[write] ", out_drv, "\n", sep="")

# ----------------------------
# Plots
# ----------------------------
p1 <- ggplot(LM[is.finite(p_cluster)], aes(p_cluster)) +
  geom_histogram(bins=50) +
  facet_wrap(~region, ncol=1) +
  theme_classic(base_size=14) +
  labs(x="p-value (mtCluster F-test)", y="count",
       title="deltaAF ~ mtCluster(manual) + treePC1 + treePC2 (withAMO)")
ggsave(file.path(OUTDIR, "Fig_pHist_mtCluster_manual.withAMO.png"),
       p1, width=7, height=6, dpi=300)

p2 <- ggplot(DRIVER, aes(driver_cluster, prop)) +
  geom_col() +
  facet_wrap(~region, ncol=1, scales="free_x") +
  theme_classic(base_size=14) +
  labs(x="driver cluster (largest |mean deltaAF|)", y="share among significant SNPs",
       title="Which mtCluster drives deltaAF signal? (withAMO)")
ggsave(file.path(OUTDIR, "Fig_driverCluster_share_sigSNP.withAMO.png"),
       p2, width=8, height=5, dpi=300)

cat("\n[OK] done. Outputs in:\n", OUTDIR, "\n", sep="")
DELTA=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz \
CLUSTER=/mnt/spareHD_2/nu_287/_assoc72_subunit/mtCluster_manual.tsv \
OUTDIR=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_withAMO_perm0 \
Rscript /mnt/spareHD_2/nu_287/q2_parallelism/09_mtClusterManual_deltaAF_SNPlevel_perm0.withAMO.v2.R



#bh
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# Gene enrichment from mtCluster-driven SNPs (BH after LM) - withAMO
#
# Inputs:
#   IN_LM: LM_perSNP_mtCluster_manual_plus_treePC12.withAMO.tsv.gz
#         (from 09_mtClusterManual_deltaAF_SNPlevel_perm0.withAMO.*.R)
#
# Outputs (OUTDIR):
#   - LM_with_q_and_driver.withAMO.tsv.gz
#   - Gene_enrichment_overall.withAMO.tsv.gz
#   - GeneCluster_enrichment.withAMO.tsv.gz
#   - Fig_topGenes_byCluster_<region>.withAMO.png
# ============================================================

# ----------------------------
# Config (env override)
# ----------------------------
IN_LM <- Sys.getenv("IN_LM",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_withAMO_perm0/LM_perSNP_mtCluster_manual_plus_treePC12.withAMO.tsv.gz"
)

OUTDIR <- Sys.getenv("OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_withAMO_perm0/gene_enrichment.withAMO"
)

Q_CUT <- as.numeric(Sys.getenv("Q_CUT", unset="0.10"))     # try 0.05 later
MIN_SNP_GENE <- as.integer(Sys.getenv("MIN_SNP_GENE", unset="20"))
REGIONS <- Sys.getenv("REGIONS", unset="")                # e.g. "AK,BC" or "" for all

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
if(!("driver_cluster" %in% names(LM))){
  if("cluster_means" %in% names(LM)){
    infer_driver <- function(x){
      if(is.na(x) || x=="") return(NA_character_)
      parts <- unlist(strsplit(x, ";", fixed=TRUE))
      cl <- sub("^([^:]+):.*$", "\\1", parts)
      mu <- sub("^[^:]+:([^\\(]+)\\(.*$", "\\1", parts)
      mu <- suppressWarnings(as.numeric(mu))
      if(all(!is.finite(mu))) return(NA_character_)
      cl[ which.max(abs(mu)) ]
    }
    LM[, driver_cluster := vapply(cluster_means, infer_driver, character(1))]
    cat("[warn] driver_cluster not found; inferred from cluster_means.\n")
  } else {
    stop("driver_cluster not found and cannot infer (no cluster_means). Re-run 09 to include driver_cluster.")
  }
}

# Save LM with q + driver
out_lm2 <- file.path(OUTDIR, "LM_with_q_and_driver.withAMO.tsv.gz")
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

REG_TOT <- LM[, .(
  N_total=.N,
  N_sig=sum(is_sig)
), by=region]

G <- LM[, .(
  n_gene=.N,
  n_sig_gene=sum(is_sig)
), by=.(region, gene)]

G <- G[n_gene >= MIN_SNP_GENE]

fisher_one <- function(a, b, c, d){
  m <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
  ft <- fisher.test(m)
  list(or=unname(ft$estimate), p=ft$p.value)
}

G[, `:=`(or_overall=NA_real_, p_overall=NA_real_)]

for(r in unique(G$region)){
  tot <- REG_TOT[region==r]
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

out_g_overall <- file.path(OUTDIR, "Gene_enrichment_overall.withAMO.tsv.gz")
fwrite(G[order(region, q_overall, p_overall)], out_g_overall, sep="\t")
cat("[write] ", out_g_overall, "\n", sep="")

# ----------------------------
# 2) Gene × driver_cluster enrichment among significant SNPs
# ----------------------------
cat("\n[run] gene × driver_cluster enrichment (within significant SNPs)\n")

SIG <- LM[is_sig==TRUE & !is.na(driver_cluster)]

if(nrow(SIG)==0){
  cat("[WARN] No significant SNPs at q<", Q_CUT, ". Try Q_CUT=0.10 or 0.20.\n", sep="")
  quit(save="no", status=0)
}

SIG_TOT <- SIG[, .N, by=.(region)]
SIG_CL  <- SIG[, .N, by=.(region, driver_cluster)]

GC <- SIG[, .N, by=.(region, gene, driver_cluster)]
setnames(GC, "N", "n_sig_gene_cl")

GS <- SIG[, .N, by=.(region, gene)]
setnames(GS, "N", "n_sig_gene")
GC <- merge(GC, GS, by=c("region","gene"), all.x=TRUE)

GC <- merge(GC, SIG_TOT, by="region", all.x=TRUE)
setnames(GC, "N", "n_sig_region")
GC <- merge(GC, SIG_CL, by=c("region","driver_cluster"), all.x=TRUE)
setnames(GC, "N", "n_sig_region_cl")

# add universe SNPs per gene for filtering
GC <- merge(GC, G[, .(region, gene, n_gene)], by=c("region","gene"), all.x=TRUE)
GC <- GC[!is.na(n_gene) & n_gene >= MIN_SNP_GENE]

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

out_gc <- file.path(OUTDIR, "GeneCluster_enrichment.withAMO.tsv.gz")
fwrite(GC[order(region, driver_cluster, q_fisher, p_fisher)], out_gc, sep="\t")
cat("[write] ", out_gc, "\n", sep="")

# ----------------------------
# Plot: top enriched genes per cluster per region
# ----------------------------
plot_top <- function(dt, region_name){
  dt <- dt[region == region_name]
  if(nrow(dt)==0) return(NULL)

  dt2 <- dt[is.finite(p_fisher) & is.finite(or_fisher)]
  if(nrow(dt2)==0) return(NULL)

  dt2 <- dt2[order(q_fisher, p_fisher)]
  top <- dt2[, head(.SD, 10), by=driver_cluster]

  top[, gene := factor(gene, levels=rev(unique(gene[order(driver_cluster, q_fisher, p_fisher)])))]

  ggplot(top, aes(x=gene, y=log2(or_fisher))) +
    geom_col() +
    coord_flip() +
    facet_wrap(~driver_cluster, scales="free_y") +
    labs(
      title=paste0("Top gene enrichments by mt driver_cluster (", region_name,
                   ") | sig SNPs: BH(q_cluster)<", Q_CUT, " (withAMO)"),
      x="gene",
      y="log2(Fisher OR) for cluster among significant SNPs"
    ) +
    theme_classic(base_size=13)
}

for(r in unique(GC$region)){
  p <- plot_top(GC, r)
  if(!is.null(p)){
    fn <- file.path(OUTDIR, paste0("Fig_topGenes_byCluster_", r, ".withAMO.png"))
    ggsave(fn, p, width=12, height=7, dpi=300)
    cat("[write] ", fn, "\n", sep="")
  }
}

# ----------------------------
# Console summaries
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


IN_LM=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_withAMO_perm0/LM_perSNP_mtCluster_manual_plus_treePC12.withAMO.tsv.gz \
OUTDIR=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_withAMO_perm0/gene_enrichment.withAMO \
Q_CUT=0.10 MIN_SNP_GENE=20 \
Rscript /mnt/spareHD_2/nu_287/q2_parallelism/10_gene_enrichment_from_cluster_LM.withAMO.R




#delta af～ delta mt
/mnt/spareHD_2/nu_287/q2_parallelism/04_nuDeltaAF_SNPlevel_vs_mtDeltaAF_popSummary.withAMO.noDloop.R
#!/usr/bin/env Rscript
# ============================================================
# 04_nuDeltaAF_SNPlevel_vs_mtDeltaAF_popSummary.withAMO.noDloop.R
#
# Goal (per region AK/BC):
#   For each nuclear SNP, regress nu ΔAF across populations on a
#   population-level summary of mt ΔAF (computed from mtDeltaAF_long).
#
# Models (per SNP, within each region):
#   1) nu ΔAF            ~ mean_mtΔAF       + treePC1 + treePC2
#   2) nu ΔAF            ~ mean_abs_mtΔAF   + treePC1 + treePC2
#   3) nu ΔAF            ~ sd_mtΔAF         + treePC1 + treePC2
#   4) |nu ΔAF| (MAG)    ~ mean_mtΔAF       + treePC1 + treePC2
#   5) |nu ΔAF| (MAG)    ~ mean_abs_mtΔAF   + treePC1 + treePC2
#   6) |nu ΔAF| (MAG)    ~ sd_mtΔAF         + treePC1 + treePC2
#
# Inputs:
#   NU_DELTA: deltaAF_long.withAMO.tsv.gz
#   MT_DELTA: mtDeltaAF_long.noDloop.tsv.gz
#   COV_FILE: covariates.treePC.tsv
#
# Outputs (OUTDIR):
#   popLevel_mtDeltaAF_summary.tsv.gz
#   snpLevel_table.withAMO.noDloop.tsv.gz
#   LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.withAMO.noDloop.tsv.gz
#   LOG_summary.withAMO.noDloop.txt
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ----------------------------
# Inputs (WITH AMO)
# ----------------------------
NU_DELTA <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz"
MT_DELTA <- "/work/cyu/poolseq/PPalign_output/mtDNA_bam/mtDeltaAF_long.noDloop.tsv.gz"
COV_FILE <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

# ----------------------------
# Outputs
# ----------------------------
OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.withAMO.noDloop"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_MTPOP <- file.path(OUTDIR, "popLevel_mtDeltaAF_summary.tsv.gz")
OUT_DT    <- file.path(OUTDIR, "snpLevel_table.withAMO.noDloop.tsv.gz")
OUT_LM    <- file.path(OUTDIR, "LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.withAMO.noDloop.tsv.gz")
OUT_LOG   <- file.path(OUTDIR, "LOG_summary.withAMO.noDloop.txt")

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

# COV fallback
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))
COV <- unique(COV[, .(pop, treePC1, treePC2)])

# Merge: NU (per SNP x pop) + MT_POP (per pop) + COV (treePC)
DT <- merge(NU, MT_POP, by=c("region","pop"), all=FALSE)
DT <- merge(DT, COV, by="pop", all.x=TRUE)
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
DT <- DT[(region=="AK" & n_in_snp >= MIN_POPS_PER_SNP_AK) |
         (region=="BC" & n_in_snp >= MIN_POPS_PER_SNP_BC)]
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
      y  <- dd[[ model_grid$y_col[i] ]]
      x  <- dd[[ model_grid$x_col[i] ]]
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
        response  = model_grid$response_name[i],
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
# 4) Log summary (includes AMO check)
# ----------------------------
amo_rows_dt <- DT[pop=="AMO", .N]
amo_pops_dt <- DT[pop=="AMO", uniqueN(pop)]
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
  sprintf("[CHECK] AMO rows in DT after filters: %d", amo_rows_dt),
  sprintf("[CHECK] AMO present as pop? (uniqueN): %d", amo_pops_dt),
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




#bh
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 08_clean_BH_stable_overlap_pipeline.R
#
# Input: two LM_perSNP tables (withAMO, noAMO)
# Output:
#   - cleaned + BH-adjusted LM tables
#   - stable overlap SNP lists for:
#       (A) magnitude coupling
#       (B) co-directional coupling (pos/neg)
#       (C) mt variability effect
#   - gene summaries (optional hit threshold)
# ============================================================

# -------- paths (edit if needed) --------
LM_WITH <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.withAMO.noDloop/LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.withAMO.noDloop.tsv.gz"
LM_NO   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.noAMO.noDloop/LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.tsv.gz"

OUTDIR  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_cleanBH_stableOverlap.noDloop"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

OUT_WITH_CLEAN <- file.path(OUTDIR, "LM_withAMO.cleanBH.tsv.gz")
OUT_NO_CLEAN   <- file.path(OUTDIR, "LM_noAMO.cleanBH.tsv.gz")

# -------- knobs --------
# choose threshold mode: "q" (recommended) or "p"
THRESH_MODE <- "q"

# thresholds
P_CUT <- 1e-3
Q_CUT <- 0.10

# minimum per-gene hits for gene-level list
MIN_HITS_PER_GENE <- 3

# ============================================================
# helpers
# ============================================================

zread <- function(f){
  if(!file.exists(f)) stop("Missing file: ", f)
  fread(cmd = paste("zcat", shQuote(f)))
}

clean_lm <- function(DT){
  # Ensure numeric
  num_cols <- c("beta","se","t","p","R2","adjR2","n")
  for(cc in intersect(num_cols, names(DT))){
    DT[, (cc) := suppressWarnings(as.numeric(get(cc)))]
  }

  # Drop "un-estimable" fits:
  # - p not finite
  # - beta not finite
  # - se <= 0 or not finite
  # - n < 5 (shouldn't happen but safe)
  DT <- DT[
    is.finite(p) &
    is.finite(beta) &
    is.finite(se) & se > 0 &
    is.finite(n) & n >= 5
  ]

  # Some pipelines accidentally keep blank region/ids
  DT <- DT[nzchar(region) & nzchar(snp) & nzchar(response) & nzchar(predictor)]

  DT
}

add_bh <- function(DT){
  # BH within each region × response × predictor
  DT[, q := p.adjust(p, method="BH"), by=.(region, response, predictor)]
  DT
}

# stable overlap merge
stable_overlap <- function(W, N, response, predictor,
                           direction = c("any","pos","neg"),
                           mode = c("q","p"),
                           p_cut = 1e-3,
                           q_cut = 0.10){

  direction <- match.arg(direction)
  mode <- match.arg(mode)

  Wt <- W[response == response & predictor == predictor]
  Nt <- N[response == response & predictor == predictor]

  if(direction %in% c("pos","neg")){
    if(direction == "pos"){
      Wt <- Wt[beta > 0]
      Nt <- Nt[beta > 0]
    } else {
      Wt <- Wt[beta < 0]
      Nt <- Nt[beta < 0]
    }
  }

  # Apply threshold
  if(mode == "p"){
    Wsig <- Wt[p < p_cut]
    Nsig <- Nt[p < p_cut]
  } else {
    Wsig <- Wt[q < q_cut]
    Nsig <- Nt[q < q_cut]
  }

  # Keep only overlap (stable)
  setkey(Wsig, region, snp)
  setkey(Nsig, region, snp)

  M <- merge(
    Wsig[, .(region, snp, gene, chr, pos,
             beta_with=beta, p_with=p, q_with=q, n_with=n)],
    Nsig[, .(region, snp, gene, chr, pos,
             beta_no=beta, p_no=p, q_no=q, n_no=n)],
    by=c("region","snp"),
    all=FALSE
  )

  # If gene/pos differs across sides, rescue
  M[, gene := fifelse(!is.na(gene.x), gene.x, gene.y)]
  M[, chr  := fifelse(!is.na(chr.x),  chr.x,  chr.y)]
  M[, pos  := fifelse(!is.na(pos.x),  pos.x,  pos.y)]
  M[, c("gene.x","gene.y","chr.x","chr.y","pos.x","pos.y") := NULL]

  setorder(M, region, p_with, p_no)

  M
}

gene_hits <- function(SNP_DT){
  if(nrow(SNP_DT) == 0) return(data.table())
  SNP_DT[, .(
    n_snp = .N,
    min_p_with = min(p_with, na.rm=TRUE),
    min_p_no   = min(p_no,   na.rm=TRUE),
    min_q_with = min(q_with, na.rm=TRUE),
    min_q_no   = min(q_no,   na.rm=TRUE),
    mean_beta_with = mean(beta_with, na.rm=TRUE),
    mean_beta_no   = mean(beta_no,   na.rm=TRUE)
  ), by=.(region, gene)][order(region, -n_snp, min_q_with, min_q_no, min_p_with, min_p_no)]
}

# ============================================================
# main
# ============================================================

cat("[read] withAMO: ", LM_WITH, "\n", sep="")
W0 <- zread(LM_WITH)

cat("[read] noAMO  : ", LM_NO, "\n", sep="")
N0 <- zread(LM_NO)

cat("[clean] dropping non-estimable fits\n")
W <- add_bh(clean_lm(copy(W0)))
N <- add_bh(clean_lm(copy(N0)))

fwrite(W, OUT_WITH_CLEAN, sep="\t", compress="gzip")
fwrite(N, OUT_NO_CLEAN,   sep="\t", compress="gzip")
cat("[write] ", OUT_WITH_CLEAN, "\n", sep="")
cat("[write] ", OUT_NO_CLEAN, "\n", sep="")

cat("\n[info] rows kept:\n")
cat("  withAMO:", nrow(W), "\n")
cat("  noAMO  :", nrow(N), "\n")

# -------- define three target tests --------
targets <- list(
  magnitude = list(response="abs_nuDeltaAF", predictor="mean_abs_mtDeltaAF", direction="any"),
  codir_pos = list(response="nuDeltaAF",     predictor="mean_mtDeltaAF",     direction="pos"),
  codir_neg = list(response="nuDeltaAF",     predictor="mean_mtDeltaAF",     direction="neg"),
  mtvar     = list(response="abs_nuDeltaAF", predictor="sd_mtDeltaAF",       direction="any")
)

# Run all targets
for(nm in names(targets)){
  tg <- targets[[nm]]

  cat("\n[run] target: ", nm,
      " | response=", tg$response,
      " predictor=", tg$predictor,
      " direction=", tg$direction,
      " | mode=", THRESH_MODE,
      "\n", sep="")

  ST <- stable_overlap(
    W, N,
    response = tg$response,
    predictor = tg$predictor,
    direction = tg$direction,
    mode = THRESH_MODE,
    p_cut = P_CUT,
    q_cut = Q_CUT
  )

  out_snp  <- file.path(OUTDIR, sprintf("stable_SNPs_%s.mode_%s.tsv", nm, THRESH_MODE))
  out_gene <- file.path(OUTDIR, sprintf("stable_genes_%s.mode_%s.minHits%d.tsv", nm, THRESH_MODE, MIN_HITS_PER_GENE))

  fwrite(ST, out_snp, sep="\t")
  cat("[write] ", out_snp, " (n=", nrow(ST), ")\n", sep="")

  G <- gene_hits(ST)
  Gk <- G[n_snp >= MIN_HITS_PER_GENE]
  fwrite(Gk, out_gene, sep="\t")
  cat("[write] ", out_gene, " (genes=", nrow(Gk), ")\n", sep="")

  # quick console
  if(nrow(ST) > 0){
    cat("[top SNPs]\n")
    print(head(ST[, .(region, snp, gene, beta_with, p_with, q_with, beta_no, p_no, q_no)], 10))
  } else {
    cat("[none]\n")
  }
}

cat("\n[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")
cat("Threshold mode: ", THRESH_MODE, " | P_CUT=", P_CUT, " | Q_CUT=", Q_CUT, "\n", sep="")




#cluster
suppressPackageStartupMessages(library(data.table))

# WITH AMO
DELTA_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz"
CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mtCluster_manual.tsv"

normalize_pop <- function(x){
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
}

# ---------- read DEL ----------
DEL <- fread(cmd=paste("zcat", shQuote(DELTA_FILE)))
DEL[, pop := normalize_pop(pop)]

# ---------- read CL (robust) ----------
CL <- fread(CLUSTER_FILE, sep="\t", header=FALSE, fill=TRUE, strip.white=TRUE)
# keep first 2 columns (pop, cluster), drop empty lines
CL <- CL[nzchar(V1) & nzchar(V2), .(pop=V1, cluster=V2)]
# drop header line if present (common when header was read as data)
CL <- CL[pop != "pop" & cluster != "mtCluster"]
CL[, pop := normalize_pop(pop)]
CL[, cluster := factor(cluster)]
CL <- unique(CL, by="pop")

cat("[info] clusters loaded:\n")
print(CL[, .N, by=cluster][order(-N)])

# ---------- merge ----------
DEL <- merge(DEL, CL, by="pop", all.x=TRUE)
DEL <- DEL[!is.na(cluster)]

# ---------- pop x snp matrix per region ----------
W <- dcast(DEL, pop + region + cluster ~ snp, value.var="deltaAF")

# ---------- permutation test ----------
perm_test <- function(dist_obj, cl, nperm=2000){
  M <- as.matrix(dist_obj)
  idx <- which(upper.tri(M), arr.ind=TRUE)
  d <- M[upper.tri(M)]
  same <- cl[idx[,1]] == cl[idx[,2]]

  if(sum(same, na.rm=TRUE) == 0 || sum(!same, na.rm=TRUE) == 0){
    return(list(
      obs=NA_real_, p=NA_real_, within=NA_real_, between=NA_real_,
      n_perm_used=0L, note="Not enough within/between pairs (cluster sizes too small?)"
    ))
  }
  within_obs  <- mean(d[same],  na.rm=TRUE)
  between_obs <- mean(d[!same], na.rm=TRUE)
  obs <- within_obs - between_obs  # <0 => within more similar

  ge <- 0L
  used <- 0L
  for(b in seq_len(nperm)){
    clp <- sample(cl)
    samep <- clp[idx[,1]] == clp[idx[,2]]
    if(sum(samep, na.rm=TRUE) == 0 || sum(!samep, na.rm=TRUE) == 0) next

    within_p  <- mean(d[samep],  na.rm=TRUE)
    between_p <- mean(d[!samep], na.rm=TRUE)
    statp <- within_p - between_p
    if(!is.finite(statp)) next

    used <- used + 1L
    if(statp <= obs) ge <- ge + 1L
  }

  p <- if(used>0) (1 + ge)/(1 + used) else NA_real_
  list(obs=obs, p=p, within=within_obs, between=between_obs, n_perm_used=used)
}

# ---------- run per region ----------
for(r in unique(W$region)){
  X <- W[region==r]
  mat <- as.matrix(X[, -(1:3)])
  rownames(mat) <- X$pop
  cl <- X$cluster

  cor_signed <- cor(t(mat), use="pairwise.complete.obs")
  out1 <- perm_test(as.dist(1 - cor_signed), cl, nperm=2000)
  cat("\n===", r, "signed ΔAF (1-cor) ===\n"); print(out1)

  cor_abs <- cor(t(abs(mat)), use="pairwise.complete.obs")
  out2 <- perm_test(as.dist(1 - cor_abs), cl, nperm=2000)
  cat("\n===", r, "abs(ΔAF) (1-cor) ===\n"); print(out2)
}