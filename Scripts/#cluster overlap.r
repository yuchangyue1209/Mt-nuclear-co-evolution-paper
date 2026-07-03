#cluster overlap
set -euo pipefail

P=0.05

WI_DIR=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_withAMO_perm0
NO_DIR=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0

# 自动找 LM 文件（匹配 treePC12 + mtCluster manual）
WI_CL=$(ls $WI_DIR/LM_perSNP*mtCluster*treePC12*.tsv.gz 2>/dev/null | head -n 1 || true)
NO_CL=$(ls $NO_DIR/LM_perSNP*mtCluster*treePC12*.tsv.gz 2>/dev/null | head -n 1 || true)

echo "[withAMO] $WI_CL"
echo "[noAMO]   $NO_CL"

# 如果任意一个没找到，就直接报错并列出目录里所有 LM 方便你定位
if [[ -z "${WI_CL}" || -z "${NO_CL}" ]]; then
  echo
  echo "[ERROR] LM file not found by glob. Listing LM-like files:"
  ls -lh "$WI_DIR" | grep -E 'LM_perSNP|LM_' || true
  ls -lh "$NO_DIR" | grep -E 'LM_perSNP|LM_' || true
  exit 1
fi

# overlap（同一 region+snp，且两边都 p_cluster < P）
join -t $'\t' -1 1 -2 1 \
  <(zcat "$WI_CL" | awk -v P="$P" 'BEGIN{FS=OFS="\t"}
     NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
     {r=$h["region"]; snp=$h["snp"]; p=$h["p_cluster"]; if(p==""||p=="NA") next; p+=0;
      if(p<P){
        key=r SUBSEP snp;
        print key,r,$h["gene"],$h["chr"],$h["pos"],snp,p,$h["driver_cluster"];
      }}' | sort -t$'\t' -k1,1) \
  <(zcat "$NO_CL" | awk -v P="$P" 'BEGIN{FS=OFS="\t"}
     NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
     {r=$h["region"]; snp=$h["snp"]; p=$h["p_cluster"]; if(p==""||p=="NA") next; p+=0;
      if(p<P){
        key=r SUBSEP snp;
        print key,r,$h["gene"],$h["chr"],$h["pos"],snp,p,$h["driver_cluster"];
      }}' | sort -t$'\t' -k1,1) \
| awk 'BEGIN{FS=OFS="\t"}
  { # key + with(2..8) + no(9..15)
    # 输出：region gene chr pos snp p_noAMO p_withAMO driver_noAMO driver_withAMO
    print $2,$3,$4,$5,$6,$14,$7,$15,$8,$16
  }' \
| (echo -e "region\tgene\tchr\tpos\tsnp\tp_noAMO\tp_withAMO\tdriver_noAMO\tdriver_withAMO"; cat) \
| column -t -s$'\t' | head -n 25


/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_withAMO_perm0/overlap_snps_p0.05_with_vs_noAMO.tsv.gz







/mnt/spareHD_2/nu_287/q2_parallelism/05_compare_candidates_withAMO_vs_noAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

LM_WITH <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.withAMO.noDloop/LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.withAMO.noDloop.tsv.gz"
LM_NO   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.noAMO.noDloop/LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.tsv.gz"

OUTDIR  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_candidate_compare_with_vs_noAMO.noDloop"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# --- choose the ONE model you should prioritize for candidates ---
TARGET_RESPONSE  <- "abs_nuDeltaAF"
TARGET_PREDICTOR <- "mean_abs_mtDeltaAF"

# --- thresholds (tune later if too many/too few) ---
P_SNP   <- 1e-3   # SNP-level candidate p cutoff
MIN_HITS_PER_GENE <- 3  # gene must have >= this many SNP candidates

read_lm <- function(path){
  if(!file.exists(path)) stop("Missing: ", path)
  fread(cmd=paste("zcat", shQuote(path)))
}

prep <- function(DT){
  DT <- DT[response == TARGET_RESPONSE & predictor == TARGET_PREDICTOR]
  DT <- DT[is.finite(p) & is.finite(beta)]
  DT
}

cat("[read] withAMO:", LM_WITH, "\n")
W <- prep(read_lm(LM_WITH))
cat("[read] noAMO  :", LM_NO, "\n")
N <- prep(read_lm(LM_NO))

# =========================
# SNP-level candidates
# =========================
W_snp <- W[p < P_SNP, .(region, gene, chr, pos, snp, beta, p, n)]
N_snp <- N[p < P_SNP, .(region, gene, chr, pos, snp, beta, p, n)]

setkey(W_snp, region, snp)
setkey(N_snp, region, snp)

SNP_M <- merge(
  W_snp, N_snp,
  by=c("region","snp"),
  all=TRUE,
  suffixes=c("_with","_no")
)

SNP_M[, status :=
  fifelse(!is.na(p_with) & !is.na(p_no), "stable",
  fifelse(!is.na(p_with) &  is.na(p_no), "with_only",
  fifelse( is.na(p_with) & !is.na(p_no), "no_only", "none")))
]

# clean gene/pos if one side missing
SNP_M[, gene := fifelse(!is.na(gene_with), gene_with, gene_no)]
SNP_M[, chr  := fifelse(!is.na(chr_with),  chr_with,  chr_no)]
SNP_M[, pos  := fifelse(!is.na(pos_with),  pos_with,  pos_no)]

SNP_M[, `:=`(gene_with=NULL, gene_no=NULL, chr_with=NULL, chr_no=NULL, pos_with=NULL, pos_no=NULL)]

fwrite(SNP_M[order(region, status, p_with, p_no)], file.path(OUTDIR, "candidate_SNPs_with_vs_noAMO.tsv"), sep="\t")

# =========================
# Gene-level candidates
# =========================
gene_sum <- function(SNP_M, which_status){
  dd <- SNP_M[status %in% which_status]
  dd[, .(
    n_snp = .N,
    min_p_with = suppressWarnings(min(p_with, na.rm=TRUE)),
    min_p_no   = suppressWarnings(min(p_no,   na.rm=TRUE)),
    mean_beta_with = mean(beta_with, na.rm=TRUE),
    mean_beta_no   = mean(beta_no,   na.rm=TRUE)
  ), by=.(region, gene)][order(region, -n_snp, min_p_with, min_p_no)]
}

G_stable   <- gene_sum(SNP_M, "stable")
G_withonly <- gene_sum(SNP_M, "with_only")
G_noonly   <- gene_sum(SNP_M, "no_only")

# apply hit-count threshold for genes
G_stable_keep   <- G_stable[n_snp >= MIN_HITS_PER_GENE]
G_withonly_keep <- G_withonly[n_snp >= MIN_HITS_PER_GENE]
G_noonly_keep   <- G_noonly[n_snp >= MIN_HITS_PER_GENE]

fwrite(G_stable_keep,   file.path(OUTDIR, "candidate_genes_STABLE.tsv"), sep="\t")
fwrite(G_withonly_keep, file.path(OUTDIR, "candidate_genes_WITH_ONLY.tsv"), sep="\t")
fwrite(G_noonly_keep,   file.path(OUTDIR, "candidate_genes_NO_ONLY.tsv"), sep="\t")

# quick console summary
cat("\n=== MODEL USED ===\n")
cat("response :", TARGET_RESPONSE, "\n")
cat("predictor:", TARGET_PREDICTOR, "\n")
cat("P_SNP    :", P_SNP, "\n")
cat("MIN_HITS :", MIN_HITS_PER_GENE, "\n")

cat("\n=== SNP candidates counts ===\n")
print(SNP_M[, .N, by=.(region, status)][order(region, status)])

cat("\n=== Gene candidates (>= MIN_HITS) counts ===\n")
cat("stable   :", nrow(G_stable_keep), "\n")
cat("with_only:", nrow(G_withonly_keep), "\n")
cat("no_only  :", nrow(G_noonly_keep), "\n")

cat("\n[OK] wrote to:\n  ", OUTDIR, "\n", sep="")









#方向一致（co-directional）脚本 05_compare_candidates_codirectional_withAMO_vs_noAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

LM_WITH <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.withAMO.noDloop/LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.withAMO.noDloop.tsv.gz"
LM_NO   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.noAMO.noDloop/LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.tsv.gz"

OUTDIR  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_candidate_compare_codirectional_with_vs_noAMO.noDloop"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# ------------------------------------------------------------
# TARGET MODEL = co-directional coupling
#   nuDeltaAF ~ mean_mtDeltaAF + treePC1 + treePC2
# ------------------------------------------------------------
TARGET_RESPONSE  <- "nuDeltaAF"
TARGET_PREDICTOR <- "mean_mtDeltaAF"

# thresholds
P_SNP <- 1e-3
MIN_HITS_PER_GENE <- 3

# direction filter options:
#   "POS"   -> beta > 0 (co-directional positive coupling)
#   "NEG"   -> beta < 0 (anti-directional / opposite coupling)
#   "BOTH"  -> keep both signs
DIRECTION_MODE <- "POS"

read_lm <- function(path){
  if(!file.exists(path)) stop("Missing: ", path)
  fread(cmd=paste("zcat", shQuote(path)))
}

prep <- function(DT){
  DT <- DT[response == TARGET_RESPONSE & predictor == TARGET_PREDICTOR]
  DT <- DT[is.finite(p) & is.finite(beta)]
  if(DIRECTION_MODE == "POS"){
    DT <- DT[beta > 0]
  } else if(DIRECTION_MODE == "NEG"){
    DT <- DT[beta < 0]
  } else if(DIRECTION_MODE == "BOTH"){
    # keep all
  } else {
    stop("Unknown DIRECTION_MODE: ", DIRECTION_MODE, " (use POS/NEG/BOTH)")
  }
  DT
}

cat("[read] withAMO:", LM_WITH, "\n")
W <- prep(read_lm(LM_WITH))
cat("[read] noAMO  :", LM_NO, "\n")
N <- prep(read_lm(LM_NO))

# SNP-level candidates
W_snp <- W[p < P_SNP, .(region, gene, chr, pos, snp, beta, p, n)]
N_snp <- N[p < P_SNP, .(region, gene, chr, pos, snp, beta, p, n)]

setkey(W_snp, region, snp)
setkey(N_snp, region, snp)

SNP_M <- merge(
  W_snp, N_snp,
  by=c("region","snp"),
  all=TRUE,
  suffixes=c("_with","_no")
)

SNP_M[, status :=
  fifelse(!is.na(p_with) & !is.na(p_no), "stable",
  fifelse(!is.na(p_with) &  is.na(p_no), "with_only",
  fifelse( is.na(p_with) & !is.na(p_no), "no_only", "none")))
]

# clean gene/pos if one side missing
SNP_M[, gene := fifelse(!is.na(gene_with), gene_with, gene_no)]
SNP_M[, chr  := fifelse(!is.na(chr_with),  chr_with,  chr_no)]
SNP_M[, pos  := fifelse(!is.na(pos_with),  pos_with,  pos_no)]

SNP_M[, `:=`(gene_with=NULL, gene_no=NULL, chr_with=NULL, chr_no=NULL, pos_with=NULL, pos_no=NULL)]

fwrite(SNP_M[order(region, status, p_with, p_no)],
       file.path(OUTDIR, "candidate_SNPs_with_vs_noAMO.tsv"),
       sep="\t")

# Gene-level candidates
gene_sum <- function(SNP_M, which_status){
  dd <- SNP_M[status %in% which_status]
  dd[, .(
    n_snp = .N,
    min_p_with = suppressWarnings(min(p_with, na.rm=TRUE)),
    min_p_no   = suppressWarnings(min(p_no,   na.rm=TRUE)),
    mean_beta_with = mean(beta_with, na.rm=TRUE),
    mean_beta_no   = mean(beta_no,   na.rm=TRUE)
  ), by=.(region, gene)][order(region, -n_snp, min_p_with, min_p_no)]
}

G_stable   <- gene_sum(SNP_M, "stable")
G_withonly <- gene_sum(SNP_M, "with_only")
G_noonly   <- gene_sum(SNP_M, "no_only")

G_stable_keep   <- G_stable[n_snp >= MIN_HITS_PER_GENE]
G_withonly_keep <- G_withonly[n_snp >= MIN_HITS_PER_GENE]
G_noonly_keep   <- G_noonly[n_snp >= MIN_HITS_PER_GENE]

fwrite(G_stable_keep,   file.path(OUTDIR, "candidate_genes_STABLE.tsv"), sep="\t")
fwrite(G_withonly_keep, file.path(OUTDIR, "candidate_genes_WITH_ONLY.tsv"), sep="\t")
fwrite(G_noonly_keep,   file.path(OUTDIR, "candidate_genes_NO_ONLY.tsv"), sep="\t")

cat("\n=== MODEL USED ===\n")
cat("response :", TARGET_RESPONSE, "\n")
cat("predictor:", TARGET_PREDICTOR, "\n")
cat("direction:", DIRECTION_MODE, "\n")
cat("P_SNP    :", P_SNP, "\n")
cat("MIN_HITS :", MIN_HITS_PER_GENE, "\n")

cat("\n=== SNP candidates counts ===\n")
print(SNP_M[, .N, by=.(region, status)][order(region, status)])

cat("\n=== Gene candidates (>= MIN_HITS) counts ===\n")
cat("stable   :", nrow(G_stable_keep), "\n")
cat("with_only:", nrow(G_withonly_keep), "\n")
cat("no_only  :", nrow(G_noonly_keep), "\n")

cat("\n[OK] wrote to:\n  ", OUTDIR, "\n", sep="")

#mt variability effect 脚本（强度版，推荐）05_compare_candidates_mtVariability_withAMO_vs_noAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

LM_WITH <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.withAMO.noDloop/LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.withAMO.noDloop.tsv.gz"
LM_NO   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_nuDeltaAF_vs_mtDeltaAF_SNPlevel.noAMO.noDloop/LM_perSNP_nuDeltaAF_vs_mtSummary_plus_treePC12.tsv.gz"

OUTDIR  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_candidate_compare_mtVariability_with_vs_noAMO.noDloop"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# ------------------------------------------------------------
# TARGET MODEL = mt variability effect (magnitude version)
#   abs_nuDeltaAF ~ sd_mtDeltaAF + treePC1 + treePC2
# ------------------------------------------------------------
TARGET_RESPONSE  <- "abs_nuDeltaAF"
TARGET_PREDICTOR <- "sd_mtDeltaAF"

# thresholds
P_SNP <- 1e-3
MIN_HITS_PER_GENE <- 3

read_lm <- function(path){
  if(!file.exists(path)) stop("Missing: ", path)
  fread(cmd=paste("zcat", shQuote(path)))
}

prep <- function(DT){
  DT <- DT[response == TARGET_RESPONSE & predictor == TARGET_PREDICTOR]
  DT <- DT[is.finite(p) & is.finite(beta)]
  DT
}

cat("[read] withAMO:", LM_WITH, "\n")
W <- prep(read_lm(LM_WITH))
cat("[read] noAMO  :", LM_NO, "\n")
N <- prep(read_lm(LM_NO))

# SNP-level candidates
W_snp <- W[p < P_SNP, .(region, gene, chr, pos, snp, beta, p, n)]
N_snp <- N[p < P_SNP, .(region, gene, chr, pos, snp, beta, p, n)]

setkey(W_snp, region, snp)
setkey(N_snp, region, snp)

SNP_M <- merge(
  W_snp, N_snp,
  by=c("region","snp"),
  all=TRUE,
  suffixes=c("_with","_no")
)

SNP_M[, status :=
  fifelse(!is.na(p_with) & !is.na(p_no), "stable",
  fifelse(!is.na(p_with) &  is.na(p_no), "with_only",
  fifelse( is.na(p_with) & !is.na(p_no), "no_only", "none")))
]

# clean gene/pos if one side missing
SNP_M[, gene := fifelse(!is.na(gene_with), gene_with, gene_no)]
SNP_M[, chr  := fifelse(!is.na(chr_with),  chr_with,  chr_no)]
SNP_M[, pos  := fifelse(!is.na(pos_with),  pos_with,  pos_no)]

SNP_M[, `:=`(gene_with=NULL, gene_no=NULL, chr_with=NULL, chr_no=NULL, pos_with=NULL, pos_no=NULL)]

fwrite(SNP_M[order(region, status, p_with, p_no)],
       file.path(OUTDIR, "candidate_SNPs_with_vs_noAMO.tsv"),
       sep="\t")

# Gene-level candidates
gene_sum <- function(SNP_M, which_status){
  dd <- SNP_M[status %in% which_status]
  dd[, .(
    n_snp = .N,
    min_p_with = suppressWarnings(min(p_with, na.rm=TRUE)),
    min_p_no   = suppressWarnings(min(p_no,   na.rm=TRUE)),
    mean_beta_with = mean(beta_with, na.rm=TRUE),
    mean_beta_no   = mean(beta_no,   na.rm=TRUE)
  ), by=.(region, gene)][order(region, -n_snp, min_p_with, min_p_no)]
}

G_stable   <- gene_sum(SNP_M, "stable")
G_withonly <- gene_sum(SNP_M, "with_only")
G_noonly   <- gene_sum(SNP_M, "no_only")

G_stable_keep   <- G_stable[n_snp >= MIN_HITS_PER_GENE]
G_withonly_keep <- G_withonly[n_snp >= MIN_HITS_PER_GENE]
G_noonly_keep   <- G_noonly[n_snp >= MIN_HITS_PER_GENE]

fwrite(G_stable_keep,   file.path(OUTDIR, "candidate_genes_STABLE.tsv"), sep="\t")
fwrite(G_withonly_keep, file.path(OUTDIR, "candidate_genes_WITH_ONLY.tsv"), sep="\t")
fwrite(G_noonly_keep,   file.path(OUTDIR, "candidate_genes_NO_ONLY.tsv"), sep="\t")

cat("\n=== MODEL USED ===\n")
cat("response :", TARGET_RESPONSE, "\n")
cat("predictor:", TARGET_PREDICTOR, "\n")
cat("P_SNP    :", P_SNP, "\n")
cat("MIN_HITS :", MIN_HITS_PER_GENE, "\n")

cat("\n=== SNP candidates counts ===\n")
print(SNP_M[, .N, by=.(region, status)][order(region, status)])

cat("\n=== Gene candidates (>= MIN_HITS) counts ===\n")
cat("stable   :", nrow(G_stable_keep), "\n")
cat("with_only:", nrow(G_withonly_keep), "\n")
cat("no_only  :", nrow(G_noonly_keep), "\n")

cat("\n[OK] wrote to:\n  ", OUTDIR, "\n", sep="")


A. magnitude coupling（你之前的主结果）

abs_nuDeltaAF ~ mean_abs_mtDeltaAF
→ mt divergence strength 与 nu divergence strength 的耦合

你之前 stable p<0.001 是 11 个 SNP，基因里 sdha 这种是核心例子。

B. co-directional coupling（刚跑的）

nuDeltaAF ~ mean_mtDeltaAF 且 beta>0
→ 方向一致的耦合存在，但较分散（AK stable=9；gene-level=0）

C. mt variability effect（刚跑的）

abs_nuDeltaAF ~ sd_mtDeltaAF
→ mt 变异度越大，nu 响应强度越大（AK strong；gene-level=1）


(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ cut -f2 q2_candidate_compare_with_vs_noAMO.noDloop/candidate_SNPs_with_vs_noAMO.tsv | tail -n +2 | sort > mag.txt

cut -f2 q2_candidate_compare_codirectional_with_vs_noAMO.noDloop/candidate_SNPs_with_vs_noAMO.tsv | tail -n +2 | sort > dir.txt

cut -f2 q2_candidate_compare_mtVariability_with_vs_noAMO.noDloop/candidate_SNPs_with_vs_noAMO.tsv | tail -n +2 | sort > var.txt

comm -12 mag.txt dir.txt | comm -12 - var.txt
chrII:22432377:ndufs3
chrIII:13902052:ndufa12
chrIX:7621946:ndufb7
chrV:14311213:uqcrc2b
chrVIII:10732320:ndufs7
chrXII:17117141:cox4i2
chrXIII:14693050:cox6a1
chrXVI:9550742:ndufs1
chrXX:14455933:sdha
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ 



abab
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ cut -f2 \
/mnt/spareHD_2/nu_287/q2_parallelism/q2_candidate_compare_with_vs_noAMO.noDloop/candidate_SNPs_with_vs_noAMO.tsv \
| sort | uniq -c
      1 chrII:22432377:ndufs3
      1 chrIII:13902052:ndufa12
      1 chrIX:7621946:ndufb7
      1 chrV:14311213:uqcrc2b
      1 chrVIII:10732320:ndufs7
      1 chrXII:17117141:cox4i2
      1 chrXIII:14693050:cox6a1
      1 chrXVI:9550742:ndufs1
      1 chrXX:14453650:sdha
      1 chrXX:14455917:sdha
      1 chrXX:14455933:sdha

delta Delta
cut -f2 \
q2_candidate_compare_codirectional_with_vs_noAMO.noDloop/candidate_SNPs_with_vs_noAMO.tsv \
| tail -n +2 \
| sort | uniq -c
      1 chrII:22432377:ndufs3
      1 chrIII:13902052:ndufa12
      1 chrIX:7621946:ndufb7
      1 chrV:14311213:uqcrc2b
      1 chrVIII:10732320:ndufs7
      1 chrXII:17117141:cox4i2
      1 chrXIII:14693050:cox6a1
      1 chrXVI:9550742:ndufs1
      1 chrXX:14455933:sdha


ab   sddelta
cut -f2 \
q2_candidate_compare_mtVariability_with_vs_noAMO.noDloop/candidate_SNPs_with_vs_noAMO.tsv \
| tail -n +2 \
| sort | uniq -c
      1 chrII:22432377:ndufs3
      1 chrIII:13902052:ndufa12
      1 chrIX:7621946:ndufb7
      1 chrV:14311213:uqcrc2b
      1 chrVIII:10732320:ndufs7
      1 chrXII:17117141:cox4i2
      1 chrXIII:14693050:cox6a1
      1 chrXVI:9550742:ndufs1
      1 chrXVII:190125:cox6c
      1 chrXVIII:7633850:ndufb1
      1 chrXX:14453650:sdha
      1 chrXX:14455917:sdha
      1 chrXX:14455933:sdha



overlap 9snps
q2_candidate_compare_codirectional_with_vs_noAMO.noDloop/candidate_SNPs_with_vs_noAMO.tsv \
| tail -n +2 \
| sort | uniq -c
      1 chrII:22432377:ndufs3
      1 chrIII:13902052:ndufa12
      1 chrIX:7621946:ndufb7
      1 chrV:14311213:uqcrc2b
      1 chrVIII:10732320:ndufs7
      1 chrXII:17117141:cox4i2
      1 chrXIII:14693050:cox6a1
      1 chrXVI:9550742:ndufs1
      1 chrXX:14455933:sdha