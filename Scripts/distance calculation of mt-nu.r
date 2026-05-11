距离
suppressPackageStartupMessages(library(data.table))

DELTA_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"

normalize_pop <- function(x){
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
}

# ---------- read DEL ----------
DEL <- fread(cmd=paste("zcat", shQuote(DELTA_FILE)))
DEL[, pop := normalize_pop(pop)]

# ---------- read CL (robust) ----------
CL <- fread(CLUSTER_FILE, sep="\t", header=FALSE, fill=TRUE, strip.white=TRUE)
# 只保留前两列（pop, cluster），去掉空行
CL <- CL[nzchar(V1) & nzchar(V2), .(pop=V1, cluster=V2)]
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

  # 观测统计量：within - between（越小越“同cluster更相似”）
  if(sum(same, na.rm=TRUE) == 0 || sum(!same, na.rm=TRUE) == 0){
    return(list(
      obs=NA_real_, p=NA_real_, within=NA_real_, between=NA_real_,
      n_perm_used=0L, note="Not enough within/between pairs (cluster sizes too small?)"
    ))
  }
  within_obs  <- mean(d[same],  na.rm=TRUE)
  between_obs <- mean(d[!same], na.rm=TRUE)
  obs <- within_obs - between_obs

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

  # signed ΔAF similarity: 1 - correlation
  cor_signed <- cor(t(mat), use="pairwise.complete.obs")
  out1 <- perm_test(as.dist(1 - cor_signed), cl, nperm=2000)
  cat("\n===", r, "signed ΔAF (1-cor) ===\n"); print(out1)

  # abs(ΔAF) similarity: 1 - correlation
  cor_abs <- cor(t(abs(mat)), use="pairwise.complete.obs")
  out2 <- perm_test(as.dist(1 - cor_abs), cl, nperm=2000)
  cat("\n===", r, "abs(ΔAF) (1-cor) ===\n"); print(out2)
}

library(ggplot2)
library(data.table)

plot_pairwise <- function(mat, cl, region_name){

  cor_mat <- cor(t(mat), use="pairwise.complete.obs")
  dist_mat <- 1 - cor_mat

  M <- as.matrix(dist_mat)
  idx <- which(upper.tri(M), arr.ind=TRUE)

  dt <- data.table(
    pop1 = rownames(M)[idx[,1]],
    pop2 = colnames(M)[idx[,2]],
    dist = M[upper.tri(M)]
  )

  dt[, type := ifelse(cl[idx[,1]] == cl[idx[,2]], "Within", "Between")]
  dt[, region := region_name]

  return(dt)
}

ALL <- list()

for(r in unique(W$region)){
  X <- W[region==r]
  mat <- as.matrix(X[, -(1:3)])
  rownames(mat) <- X$pop
  cl <- X$cluster

  ALL[[r]] <- plot_pairwise(mat, cl, r)
}

PLOT_DT <- rbindlist(ALL)

p <- ggplot(PLOT_DT, aes(x=type, y=dist, fill=type)) +
  geom_boxplot(outlier.shape=NA, width=0.6) +
  geom_jitter(width=0.1, alpha=0.3, size=1) +
  facet_wrap(~region, scales="free_y") +
  theme_classic(base_size=14) +
  labs(
    x=NULL,
    y="1 − correlation (ΔAF similarity)",
    title="Mitochondrial lineage structures nuclear evolutionary trajectories"
  )

ggsave("Fig_pairwise_similarity_boxplot.png", p, width=6, height=4, dpi=300)

#plot
library(data.table)
library(ggplot2)

PAIRWISE <- fread("/mnt/spareHD_2/nu_287/q2_parallelism/q2_pairwise_deltaAF_similarity_ldpruned_r2_0.2/pairwise_deltaAF_similarity_within_between.tsv")
SUMMARY  <- fread("/mnt/spareHD_2/nu_287/q2_parallelism/q2_pairwise_deltaAF_similarity_ldpruned_r2_0.2/pairwise_deltaAF_similarity_summary.tsv")

# only signed ΔAF
P <- PAIRWISE[metric == "Signed_deltaAF"]
P[, type := factor(type, levels=c("Within", "Between"))]

# mean ± SE for bars
BAR <- P[, .(
  mean_dist = mean(dist, na.rm=TRUE),
  se_dist = sd(dist, na.rm=TRUE) / sqrt(.N),
  n_pairs = .N
), by=.(region, type)]

# p labels
LAB <- SUMMARY[metric == "Signed_deltaAF"]
LAB[, label := ifelse(p_perm < 0.05,
                      paste0("p = ", signif(p_perm, 2), " *"),
                      paste0("p = ", signif(p_perm, 2)))]
ypos <- BAR[, .(y=max(mean_dist + se_dist, na.rm=TRUE) * 1.08), by=region]
LAB <- merge(LAB, ypos, by="region")

p <- ggplot(BAR, aes(x=type, y=mean_dist, fill=type)) +
  geom_col(width=0.6, alpha=0.85) +
  geom_errorbar(
    aes(ymin=mean_dist - se_dist, ymax=mean_dist + se_dist),
    width=0.18,
    linewidth=0.7
  ) +
  geom_text(
    data=LAB,
    aes(x=1.5, y=y, label=label),
    inherit.aes=FALSE,
    size=4
  ) +
  facet_wrap(~region, scales="free_y") +
  theme_classic(base_size=14) +
  labs(
    x=NULL,
    y="1 - correlation of ΔAF profiles",
    
  ) +
  theme(
    legend.position="none",
    strip.text=element_text(face="bold"),
    plot.title=element_text(face="bold")
  )

OUTFIG <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_pairwise_deltaAF_similarity_ldpruned_r2_0.2/Fig_pairwise_similarity_barplot_signed_only_SE.png"
ggsave(OUTFIG, p, width=6, height=4, dpi=300)

cat("saved:", OUTFIG, "\n")
print(p)

library(ggplot2)

DELTA_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_withAMO/deltaAF_long.withAMO.tsv.gz"
OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/figs_deltaAF_pop_extreme"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

DEL <- read.table(gzfile(DELTA_FILE), header=TRUE, sep="\t")

POPSTAT <- aggregate(deltaAF ~ region + pop,
                     data=DEL,
                     FUN=function(x){
                       c(mean=mean(x, na.rm=TRUE),
                         mean_abs=mean(abs(x), na.rm=TRUE),
                         sd=sd(x, na.rm=TRUE),
                         n=length(x))
                     })

# 展开 list 列
POPSTAT <- do.call(data.frame, POPSTAT)
colnames(POPSTAT) <- c("region","pop",
                       "mean_deltaAF",
                       "mean_abs_deltaAF",
                       "sd_deltaAF",
                       "n_snps")

POPSTAT <- POPSTAT[order(-POPSTAT$mean_abs_deltaAF), ]

write.table(POPSTAT,
            file=file.path(OUTDIR,"POP_deltaAF_summary.tsv"),
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

head(POPSTAT)



head(POPSTAT)
   region  pop  mean_deltaAF mean_abs_deltaAF sd_deltaAF n_snps
7      BC  JOE -0.0002997228      0.002756867 0.03623437  42634
4      BC ECHO -0.0002287575      0.002599823 0.03662072  42592
18     AK   TL  0.0004880247      0.002376227 0.02936164  42556
21     AK   WT  0.0004924862      0.002365655 0.02939016  42649
16     BC  SWA -0.0001374537      0.002360428 0.03267996  42549
14     AK   SL  0.0006115779      0.002333856 0.02884917  42668
ΔAF variation is continuous, not driven by a single outlier population.




#gene divergent
suppressPackageStartupMessages({
  library(data.table)
})

# ==========================================
# 1. 路径设置 (Paths)
# ==========================================
DELTA_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
OUT_FILE     <- "/mnt/spareHD_2/nu_287/q2_parallelism/Mitonuclear_Target_Genes_Final.tsv"

# ==========================================
# 2. 数据读取与清洗 (Data Loading)
# ==========================================
normalize_pop <- function(x){
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
}

# 读取 ΔAF
DEL <- fread(cmd=paste("zcat", shQuote(DELTA_FILE)))
DEL[, pop := normalize_pop(pop)]

# 读取 Cluster
CL <- fread(CLUSTER_FILE, sep="\t", header=FALSE, fill=TRUE, strip.white=TRUE)
CL <- CL[nzchar(V1) & nzchar(V2), .(pop=V1, cluster=V2)]
CL <- CL[!(pop %in% c("pop", "POP", "sample"))] 
CL[, pop := normalize_pop(pop)]
CL <- unique(CL, by="pop")

# 合并
DEL <- merge(DEL, CL, by="pop", all.x=TRUE)
DEL <- DEL[!is.na(cluster)]

# ==========================================
# 3. 统计计算 (Statistics)
# ==========================================
# 计算每个基因在每个 Region、每个 Cluster 下的平均 |ΔAF|
gene_stat <- DEL[region == "BC", .(
  mean_abs_daf = mean(abs(deltaAF), na.rm = TRUE),
  n_snps = .N
), by = .(cluster, gene)]

# 转换为宽表
bc_wide <- dcast(gene_stat, gene + n_snps ~ cluster, value.var = "mean_abs_daf")
bc_wide <- bc_wide[n_snps >= 3] # 基础过滤

# 定义 Cluster 列名
cl_cols <- setdiff(names(bc_wide), c("gene", "n_snps"))

# 计算指标
bc_wide[, global_mean := rowMeans(.SD, na.rm=TRUE), .SDcols = cl_cols]
bc_wide[, cv := apply(.SD, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)), .SDcols = cl_cols]

# ==========================================
# 4. 双轨筛选逻辑 (Two-track Filtering)
# ==========================================

# --- 轨迹 A: 全谱系通用型 (Shared/Global Targets) ---
# 逻辑：分化极强且在不同 mt 背景下表现一致 (低 CV)
global_targets <- bc_wide[
  global_mean > quantile(global_mean, 0.95, na.rm=TRUE) & 
  cv < quantile(cv, 0.4, na.rm=TRUE)
][order(-global_mean)]
global_targets[, category := "Shared_Adaptive"]

# --- 轨迹 B: 谱系定制型 (Lineage-specific Targets) ---
# 逻辑：在特定 Cluster 表现极强且具有高特异性
specific_targets <- data.table()

for(cl in cl_cols){
  others <- setdiff(cl_cols, cl)
  spec_col <- paste0("spec_", cl)
  
  # 计算特异性倍数 (处理 NA)
  bc_wide[[spec_col]] <- apply(bc_wide, 1, function(row){
    target_val <- as.numeric(row[cl])
    other_vals <- as.numeric(row[others])
    other_vals <- other_vals[!is.na(other_vals)]
    if(is.na(target_val) || length(other_vals) == 0) return(NA_real_)
    return(target_val / mean(other_vals))
  })
  
  # 筛选
  res <- bc_wide[get(cl) > quantile(get(cl), 0.85, na.rm=TRUE) & get(spec_col) > 1.5]
  if(nrow(res) > 0){
    res_sub <- res[, .(gene, n_snps, intensity = get(cl), specificity = get(spec_col))]
    res_sub[, category := paste0("Specific_", cl)]
    specific_targets <- rbind(specific_targets, res_sub)
  }
}

# ==========================================
# 5. 合并并保存结果 (Export)
# ==========================================
# 整理最终表格
final_table <- rbind(
  global_targets[, .(gene, n_snps, score = global_mean, specificity = 1/cv, category)],
  specific_targets[, .(gene, n_snps, score = intensity, specificity, category)]
)

# 按照分数排序
final_table <- final_table[order(-score)]

write.table(final_table, OUT_FILE, sep="\t", row.names=FALSE, quote=FALSE)

cat("\nDone! 结果已保存至:", OUT_FILE, "\n")
print(final_table[gene %in% c("ndufa4", "cox4i1", "hccsb", "cox6a1")])


#

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 1. Paths
# ============================================================

# 改成你的 r2 < 0.2 LD-pruned deltaAF 文件
DELTA_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0_ldPruned/LM_perSNP_mtCluster_manual_plus_treePC12.ldPruned.tsv.gz"

CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"

OUT_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/Mitonuclear_Target_Genes_LDpruned_r02.tsv"

# 只分析 BC；如果要 AK，改成 "AK"
TARGET_REGION <- "BC"

# gene-level 最少 SNP 数
MIN_SNPS <- 5

# ============================================================
# 2. Functions
# ============================================================

normalize_pop <- function(x) {
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl = TRUE)
}

safe_cv <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(m) || m == 0) return(NA_real_)
  sd(x, na.rm = TRUE) / m
}

# ============================================================
# 3. Load deltaAF
# ============================================================

cat("[1] Loading deltaAF file...\n")

if (grepl("\\.gz$", DELTA_FILE)) {
  DEL <- fread(cmd = paste("zcat", shQuote(DELTA_FILE)))
} else {
  DEL <- fread(DELTA_FILE)
}

DEL[, pop := normalize_pop(pop)]

# 如果你的 LD-pruned 文件里 deltaAF 列名不同，在这里统一
if (!"deltaAF" %in% names(DEL)) {
  possible_cols <- c("delta_af", "daf", "DeltaAF", "delta")
  hit <- intersect(possible_cols, names(DEL))
  if (length(hit) == 1) {
    setnames(DEL, hit, "deltaAF")
  } else {
    stop("Cannot find deltaAF column. Please check column names.")
  }
}

needed_cols <- c("pop", "region", "gene", "deltaAF")
missing_cols <- setdiff(needed_cols, names(DEL))
if (length(missing_cols) > 0) {
  stop("Missing columns in DELTA_FILE: ", paste(missing_cols, collapse = ", "))
}

# ============================================================
# 4. Load mt cluster
# ============================================================

cat("[2] Loading mt cluster file...\n")

CL <- fread(CLUSTER_FILE, sep = "\t", header = FALSE, fill = TRUE, strip.white = TRUE)
CL <- CL[nzchar(V1) & nzchar(V2), .(pop = V1, cluster = V2)]
CL <- CL[!(pop %in% c("pop", "POP", "sample", "Sample"))]
CL[, pop := normalize_pop(pop)]
CL <- unique(CL, by = "pop")

# ============================================================
# 5. Merge
# ============================================================

cat("[3] Merging deltaAF with mt cluster...\n")

DEL <- merge(DEL, CL, by = "pop", all.x = TRUE)
DEL <- DEL[!is.na(cluster)]
DEL <- DEL[region == TARGET_REGION]
DEL <- DEL[is.finite(deltaAF)]

cat("Region:", TARGET_REGION, "\n")
cat("Number of rows:", nrow(DEL), "\n")
cat("Number of genes:", length(unique(DEL$gene)), "\n")
cat("Number of clusters:", length(unique(DEL$cluster)), "\n")

# ============================================================
# 6. Gene x cluster summary
# ============================================================

cat("[4] Calculating gene-level statistics...\n")

gene_stat <- DEL[, .(
  median_abs_daf = median(abs(deltaAF), na.rm = TRUE),
  mean_abs_daf   = mean(abs(deltaAF), na.rm = TRUE),
  n_snps_cluster = .N
), by = .(gene, cluster)]

# gene total SNP number
gene_n <- DEL[, .(
  n_snps = uniqueN(paste(chr, pos, sep = ":"))
), by = gene]

# 如果没有 chr/pos，就用行数
if (!all(c("chr", "pos") %in% names(DEL))) {
  gene_n <- DEL[, .(n_snps = .N), by = gene]
}

bc_wide <- dcast(
  gene_stat,
  gene ~ cluster,
  value.var = "median_abs_daf"
)

bc_wide <- merge(bc_wide, gene_n, by = "gene", all.x = TRUE)

# filter after LD pruning
bc_wide <- bc_wide[n_snps >= MIN_SNPS]

cl_cols <- setdiff(names(bc_wide), c("gene", "n_snps"))

cat("Genes retained after n_snps >=", MIN_SNPS, ":", nrow(bc_wide), "\n")

# ============================================================
# 7. Global/shared target statistics
# ============================================================

bc_wide[, global_score := rowMeans(.SD, na.rm = TRUE), .SDcols = cl_cols]

bc_wide[, cv := apply(.SD, 1, safe_cv), .SDcols = cl_cols]

# avoid Inf
bc_wide[, shared_specificity := 1 / (cv + 1e-6)]

# ============================================================
# 8. Track A: shared/global adaptive targets
# ============================================================

cat("[5] Identifying shared adaptive targets...\n")

global_cutoff <- quantile(bc_wide$global_score, 0.95, na.rm = TRUE)
cv_cutoff     <- quantile(bc_wide$cv, 0.40, na.rm = TRUE)

global_targets <- bc_wide[
  global_score >= global_cutoff &
    cv <= cv_cutoff
]

global_targets <- global_targets[, .(
  gene,
  n_snps,
  score = global_score,
  specificity = shared_specificity,
  best_cluster = NA_character_,
  category = "Shared_Adaptive"
)]

# ============================================================
# 9. Track B: lineage-specific targets
# ============================================================

cat("[6] Identifying lineage-specific targets...\n")

specific_targets <- data.table()

for (cl in cl_cols) {
  
  others <- setdiff(cl_cols, cl)
  spec_col <- paste0("specificity_", cl)
  
  bc_wide[, (spec_col) := {
    target_val <- get(cl)
    other_mean <- rowMeans(.SD, na.rm = TRUE)
    target_val / (other_mean + 1e-6)
  }, .SDcols = others]
  
  cluster_cutoff <- quantile(bc_wide[[cl]], 0.85, na.rm = TRUE)
  
  res <- bc_wide[
    get(cl) >= cluster_cutoff &
      get(spec_col) >= 1.5
  ]
  
  if (nrow(res) > 0) {
    res_sub <- res[, .(
      gene,
      n_snps,
      score = get(cl),
      specificity = get(spec_col),
      best_cluster = cl,
      category = paste0("Specific_", cl)
    )]
    
    specific_targets <- rbind(specific_targets, res_sub, fill = TRUE)
  }
}

# ============================================================
# 10. Combine final table
# ============================================================

cat("[7] Combining final results...\n")

final_table <- rbind(
  global_targets,
  specific_targets,
  fill = TRUE
)

final_table <- final_table[order(-score)]

# remove exact duplicates if any
final_table <- unique(final_table)

# ============================================================
# 11. Export
# ============================================================

fwrite(final_table, OUT_FILE, sep = "\t")

cat("\nDone!\n")
cat("Output saved to:\n", OUT_FILE, "\n\n")

cat("Top candidates:\n")
print(head(final_table, 30))

cat("\nCheck focal genes:\n")
print(final_table[gene %in% c(
  "ndufa4", "cox4i1", "hccsb", "cox6a1",
  "ndufb4", "ndufs5", "cox7a1", "ndufab1b"
)])


/mnt/spareHD_2/nu_287/q2_parallelism/Table6_gene_candidates_LDpruned.tsv 

    region     gene n_SNPs     Score Best_cluster Specificity
    <char>   <char>  <int>     <num>       <char>       <num>
 1:     AK   cox7a1      2 65.564988       C4_GOS   1.0000000
 2:     AK    hccsb     39 46.934795       C4_GOS   0.7179487
 3:     AK ndufab1b     17 34.905428       C4_GOS   0.7058824
 4:     AK   cox6a1     15 34.791622       C4_GOS   0.4666667
 5:     AK   ndufs5     10 33.110967       C4_GOS   0.5000000
 6:     AK   cox4i1     14 28.845777       C4_GOS   0.6428571
 7:     AK   ndufa4      7 19.151667        C1_AK   0.4285714
 8:     AK   ndufb4      9 15.841648        C1_AK   0.5555556
 9:     BC   cox6a1      1  9.710942    C2_Recent   1.0000000
10:     BC   ndufs5      1  3.316732    C2_Recent   1.0000000
11:     BC    hccsb      4  3.244737    C2_Recent   0.7500000
12:     BC   cox7a1      1  1.643157        C1_AK   1.0000000
               Pattern
                <char>
 1:    Specific_C4_GOS
 2:    Specific_C4_GOS
 3:    Specific_C4_GOS
 4:             Shared
 5:    Specific_C4_GOS
 6:    Specific_C4_GOS
 7:             Shared
 8:     Specific_C1_AK
 9: Specific_C2_Recent
10: Specific_C2_Recent
11: Specific_C2_Recent
12:     Specific_C1_AK


 /mnt/spareHD_2/nu_287/q2_parallelism/Table6_gene_candidates_LDpruned.tsv 