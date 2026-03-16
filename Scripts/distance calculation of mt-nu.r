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