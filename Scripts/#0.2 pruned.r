#0.2 r
#!/usr/bin/env python3

import os
import pandas as pd
from collections import defaultdict

# =========================
# 1. 参数
# =========================
ld_dir = "/work/cyu/ldx_all_subunits/ld"
delta_file = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"

# Bayenv-style pruning:
# keep SNPs with low LD; remove redundant SNPs connected by r2 >= 0.2
R2_THRESHOLD = 0.2

out_dir = "/work/cyu/ldx_all_subunits/ld/bayenv_style_r2_0.2"
os.makedirs(out_dir, exist_ok=True)

freshwater = {
    "BEA","BOOT","ECHO","FG","GOS","JOE","LAW","LB","LG",
    "MUC","PYE","ROB","SL","SR","SWA","THE","TL","WB","WK","WT"
}

cluster_map = {
    "SR":"C1_AK","TL":"C1_AK","WK":"C1_AK","LB":"C1_AK","MUC":"C1_AK","SWA":"C1_AK",
    "BEA":"C2_Recent","THE":"C2_Recent",
    "ROB":"C3_MarineLike","WB":"C3_MarineLike","LG":"C3_MarineLike","SL":"C3_MarineLike",
    "LAW":"C3_MarineLike","BOOT":"C3_MarineLike","JOE":"C3_MarineLike","FG":"C3_MarineLike",
    "GOS":"C4_GOS","PYE":"C4_GOS","ECHO":"C4_GOS","WT":"C4_GOS"
}

# =========================
# 2. 读入 deltaAF 表
# =========================
df = pd.read_csv(delta_file, sep="\t", compression="gzip")

df["pop"] = df["pop"].astype(str).str.upper()
df["pos"] = df["pos"].astype(str)

df = df[df["pop"].isin(freshwater)].copy()

val_map = df.set_index(["pop", "pos"])["deltaAF"].abs().to_dict()
delta_pos_set = set(zip(df["pop"], df["pos"]))

# =========================
# 3. 容器
# =========================
block_rows = []
summary_rows = []

redundant_pairs = set()
leader_pairs = set()

# =========================
# 4. 遍历 LDx .ld 文件
# =========================
for ld_file in os.listdir(ld_dir):

    if not ld_file.endswith(".ld"):
        continue

    prefix = ld_file.replace(".ld", "")

    if "__" not in prefix:
        continue

    pop, gene = prefix.split("__", 1)
    pop = pop.upper()

    if pop not in freshwater:
        continue

    graph = defaultdict(set)
    nodes = set()

    with open(os.path.join(ld_dir, ld_file)) as f:
        for line in f:
            cols = line.strip().split()

            if len(cols) < 15:
                continue

            p1, p2 = cols[0], cols[1]

            try:
                r2 = float(cols[14])
            except ValueError:
                continue

            # Bayenv-style LD threshold
            if r2 >= R2_THRESHOLD:
                graph[p1].add(p2)
                graph[p2].add(p1)
                nodes.update([p1, p2])

    # =========================
    # 5. connected components = LD blocks
    # =========================
    visited = set()
    blocks = []

    for node in nodes:
        if node in visited:
            continue

        block = []
        stack = [node]

        while stack:
            curr = stack.pop()

            if curr in visited:
                continue

            visited.add(curr)
            block.append(curr)
            stack.extend(graph[curr] - visited)

        blocks.append(sorted(block, key=lambda x: int(x)))

    # =========================
    # 6. 每个 LD block 选 leader SNP
    # =========================
    for i, block in enumerate(blocks):

        # 只 pruning 出现在 deltaAF 表里的 SNP
        block_hits = [x for x in block if (pop, x) in delta_pos_set]

        if len(block_hits) > 0:
            # leader = 最大 |deltaAF| SNP
            leader = sorted(
                block_hits,
                key=lambda x: (-val_map[(pop, x)], int(x))
            )[0]

            leader_source = "max_abs_deltaAF"
            leader_deltaAF = val_map[(pop, leader)]
            leader_has_deltaAF = True

            leader_pairs.add((pop, leader))

            # 其余 deltaAF SNP 标记为 redundant
            for snp_pos in block_hits:
                if snp_pos != leader:
                    redundant_pairs.add((pop, snp_pos))

        else:
            leader = block[0]
            leader_source = "min_pos_no_deltaAF_overlap"
            leader_deltaAF = None
            leader_has_deltaAF = False

        block_rows.append({
            "Pop": pop,
            "Gene": gene,
            "Cluster": cluster_map.get(pop, "NA"),
            "Block_ID": f"{pop}_{gene}_{i}",
            "Block_size": len(block),
            "n_deltaAF_snps_in_block": len(block_hits),
            "block_has_deltaAF": len(block_hits) > 0,
            "Leader_SNP": leader,
            "Leader_source": leader_source,
            "Leader_deltaAF": leader_deltaAF,
            "Leader_has_deltaAF": leader_has_deltaAF
        })

    summary_rows.append({
        "Pop": pop,
        "Gene": gene,
        "Cluster": cluster_map.get(pop, "NA"),
        "n_blocks": len(blocks),
        "total_snps_in_blocks": sum(len(b) for b in blocks) if blocks else 0,
        "mean_block_size": sum(len(b) for b in blocks) / len(blocks) if blocks else 0,
        "max_block_size": max(len(b) for b in blocks) if blocks else 0
    })

# =========================
# 7. 保存 block 信息
# =========================
block_df = pd.DataFrame(block_rows)
summary_df = pd.DataFrame(summary_rows)

block_out = os.path.join(out_dir, "ld_blocks_bayenvStyle_r2_0.2.csv")
summary_out = os.path.join(out_dir, "ld_block_pop_gene_summary_bayenvStyle_r2_0.2.csv")

block_df.to_csv(block_out, index=False)
summary_df.to_csv(summary_out, index=False)

# =========================
# 8. 对 deltaAF 表 pruning
# =========================
df["is_ld_leader"] = df.apply(lambda r: (r["pop"], r["pos"]) in leader_pairs, axis=1)
df["is_ld_redundant"] = df.apply(lambda r: (r["pop"], r["pos"]) in redundant_pairs, axis=1)

df_keep = df[~df["is_ld_redundant"]].copy()

masked_out = os.path.join(out_dir, "deltaAF_long.noAMO.bayenvStyle_r2_0.2_ldPruned_masked.tsv.gz")
kept_out = os.path.join(out_dir, "deltaAF_long.noAMO.bayenvStyle_r2_0.2_ldPruned_kept.tsv.gz")

df.to_csv(masked_out, sep="\t", index=False, compression="gzip")
df_keep.to_csv(kept_out, sep="\t", index=False, compression="gzip")

# =========================
# 9. Summary
# =========================
print("✅ Bayenv-style LD pruning complete")
print("R2 threshold:", R2_THRESHOLD)
print("Output dir:", out_dir)
print("Total LD blocks:", len(block_df))
print("Blocks with >=1 deltaAF SNP:", int(block_df["block_has_deltaAF"].sum()) if len(block_df) else 0)
print("Leader deltaAF SNPs kept:", int(df["is_ld_leader"].sum()))
print("Redundant deltaAF SNPs removed:", int(df["is_ld_redundant"].sum()))
print("Rows in original deltaAF table:", len(df))
print("Rows in pruned-kept table:", len(df_keep))
print("Saved:", kept_out)



#
library(data.table)

IN_LM <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm/LM_perSNP_mtCluster_plus_treePC12.LDpruned_r2_0.2.tsv.gz"

LM <- fread(cmd = paste("zcat", shQuote(IN_LM)))

# ----------------------------
# summary（和你旧的一样）
# ----------------------------
SUM <- LM[, .(
  n_snps = .N,
  prop_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
  prop_p10 = mean(p_cluster < 0.10, na.rm=TRUE),
  min_p = suppressWarnings(min(p_cluster, na.rm=TRUE)),
  med_p = suppressWarnings(median(p_cluster, na.rm=TRUE))
), by=region][order(region)]

cat("\n=== summary ===\n")
print(SUM)

# ----------------------------
# driver cluster（和旧的一样）
# ----------------------------
SIG <- LM[is.finite(p_cluster) & p_cluster < 0.05]

DRIVER <- SIG[, .N, by=.(region, driver_cluster)]
DRIVER[, prop := N / sum(N), by=region]
setorder(DRIVER, region, -prop)

cat("\n=== driver cluster share among significant SNPs p<0.05 ===\n")
print(DRIVER)


cat 09_mtCluster_deltaAF_SNPlevel_LM_perm_LDpruned_r2_0.2.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
})

# ----------------------------
# Config
# ----------------------------
DELTA_FILE <- Sys.getenv(
  "DELTA",
  unset = "/work/cyu/ldx_all_subunits/ld/bayenv_style_r2_0.2/deltaAF_long.noAMO.bayenvStyle_r2_0.2_ldPruned_kept.tsv.gz"
)

OUT_DIR <- Sys.getenv(
  "OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm"
)

MT_TREE <- Sys.getenv(
  "MT_TREE",
  unset = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_noDloop_noAMO.iqtree.treefile"
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

lm_cluster_Ftest <- function(dt){
  dt <- dt[
    is.finite(deltaAF) &
      !is.na(cluster) &
      is.finite(treePC1) &
      is.finite(treePC2)
  ]

  if(nrow(dt) < 6) {
    return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))
  }

  if(length(unique(dt$cluster)) < 2) {
    return(list(p=NA_real_, F=NA_real_, df1=NA_real_, df2=NA_real_, n=nrow(dt)))
  }

  fit_full <- try(lm(deltaAF ~ cluster + treePC1 + treePC2, data=dt), silent=TRUE)
  fit_red  <- try(lm(deltaAF ~ treePC1 + treePC2, data=dt), silent=TRUE)

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

permute_clusters_within_region <- function(pop_cluster_dt, region_name){
  sub <- pop_cluster_dt[region == region_name]
  if(nrow(sub) == 0) {
    return(data.table(pop=character(), cluster_perm=factor()))
  }

  sub[, cluster_perm := sample(cluster)]
  sub[, .(pop, cluster_perm)]
}

# ----------------------------
# Load LD-pruned deltaAF
# ----------------------------
cat("[read] LD-pruned deltaAF: ", DELTA_FILE, "\n", sep="")
DEL <- safe_zread(DELTA_FILE)

req <- c("region","snp","chr","pos","gene","pop","deltaAF","treePC1","treePC2")
miss <- setdiff(req, names(DEL))
if(length(miss) > 0){
  stop("Missing columns in DELTA_FILE: ", paste(miss, collapse=", "))
}

DEL[, pop := normalize_pop(pop)]
DEL <- DEL[pop != "AMO"]
DEL <- DEL[is.finite(deltaAF) & is.finite(treePC1) & is.finite(treePC2)]

cat("[info] rows after filter: ", nrow(DEL), "\n", sep="")
cat("[info] unique SNP-pop observations: ", uniqueN(DEL[, .(region, snp, pop)]), "\n", sep="")

# ----------------------------
# Load mt tree and define tree-based clusters
# ----------------------------
if(!file.exists(MT_TREE)){
  stop("MT_TREE not found: ", MT_TREE)
}

cat("[read tree] ", MT_TREE, "\n", sep="")
tr <- read.tree(MT_TREE)

D <- cophenetic(tr)
hc <- hclust(as.dist(D), method="average")
cl <- cutree(hc, k=K_CLUST)

mt_cluster <- data.table(
  pop = normalize_pop(names(cl)),
  mt_cluster = as.integer(cl)
)
mt_cluster <- unique(mt_cluster, by="pop")
mt_cluster <- mt_cluster[pop %in% unique(DEL$pop)]

DEL <- merge(DEL, mt_cluster, by="pop", all.x=TRUE)

cat("[info] pops missing mt_cluster: ", sum(is.na(DEL$mt_cluster)), "\n", sep="")
if(sum(is.na(DEL$mt_cluster)) > 0){
  print(unique(DEL[is.na(mt_cluster), pop]))
}

DEL <- DEL[!is.na(mt_cluster)]
DEL[, cluster := factor(paste0("C", mt_cluster))]

pop_region <- unique(DEL[, .(pop, region)], by="pop")
pop_cluster_out <- merge(
  pop_region,
  unique(DEL[, .(pop, mt_cluster)], by="pop"),
  by="pop",
  all.x=TRUE
)

fwrite(
  pop_cluster_out[order(region, mt_cluster, pop)],
  file.path(OUT_DIR, "mtCluster_by_pop.LDpruned.noAMO.tsv"),
  sep="\t"
)

cat("[info] pops per cluster:\n")
print(unique(DEL[, .(pop, region, cluster)])[order(region, cluster, pop)])

# ----------------------------
# SNP-level LM
# ----------------------------
cat("\n[info] rows by region:\n")
print(DEL[, .N, by=region][order(region)])

setkey(DEL, region, snp)

regions <- sort(unique(DEL$region))
snps_by_region <- lapply(regions, function(r) unique(DEL[region == r, snp]))
names(snps_by_region) <- regions

cat("\n[info] SNPs per region:\n")
for(r in regions){
  cat("  ", r, ": ", length(snps_by_region[[r]]), "\n", sep="")
}

lm_list <- vector("list", 100000)
idx <- 0L

for(r in regions){
  cat("[LM] region: ", r, "\n", sep="")
  snps <- snps_by_region[[r]]

  for(s in snps){
    dt <- DEL[list(r, s)]

    out <- lm_cluster_Ftest(dt)

    cm <- dt[, .(
      mu = mean(deltaAF, na.rm=TRUE),
      n = .N
    ), by=cluster][order(cluster)]

    cm[, abs_mu := abs(mu)]
    driver_i <- which.max(cm$abs_mu)

    driver <- if(length(driver_i) == 1) as.character(cm$cluster[driver_i]) else NA_character_
    driver_mu <- if(length(driver_i) == 1) cm$mu[driver_i] else NA_real_

    cm_str <- paste0(
      cm$cluster, ":",
      sprintf("%.6g", cm$mu),
      "(n=", cm$n, ")",
      collapse=";"
    )

    idx <- idx + 1L
    lm_list[[idx]] <- data.table(
      region = r,
      snp = s,
      chr = dt$chr[1],
      pos = dt$pos[1],
      gene = dt$gene[1],
      n = out$n,
      n_clusters = length(unique(dt$cluster)),
      F_cluster = out$F,
      df1 = out$df1,
      df2 = out$df2,
      p_cluster = out$p,
      driver_cluster = driver,
      driver_mu = driver_mu,
      cluster_means = cm_str
    )
  }
}

LM <- rbindlist(lm_list[seq_len(idx)], use.names=TRUE, fill=TRUE)

LM[, q_cluster := p.adjust(p_cluster, method="BH"), by=region]

lm_out <- file.path(OUT_DIR, "LM_perSNP_mtCluster_plus_treePC12.LDpruned_r2_0.2.tsv.gz")
fwrite(LM, lm_out, sep="\t")

cat("\n[write] ", lm_out, "\n", sep="")

SUM <- LM[, .(
  n_snps = .N,
  prop_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
  prop_p10 = mean(p_cluster < 0.10, na.rm=TRUE),
  prop_q10 = mean(q_cluster < 0.10, na.rm=TRUE),
  min_p = suppressWarnings(min(p_cluster, na.rm=TRUE)),
  min_q = suppressWarnings(min(q_cluster, na.rm=TRUE)),
  med_p = suppressWarnings(median(p_cluster, na.rm=TRUE))
), by=region][order(region)]

fwrite(SUM, file.path(OUT_DIR, "LM_cluster_summary_by_region.tsv"), sep="\t")

cat("\n=== LM summary ===\n")
print(SUM)

SIG <- LM[is.finite(p_cluster) & p_cluster < 0.05]
DRIVER <- SIG[, .N, by=.(region, driver_cluster)]
DRIVER[, prop := N / sum(N), by=region]
setorder(DRIVER, region, -prop)

fwrite(DRIVER, file.path(OUT_DIR, "LM_driverCluster_share_sigSNP.tsv"), sep="\t")

cat("\n=== driver cluster share among p<0.05 SNPs ===\n")
print(DRIVER)

# ----------------------------
# Permutation
# ----------------------------
perm_out_file <- file.path(OUT_DIR, "Perm_perSNP_mtCluster_plus_treePC12.LDpruned_r2_0.2.tsv.gz")
perm_global_file <- file.path(OUT_DIR, "Perm_global_summary.tsv")

if(N_PERM <= 0){

  cat("\n[perm] N_PERM=0, skip permutation.\n")

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

    if(b %% 25 == 0) {
      cat("[perm] ", b, "/", N_PERM, "\n", sep="")
    }

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

        row_i <- which(OBS$region == r & OBS$snp == s)

        if(length(row_i) == 1){
          OBS$n_perm_used[row_i] <- OBS$n_perm_used[row_i] + 1L
          if(out$F >= f_obs){
            OBS$ge_count[row_i] <- OBS$ge_count[row_i] + 1L
          }
        }
      }
    }
  }

  OBS[, p_perm := (1 + ge_count) / (1 + n_perm_used)]

  PERM <- merge(
    LM,
    OBS[, .(region, snp, n_perm_used, p_perm)],
    by=c("region","snp"),
    all.x=TRUE
  )

  PERM[, q_perm := p.adjust(p_perm, method="BH"), by=region]

  fwrite(PERM, perm_out_file, sep="\t")

  cat("\n[write] ", perm_out_file, "\n", sep="")

  G <- PERM[, .(
    n_snps = .N,
    prop_LM_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
    prop_LM_q10 = mean(q_cluster < 0.10, na.rm=TRUE),
    prop_perm_p05 = mean(p_perm < 0.05, na.rm=TRUE),
    prop_perm_q10 = mean(q_perm < 0.10, na.rm=TRUE),
    min_perm_p = suppressWarnings(min(p_perm, na.rm=TRUE)),
    min_perm_q = suppressWarnings(min(q_perm, na.rm=TRUE)),
    med_perm_p = suppressWarnings(median(p_perm, na.rm=TRUE)),
    n_perm_used_med = suppressWarnings(median(n_perm_used, na.rm=TRUE))
  ), by=region][order(region)]

  fwrite(G, perm_global_file, sep="\t")

  cat("\n=== permutation summary ===\n")
  print(G)
}

# ----------------------------
# Plots
# ----------------------------
p_hist <- ggplot(LM[is.finite(p_cluster)], aes(p_cluster)) +
  geom_histogram(bins=50) +
  facet_wrap(~region, ncol=1, scales="free_y") +
  labs(
    x="LM p-value for mtCluster term",
    y="count",
    title="LD-pruned SNP-level association: deltaAF ~ mtCluster + treePC1 + treePC2"
  ) +
  theme_classic(base_size=14)

ggsave(
  file.path(OUT_DIR, "Fig_pHist_mtCluster_LM_LDpruned_r2_0.2.png"),
  p_hist,
  width=7.2,
  height=6.5,
  dpi=300
)

LMq <- LM[is.finite(p_cluster) & p_cluster > 0 & p_cluster <= 1]

if(nrow(LMq) > 0){
  QQ <- LMq[, .(p_cluster = sort(p_cluster)), by=region]
  QQ[, obs := -log10(p_cluster), by=region]
  QQ[, exp := -log10(ppoints(.N)), by=region]

  p_qq <- ggplot(QQ, aes(exp, obs)) +
    geom_point(alpha=0.4, size=0.8) +
    geom_abline(slope=1, intercept=0) +
    facet_wrap(~region) +
    labs(
      x="Expected -log10(p)",
      y="Observed -log10(p)",
      title="QQ plot: LD-pruned mtCluster p-values"
    ) +
    theme_classic(base_size=14)

  ggsave(
    file.path(OUT_DIR, "Fig_QQ_mtCluster_LM_LDpruned_r2_0.2.png"),
    p_qq,
    width=7.2,
    height=4.2,
    dpi=300
  )
}

if(N_PERM > 0 && file.exists(perm_out_file)){
  PERM <- safe_zread(perm_out_file)

  p_hist2 <- ggplot(PERM[is.finite(p_perm)], aes(p_perm)) +
    geom_histogram(bins=50) +
    facet_wrap(~region, ncol=1, scales="free_y") +
    labs(
      x="Permutation p-value for mtCluster term",
      y="count",
      title=paste0("LD-pruned permutation test N_PERM=", N_PERM)
    ) +
    theme_classic(base_size=14)

  ggsave(
    file.path(OUT_DIR, "Fig_pHist_mtCluster_perm_LDpruned_r2_0.2.png"),
    p_hist2,
    width=7.2,
    height=6.5,
    dpi=300
  )
}

cat("\n=== Top 20 SNPs by p_cluster per region ===\n")
LM_top <- LM[is.finite(p_cluster)][order(p_cluster)][, head(.SD, 20), by=region]
print(LM_top[, .(region, snp, gene, chr, pos, n, n_clusters, F_cluster, p_cluster, q_cluster, driver_cluster)])

cat("\n[OK] done. Outputs in:\n", OUT_DIR, "\n", sep="")







library(data.table)
library(ggplot2)

INFILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm/gene_enrichment/Gene_enrichment_overall.tsv.gz"

OUTFIG <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm/gene_enrichment/Fig_gene_enrichment_heatmap_OR.png"

G <- fread(INFILE)

# log2 OR
G[, log2OR := log2(or_overall)]
G[!is.finite(log2OR), log2OR := NA_real_]

# significance label
G[, sig_label := fifelse(q_overall < 0.1, "+",
                  
                                 fifelse(p_overall < 0.05, "*", ""))]

# top genes
top_genes <- unique(G[order(p_overall), head(gene, 20), by=region]$V1)
P <- G[gene %in% top_genes]

# order genes
gene_order <- P[, .(best_p = min(p_overall, na.rm=TRUE)), by=gene][order(best_p)]$gene
P[, gene := factor(gene, levels=rev(gene_order))]

p <- ggplot(P, aes(x=region, y=gene, fill=log2OR)) +
  geom_tile(color="white") +
  geom_text(aes(label=sig_label), size=5) +
  scale_fill_gradient2(
    low="blue",
    mid="white",
    high="red",
    midpoint=0,
    name="log2(OR)"
  ) +
  theme_classic(base_size=14) +
  labs(
    x=NULL,
    y="Gene",
    title="Gene enrichment heatmap",
    subtitle="* p<0.05, + q<0.10 "
  )

ggsave(OUTFIG, p, width=5, height=8, dpi=300)

cat("saved:", OUTFIG, "\n")
print(p)



= summary === region n_snps prop_p05 prop_p10 min_p med_p <char> <int> <num> <num> <num> <num> 1: AK 42733 0.06877605 0.07079702 5.009329e-70 0.7044233 2: BC 42695 0.01165226 0.01872218 4.832982e-64 0.2948301 === driver cluster share among significant SNPs p<0.05 === region driver_cluster N prop <char> <char> <int> <num> 1: AK C4_GOS 572 0.52525253 2: AK C1_AK 364 0.33425161 3: AK C3_MarineLike 153 0.14049587 4: BC C2_Recent 61 0.68539326 5: BC C1_AK 18 0.20224719 6: BC C3_MarineLike 6 0.06741573 7: BC C4_GOS 4 0.0449438





#vol plot

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

IN_LM <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm/LM_perSNP_mtCluster_plus_treePC12.LDpruned_r2_0.2.tsv.gz"

OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_treeBased_bayenvStyle_r2_0.2_perm/figures_mtCluster"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

while (dev.cur() > 1) dev.off()

LM <- fread(cmd = paste("zcat", shQuote(IN_LM)))

LM <- LM[
  is.finite(p_cluster) &
    p_cluster > 0 &
    p_cluster <= 1
]

LM[, logp := -log10(p_cluster)]
LM[, pos := as.numeric(pos)]
LM[, driver_mu := as.numeric(driver_mu)]
LM[, F_cluster := as.numeric(F_cluster)]

LM <- LM[
  !is.na(chr) &
    is.finite(pos) &
    is.finite(logp)
]

cat("Rows used:", nrow(LM), "\n")
cat("Regions:", paste(unique(LM$region), collapse = ", "), "\n")

# =========================
# Manhattan
# =========================

chr_order <- unique(LM[order(chr), chr])
LM[, chr := factor(chr, levels = chr_order)]

chr_info <- LM[, .(
  chr_len = max(pos, na.rm = TRUE)
), by = chr][order(chr)]

chr_info[, offset := shift(cumsum(chr_len), fill = 0)]

LM <- merge(
  LM,
  chr_info[, .(chr, offset)],
  by = "chr",
  all.x = TRUE
)

LM[, pos_cum := pos + offset]

axis_df <- LM[, .(
  center = mean(range(pos_cum, na.rm = TRUE))
), by = chr][order(chr)]

TOP_MAN <- LM[order(p_cluster), head(.SD, 10), by = region]
TOP_MAN[, label := gene]

p_manhattan <- ggplot(
  LM,
  aes(x = pos_cum, y = logp, color = chr)
) +
  geom_point(alpha = 0.65, size = 0.7) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.10),
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = TOP_MAN,
    aes(label = label),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.3,
    min.segment.length = 0
  ) +
  facet_wrap(~ region, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = as.character(axis_df$chr),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_color_manual(
    values = rep(c("grey25", "grey65"), length.out = length(chr_order))
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](p)),
    title = "mtCluster association Manhattan plot",
    subtitle = expression(delta*AF~"~ mtCluster + treePC1 + treePC2")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(OUTDIR, "Fig_mtCluster_association_Manhattan.png"),
  p_manhattan,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig_mtCluster_association_Manhattan.pdf"),
  p_manhattan,
  width = 12,
  height = 7
)

# =========================
# Volcano / association-strength plot
# =========================

VOL <- LM[
  is.finite(F_cluster) &
    is.finite(logp)
]

VOL[, sig := fifelse(p_cluster < 0.05, "p < 0.05", "NS")]

# label strongest non-duplicated genes per region
TOP_VOL <- VOL[logp > 5]
TOP_VOL <- TOP_VOL[order(region, -logp)]
TOP_VOL <- TOP_VOL[, .SD[!duplicated(gene)], by = region]
TOP_VOL <- TOP_VOL[, head(.SD, 8), by = region]

p_volcano_F <- ggplot(
  VOL,
  aes(x = F_cluster, y = logp)
) +
  geom_point(
    aes(color = driver_cluster, shape = sig),
    alpha = 0.75,
    size = 1.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = TOP_VOL,
    aes(label = gene),
    size = 3.2,
    max.overlaps = 30,
    box.padding = 0.35,
    min.segment.length = 0
  ) +
  facet_wrap(~ region, scales = "free_x") +
  coord_cartesian(ylim = c(0, 15)) +
  labs(
    x = "F statistic for mtCluster association",
    y = expression(-log[10](p)),
    color = "Driver cluster",
    shape = NULL,
    title = "mtCluster-associated nuclear SNPs",
    subtitle = expression(delta*AF~"~ mtCluster + treePC1 + treePC2")
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(OUTDIR, "Fig_mtCluster_association_Fstat_Volcano_cappedY15.png"),
  p_volcano_F,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(OUTDIR, "Fig_mtCluster_association_Fstat_Volcano_cappedY15.pdf"),
  p_volcano_F,
  width = 9,
  height = 5
)

# =========================
# Optional old driver_mu volcano
# =========================

VOL2 <- LM[
  is.finite(driver_mu) &
    is.finite(logp)
]

VOL2[, sig := fifelse(p_cluster < 0.05, "p < 0.05", "NS")]

TOP_VOL2 <- VOL2[logp > 5]
TOP_VOL2 <- TOP_VOL2[order(region, -logp)]
TOP_VOL2 <- TOP_VOL2[, .SD[!duplicated(gene)], by = region]
TOP_VOL2 <- TOP_VOL2[, head(.SD, 8), by = region]

p_volcano_mu <- ggplot(
  VOL2,
  aes(x = driver_mu, y = logp)
) +
  geom_point(
    aes(color = driver_cluster, shape = sig),
    alpha = 0.75,
    size = 1.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = TOP_VOL2,
    aes(label = gene),
    size = 3.2,
    max.overlaps = 30,
    box.padding = 0.35,
    min.segment.length = 0
  ) +
  facet_wrap(~ region, scales = "free") +
  coord_cartesian(ylim = c(0, 15)) +
  labs(
    x = "Driver-cluster mean ΔAF",
    y = expression(-log[10](p)),
    color = "Driver cluster",
    shape = NULL,
    title = "mtCluster association volcano plot",
    subtitle = "Y-axis capped at 15 for visualization"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(OUTDIR, "Fig_mtCluster_association_driverMu_Volcano_cappedY15.png"),
  p_volcano_mu,
  width = 9,
  height = 5,
  dpi = 300
)

# =========================
# Top table
# =========================

TOP_TABLE <- LM[order(p_cluster), head(.SD, 20), by = region]

fwrite(
  TOP_TABLE[, .(
    region, snp, gene, chr, pos,
    n, n_clusters,
    F_cluster, p_cluster, q_cluster,
    driver_cluster, driver_mu
  )],
  file.path(OUTDIR, "Top20_mtCluster_association_SNPs_by_region.tsv"),
  sep = "\t"
)

cat("\nDone. Files saved in:\n")
cat(OUTDIR, "\n")

cat("\nCheck with:\n")
cat("ls -lh ", OUTDIR, "\n", sep = "")


