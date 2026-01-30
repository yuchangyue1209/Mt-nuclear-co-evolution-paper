#q2 enviornment factors
/mnt/spareHD_2/nu_287/oldsig_subunit_SNP_rows_plusSI.tsv
oldsig_subunit_gene_table.tsv



#287 跑
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ========= 输入 =========
af_file  <- "/mnt/spareHD_2/nu_287/af_long.tsv.gz"          # chr pos gene pop af depth
cov_file <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"   # pop region mitoPC* treePC*
out_dir  <- "/mnt/spareHD_2/nu_287/q2_parallelism"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========= 参数 =========
min_depth <- 10        # 每个pop的最低depth过滤（你也可以设 5 或 20）
min_n_fresh_AK <- 5    # AK freshwater 至少多少个群体有数据
min_n_fresh_BC <- 8    # BC freshwater 至少多少个群体有数据

# ========= 你的分组（按你固定那套）=========
AK_fresh <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_fresh <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")
AK_marine <- "RS"
BC_marine <- "SAY"

# ========= 读数据 =========
AF <- fread(cmd = paste("zcat", shQuote(af_file)))
COV <- fread(cov_file)

# pop 名标准化：10_THE_S10 -> THE
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}
AF[, pop := normalize_pop(pop)]
COV[, pop := toupper(pop)]

# depth过滤
AF <- AF[depth >= min_depth]
AF[, snp := paste(chr, pos, gene, sep=":")]

# 只保留我们关心的 pops（AK freshwater + RS + BC freshwater + SAY）
keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]

# ========= helper：给定一个 SNP，算某组 pops 的均值 AF（并记录n）=========
mean_af <- function(dt, pops){
  sub <- dt[pop %in% pops]
  if(nrow(sub)==0) return(list(mu=NA_real_, n=0L))
  list(mu=mean(sub$af, na.rm=TRUE), n=as.integer(sum(!is.na(sub$af))))
}

# ========= SNP层面的 ΔAF =========
snps <- unique(AF$snp)

res <- vector("list", length(snps))
names(res) <- snps

for(i in seq_along(snps)){
  s <- snps[i]
  dt <- AF[snp == s, .(pop, af)]
  # marine
  mu_RS  <- mean_af(dt, AK_marine)
  mu_SAY <- mean_af(dt, BC_marine)
  # freshwater mean
  mu_AK  <- mean_af(dt, AK_fresh)
  mu_BC  <- mean_af(dt, BC_fresh)

  if(mu_RS$n < 1 || mu_SAY$n < 1) next
  if(mu_AK$n < min_n_fresh_AK || mu_BC$n < min_n_fresh_BC) next

  dAK <- mu_AK$mu - mu_RS$mu
  dBC <- mu_BC$mu - mu_SAY$mu

  # 并行性指标
  concordant <- sign(dAK) == sign(dBC) && sign(dAK)!=0 && sign(dBC)!=0
  prodv <- dAK * dBC   # >0 同向；<0 反向

  # gene/pos拆出来
  parts <- strsplit(s, ":", fixed=TRUE)[[1]]
  chr_ <- parts[1]; pos_ <- as.integer(parts[2]); gene_ <- parts[3]

  res[[i]] <- data.frame(
    chr=chr_, pos=pos_, gene=gene_, snp=s,
    AF_RS=mu_RS$mu,  n_RS=mu_RS$n,
    AF_SAY=mu_SAY$mu, n_SAY=mu_SAY$n,
    AF_AKfresh=mu_AK$mu, n_AKfresh=mu_AK$n,
    AF_BCfresh=mu_BC$mu, n_BCfresh=mu_BC$n,
    dAF_AK=dAK, dAF_BC=dBC,
    prod=prodv,
    concordant=concordant
  )
}

SNP <- bind_rows(res)
SNP <- SNP[is.finite(SNP$dAF_AK) & is.finite(SNP$dAF_BC), ]

# ========= 全局并行性统计 =========
# 1) 相关（越强越支持“环境驱动的平行响应”）
cor_pearson <- cor(SNP$dAF_AK, SNP$dAF_BC, method="pearson")
cor_spearman <- cor(SNP$dAF_AK, SNP$dAF_BC, method="spearman")

# 2) 同向比例 vs 50%（符号检验：二项检验）
ok <- SNP$dAF_AK!=0 & SNP$dAF_BC!=0
k <- sum(SNP$concordant[ok])
n <- sum(ok)
binom_p <- binom.test(k, n, p=0.5)$p.value

stats <- data.frame(
  n_snps=nrow(SNP),
  n_nonzero=n,
  concordant=k,
  concordant_prop=k/n,
  cor_pearson=cor_pearson,
  cor_spearman=cor_spearman,
  binom_p=binom_p
)

# ========= 输出 =========
fwrite(SNP, file.path(out_dir, "q2_parallelism_perSNP.tsv.gz"), sep="\t")
fwrite(stats, file.path(out_dir, "q2_parallelism_summary.tsv"), sep="\t")

# ========= 图：ΔAF散点 =========
p1 <- ggplot(SNP, aes(dAF_AK, dAF_BC)) +
  geom_point(alpha=0.35, size=0.8) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_smooth(method="lm", se=FALSE) +
  labs(
    x=expression(Delta~AF[AK]~" (AK freshwater mean - RS)"),
    y=expression(Delta~AF[BC]~" (BC freshwater mean - SAY)"),
    title="Q2 parallelism: per-SNP ΔAF concordance across AK & BC"
  )

ggsave(file.path(out_dir, "q2_parallelism_scatter_perSNP.png"), p1, width=6.5, height=5, dpi=300)

# ========= 图：prod分布（>0 同向）=========
p2 <- ggplot(SNP, aes(prod)) +
  geom_histogram(bins=80) +
  geom_vline(xintercept=0) +
  labs(x="ΔAF_AK * ΔAF_BC", title="Positive values indicate same-direction shifts")

ggsave(file.path(out_dir, "q2_parallelism_prod_hist.png"), p2, width=6.5, height=4.5, dpi=300)

# ========= gene层面汇总：每基因用 SNP ΔAF 的中位数（更稳健）=========
GENE <- SNP %>%
  group_by(gene) %>%
  summarise(
    n_snps=n(),
    dAK_med=median(dAF_AK, na.rm=TRUE),
    dBC_med=median(dAF_BC, na.rm=TRUE),
    concordant_prop=mean(concordant[dAF_AK!=0 & dAF_BC!=0], na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(prod_med = dAK_med * dBC_med)

fwrite(GENE, file.path(out_dir, "q2_parallelism_perGene.tsv"), sep="\t")

p3 <- ggplot(GENE, aes(dAK_med, dBC_med)) +
  geom_point(alpha=0.7) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_smooth(method="lm", se=FALSE) +
  labs(
    x="gene median ΔAF_AK",
    y="gene median ΔAF_BC",
    title="Q2 parallelism: gene-level median ΔAF"
  )

ggsave(file.path(out_dir, "q2_parallelism_scatter_perGene.png"), p3, width=6.5, height=5, dpi=300)

cat("[OK] wrote outputs to:", out_dir, "\n")
print(stats)


[OK] wrote outputs to: /mnt/spareHD_2/nu_287/q2_parallelism n_snps n_nonzero concordant concordant_prop cor_pearson cor_spearman 1 226456 21253 7324 0.3446102 0.5635678 -0.01660664 binom_p 1 4.940656e-324


#q2 step1 使用 subunit
#build subunit only af tsv
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ========= 输入 =========
af_file   <- "/mnt/spareHD_2/nu_287/af_long.tsv.gz"          
role_file <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/codeml_sites_summary.merged.tsv"
out_dir   <- "/mnt/spareHD_2/nu_287/q2_parallelism"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_af   <- file.path(out_dir, "af_long_subunit_plus_si.tsv.gz")
out_sum  <- file.path(out_dir, "subunit_plus_si_coverage_summary.tsv")
out_si   <- file.path(out_dir, "si_genes_detected_from_af.txt")
out_miss <- file.path(out_dir, "missing_genes_in_AF.subunit_plus_si.txt")

# ========= 1) subunit genes =========
role_dt <- fread(role_file)
subunit_genes <- unique(role_dt[role == "subunit" & !is.na(gene) & gene != "", gene])
cat("[INFO] subunit genes in role table:", length(subunit_genes), "\n")

# ========= 2) AF data =========
af_dt <- fread(af_file)
stopifnot("gene" %in% names(af_dt))

# ========= 3) 从 AF 里识别 si genes =========
# 默认规则：基因名以 si_ 开头，或包含 |si| 这种 token（你可按实际命名改这里）
si_pattern <- "(^si[_\\.\\-])|(\\bsi\\b)|([_\\.\\-]si[_\\.\\-])|([_\\.\\-]si$)"
si_genes <- sort(unique(af_dt[grepl(si_pattern, gene, ignore.case = TRUE), gene]))

cat("[INFO] si genes detected from AF:", length(si_genes), "\n")
writeLines(si_genes, out_si)
cat("[OK] wrote si gene list to:", out_si, "\n")

# ========= 4) 合并保留 gene 集合并过滤 =========
keep_genes <- sort(unique(c(subunit_genes, si_genes)))
af_keep <- af_dt[gene %in% keep_genes]
fwrite(af_keep, out_af)

# ========= 5) 覆盖统计 + 缺失 =========
genes_in_af <- sort(unique(af_keep$gene))
missing_genes <- setdiff(keep_genes, genes_in_af)

sum_dt <- data.table(
  n_subunit_in_role_table = length(subunit_genes),
  n_si_detected_from_AF   = length(si_genes),
  n_keep_union            = length(keep_genes),
  n_keep_with_AF          = length(genes_in_af),
  n_keep_missing_AF       = length(missing_genes),
  n_rows_output           = nrow(af_keep)
)
fwrite(sum_dt, out_sum)
writeLines(missing_genes, out_miss)

cat("[OK] wrote filtered AF to:", out_af, "\n")
cat("[OK] wrote summary to:", out_sum, "\n")
cat("[OK] wrote missing genes list to:", out_miss, "\n")
cat("\nKeep genes present in AF:", length(genes_in_af), "\n")













72 subunit af files
/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz

q2_parallelism_from_af72.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ========= 输入 =========
af_file  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
cov_file <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
out_dir  <- "/mnt/spareHD_2/nu_287/q2_parallelism"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========= 参数 =========
min_depth <- 10
min_n_fresh_AK <- 5
min_n_fresh_BC <- 8

# ========= 分组 =========
AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")
AK_marine <- "RS"
BC_marine <- "SAY"

# ========= 读数据（注意：AF 是 CSV！） =========
AF  <- fread(cmd = paste("zcat", shQuote(af_file)), sep = ",", header = TRUE)
COV <- fread(cov_file)

# pop 名标准化：10_THE_S10 -> THE
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}
AF[, pop := normalize_pop(pop)]
COV[, pop := toupper(pop)]

# depth过滤
AF <- AF[depth >= min_depth]
AF[, snp := paste(chr, pos, gene, sep=":")]

# 只保留我们关心的 pops
keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]

# helper：给定一个 SNP，算某组 pops 的均值 AF（并记录n）
mean_af <- function(dt, pops){
  sub <- dt[pop %in% pops]
  if(nrow(sub)==0) return(list(mu=NA_real_, n=0L))
  list(mu=mean(sub$af, na.rm=TRUE), n=as.integer(sum(!is.na(sub$af))))
}

# ========= SNP层面的 ΔAF =========
snps <- unique(AF$snp)

res <- vector("list", length(snps))
names(res) <- snps

for(i in seq_along(snps)){
  s <- snps[i]
  dt <- AF[snp == s, .(pop, af)]

  mu_RS  <- mean_af(dt, AK_marine)
  mu_SAY <- mean_af(dt, BC_marine)
  mu_AK  <- mean_af(dt, AK_fresh)
  mu_BC  <- mean_af(dt, BC_fresh)

  if(mu_RS$n < 1 || mu_SAY$n < 1) next
  if(mu_AK$n < min_n_fresh_AK || mu_BC$n < min_n_fresh_BC) next

  dAK <- mu_AK$mu - mu_RS$mu
  dBC <- mu_BC$mu - mu_SAY$mu

  concordant <- sign(dAK) == sign(dBC) && sign(dAK)!=0 && sign(dBC)!=0
  prodv <- dAK * dBC

  parts <- strsplit(s, ":", fixed=TRUE)[[1]]
  chr_ <- parts[1]; pos_ <- as.integer(parts[2]); gene_ <- parts[3]

  res[[i]] <- data.frame(
    chr=chr_, pos=pos_, gene=gene_, snp=s,
    AF_RS=mu_RS$mu,  n_RS=mu_RS$n,
    AF_SAY=mu_SAY$mu, n_SAY=mu_SAY$n,
    AF_AKfresh=mu_AK$mu, n_AKfresh=mu_AK$n,
    AF_BCfresh=mu_BC$mu, n_BCfresh=mu_BC$n,
    dAF_AK=dAK, dAF_BC=dBC,
    prod=prodv,
    concordant=concordant
  )
}

SNP <- bind_rows(res)
SNP <- SNP[is.finite(SNP$dAF_AK) & is.finite(SNP$dAF_BC), ]

# ========= 全局并行性统计 =========
cor_pearson  <- cor(SNP$dAF_AK, SNP$dAF_BC, method="pearson")
cor_spearman <- cor(SNP$dAF_AK, SNP$dAF_BC, method="spearman")

ok <- SNP$dAF_AK!=0 & SNP$dAF_BC!=0
k <- sum(SNP$concordant[ok])
n <- sum(ok)
binom_p <- binom.test(k, n, p=0.5)$p.value

stats <- data.frame(
  n_snps=nrow(SNP),
  n_nonzero=n,
  concordant=k,
  concordant_prop=k/n,
  cor_pearson=cor_pearson,
  cor_spearman=cor_spearman,
  binom_p=binom_p
)

# ========= 输出 =========
fwrite(SNP,   file.path(out_dir, "q2_parallelism_perSNP.tsv.gz"), sep="\t")
fwrite(stats, file.path(out_dir, "q2_parallelism_summary.tsv"), sep="\t")

# 图：ΔAF散点
p1 <- ggplot(SNP, aes(dAF_AK, dAF_BC)) +
  geom_point(alpha=0.35, size=0.8) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_smooth(method="lm", se=FALSE) +
  labs(
    x=expression(Delta~AF[AK]~" (AK freshwater mean - RS)"),
    y=expression(Delta~AF[BC]~" (BC freshwater mean - SAY)"),
    title="Q2 parallelism (72 features): per-SNP ΔAF concordance across AK & BC"
  )
ggsave(file.path(out_dir, "q2_parallelism_scatter_perSNP.png"), p1, width=6.5, height=5, dpi=300)

# 图：prod分布
p2 <- ggplot(SNP, aes(prod)) +
  geom_histogram(bins=80) +
  geom_vline(xintercept=0) +
  labs(x="ΔAF_AK * ΔAF_BC", title="Positive values indicate same-direction shifts")
ggsave(file.path(out_dir, "q2_parallelism_prod_hist.png"), p2, width=6.5, height=4.5, dpi=300)

# gene层面汇总
GENE <- SNP %>%
  group_by(gene) %>%
  summarise(
    n_snps=n(),
    dAK_med=median(dAF_AK, na.rm=TRUE),
    dBC_med=median(dAF_BC, na.rm=TRUE),
    concordant_prop=mean(concordant[dAF_AK!=0 & dAF_BC!=0], na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(prod_med = dAK_med * dBC_med)

fwrite(GENE, file.path(out_dir, "q2_parallelism_perGene.tsv"), sep="\t")

p3 <- ggplot(GENE, aes(dAK_med, dBC_med)) +
  geom_point(alpha=0.7) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_smooth(method="lm", se=FALSE) +
  labs(
    x="gene median ΔAF_AK",
    y="gene median ΔAF_BC",
    title="Q2 parallelism (72 features): gene-level median ΔAF"
  )
ggsave(file.path(out_dir, "q2_parallelism_scatter_perGene.png"), p3, width=6.5, height=5, dpi=300)

cat("[OK] wrote outputs to:", out_dir, "\n")
print(stats)




#within
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ========= 输入 =========
af_file <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"  # CSV gz
out_dir <- "/mnt/spareHD_2/nu_287/q2_parallelism/within_region_strict_rawGene_min10"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========= 参数（主分析推荐版） =========
min_depth <- 10

# 淡水群体“有AF数据”的最小数量
min_n_AK <- 5   # out of 9
min_n_BC <- 8   # out of 13

# 淡水群体中“有明确方向(非零ΔAF)”的最小数量
min_nonzero_AK <- 5
min_nonzero_BC <- 8

# 把小效应当成 0（无方向）
eps <- 0.01  # 可试 0.005 / 0.02

# SNP 层面的“平行”阈值：同向比例
parallel_prop_cutoff <- 0.70

# ✅ gene 层面：strict 版最低 SNP 数（你要的：先用 10）
min_snps_per_gene <- 10

# gene 层面：定义“平行基因”的阈值（在 strict gene 上做）
gene_parallel_snp_fraction_cutoff <- 0.30

# ========= 分组 =========
AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")
AK_marine <- "RS"
BC_marine <- "SAY"

# ========= 读 AF（CSV） =========
AF <- fread(cmd = paste("zcat", shQuote(af_file)), sep = ",", header = TRUE)
stopifnot(all(c("chr","pos","gene","pop","af","depth") %in% names(AF)))

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}
AF[, pop := normalize_pop(pop)]

# depth 过滤
AF <- AF[depth >= min_depth]
AF[, snp := paste(chr, pos, gene, sep=":")]

# 只保留关心的 pops
keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]

# ========= reshape：每个 snp 一行，每个 pop 一列 =========
W <- dcast(AF, chr + pos + gene + snp ~ pop, value.var = "af",
           fun.aggregate = mean, fill = NA_real_)

# ========= helper：地区内部平行（每 SNP） =========
within_region_per_snp <- function(W, fresh_pops, marine_pop, region_name,
                                  min_n_fresh, min_nonzero, eps){

  if(!(marine_pop %in% colnames(W))) stop("missing marine pop column: ", marine_pop)
  marine <- W[[marine_pop]]

  # ΔAF 矩阵：每行 SNP，每列 淡水pop
  deltas <- lapply(fresh_pops, function(p){
    if(!(p %in% colnames(W))) rep(NA_real_, nrow(W)) else W[[p]] - marine
  })
  D <- do.call(cbind, deltas)
  colnames(D) <- fresh_pops

  # 有AF数据的淡水群体数（ΔAF非NA）
  n_fresh <- rowSums(!is.na(D))

  # 把小效应当成 0（无方向）
  D2 <- D
  D2[!is.na(D2) & abs(D2) < eps] <- 0

  # 方向符号（-1,0,1）
  S <- sign(D2)

  # 有明确方向（非零）的淡水群体数
  n_nonzero <- rowSums(!is.na(S) & S != 0)

  # 多数方向：对符号求和再取 sign
  sum_sign <- rowSums(ifelse(is.na(S), 0, S))
  maj_sign <- sign(sum_sign)

  # 同向计数：符号与 maj_sign 相同，且非零
  concord <- rowSums(!is.na(S) & S != 0 & (S == maj_sign))

  # 同向比例（只在非零淡水群体上算）
  concord_prop <- ifelse(n_nonzero > 0, concord / n_nonzero, NA_real_)

  # 连续一致性：非零符号的均值绝对值
  mean_sign <- apply(S, 1, function(x){
    x <- x[!is.na(x) & x != 0]
    if(length(x) == 0) return(NA_real_)
    mean(x)
  })
  sign_consistency <- abs(mean_sign)

  # 效应大小：淡水 ΔAF 的中位数（用原始 D）
  dAF_med <- apply(D, 1, function(x) median(x, na.rm = TRUE))

  out <- data.table(
    chr = W$chr,
    pos = W$pos,
    gene = W$gene,
    snp = W$snp,
    region = region_name,
    n_fresh = n_fresh,
    n_nonzero = n_nonzero,
    maj_sign = maj_sign,
    concord_prop = concord_prop,
    sign_consistency = sign_consistency,
    dAF_med = dAF_med
  )

  # ✅ 关键过滤：既要有足够数据，也要有足够“非零方向信息”
  out <- out[n_fresh >= min_n_fresh & n_nonzero >= min_nonzero]
  out
}

AK_SNP <- within_region_per_snp(W, AK_fresh, AK_marine, "AK",
                                min_n_AK, min_nonzero_AK, eps)
BC_SNP <- within_region_per_snp(W, BC_fresh, BC_marine, "BC",
                                min_n_BC, min_nonzero_BC, eps)

# ========= 输出 SNP 层面 =========
fwrite(AK_SNP, file.path(out_dir, "within_AK_perSNP.strict.tsv.gz"), sep="\t")
fwrite(BC_SNP, file.path(out_dir, "within_BC_perSNP.strict.tsv.gz"), sep="\t")

# ========= gene 汇总函数：不过滤 n_snps（raw） =========
gene_summarise_raw <- function(SNP_dt, cutoff){
  SNP_dt %>%
    group_by(gene) %>%
    summarise(
      n_snps = n(),
      concord_prop_med = median(concord_prop, na.rm = TRUE),
      sign_consistency_med = median(sign_consistency, na.rm = TRUE),
      dAF_med_med = median(dAF_med, na.rm = TRUE),
      prop_snps_parallel = mean(concord_prop >= cutoff, na.rm = TRUE),
      prop_snps_consistent = mean(sign_consistency >= 0.8, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.table()
}

# ========= 1) Raw gene 表：不过滤 n_snps =========
AK_GENE_raw <- gene_summarise_raw(AK_SNP, parallel_prop_cutoff)
BC_GENE_raw <- gene_summarise_raw(BC_SNP, parallel_prop_cutoff)

fwrite(AK_GENE_raw, file.path(out_dir, "within_AK_perGene.raw.tsv"), sep="\t")
fwrite(BC_GENE_raw, file.path(out_dir, "within_BC_perGene.raw.tsv"), sep="\t")

# ========= 2) Strict gene 表：过滤 n_snps（min_snps_per_gene=10） =========
AK_GENE <- AK_GENE_raw[n_snps >= min_snps_per_gene]
BC_GENE <- BC_GENE_raw[n_snps >= min_snps_per_gene]

fwrite(AK_GENE, file.path(out_dir, "within_AK_perGene.strict.tsv"), sep="\t")
fwrite(BC_GENE, file.path(out_dir, "within_BC_perGene.strict.tsv"), sep="\t")

# ========= 定义“平行基因”集合（strict 上做） =========
AK_parallel_genes <- AK_GENE[prop_snps_parallel >= gene_parallel_snp_fraction_cutoff, gene]
BC_parallel_genes <- BC_GENE[prop_snps_parallel >= gene_parallel_snp_fraction_cutoff, gene]
INTERSECT <- sort(intersect(AK_parallel_genes, BC_parallel_genes))

writeLines(AK_parallel_genes, file.path(out_dir, "AK_parallel_genes.strict.txt"))
writeLines(BC_parallel_genes, file.path(out_dir, "BC_parallel_genes.strict.txt"))
writeLines(INTERSECT,        file.path(out_dir, "AK_BC_parallel_genes_intersect.strict.txt"))

# ========= 合并 gene 表（raw + strict） =========
AK_GENE_raw2 <- copy(AK_GENE_raw); AK_GENE_raw2[, region := "AK"]
BC_GENE_raw2 <- copy(BC_GENE_raw); BC_GENE_raw2[, region := "BC"]
GENE_ALL_raw <- rbindlist(list(AK_GENE_raw2, BC_GENE_raw2), use.names=TRUE, fill=TRUE)
fwrite(GENE_ALL_raw, file.path(out_dir, "within_regions_perGene_merged.raw.tsv"), sep="\t")

AK_GENE2 <- copy(AK_GENE); AK_GENE2[, region := "AK"]
BC_GENE2 <- copy(BC_GENE); BC_GENE2[, region := "BC"]
GENE_ALL <- rbindlist(list(AK_GENE2, BC_GENE2), use.names=TRUE, fill=TRUE)
fwrite(GENE_ALL, file.path(out_dir, "within_regions_perGene_merged.strict.tsv"), sep="\t")

# ========= 图1：AK vs BC gene-level（raw，永远有） =========
GENE_WIDE_raw <- dcast(GENE_ALL_raw, gene ~ region, value.var = "prop_snps_parallel")
if(all(c("AK","BC") %in% colnames(GENE_WIDE_raw))) {
  p1 <- ggplot(GENE_WIDE_raw, aes(AK, BC)) +
    geom_point(alpha=0.8) +
    geom_hline(yintercept = gene_parallel_snp_fraction_cutoff) +
    geom_vline(xintercept = gene_parallel_snp_fraction_cutoff) +
    labs(
      x="AK: fraction of parallel SNPs per gene (raw)",
      y="BC: fraction of parallel SNPs per gene (raw)",
      title=paste0("Within-region parallelism (raw genes): eps=", eps,
                   ", SNP cutoff=", parallel_prop_cutoff)
    )
  ggsave(file.path(out_dir, "gene_scatter.propParallel.raw.png"),
         p1, width=6.5, height=5, dpi=300)
}

# ========= 图2：AK vs BC gene-level（strict，可选，自动跳过缺列） =========
GENE_WIDE <- dcast(GENE_ALL, gene ~ region, value.var = "prop_snps_parallel")
if(all(c("AK","BC") %in% colnames(GENE_WIDE))) {
  p2 <- ggplot(GENE_WIDE, aes(AK, BC)) +
    geom_point(alpha=0.8) +
    geom_hline(yintercept = gene_parallel_snp_fraction_cutoff) +
    geom_vline(xintercept = gene_parallel_snp_fraction_cutoff) +
    labs(
      x="AK: fraction of parallel SNPs per gene (strict)",
      y="BC: fraction of parallel SNPs per gene (strict)",
      title=paste0("Within-region parallelism (strict genes): min_snps=", min_snps_per_gene,
                   ", eps=", eps, ", SNP cutoff=", parallel_prop_cutoff)
    )
  ggsave(file.path(out_dir, "gene_scatter.propParallel.strict.png"),
         p2, width=6.5, height=5, dpi=300)
} else {
  cat("[WARN] Skip strict gene scatter plot: missing AK/BC columns.\n")
  cat("[WARN] strict gene rows: AK=", nrow(AK_GENE), " BC=", nrow(BC_GENE), "\n", sep="")
  cat("[WARN] Try lowering min_snps_per_gene further (e.g., 5) if AK still 0.\n")
}

# ========= 图3：SNP 层面 concord_prop 分布 =========
SNP_ALL <- rbindlist(list(AK_SNP, BC_SNP), use.names=TRUE, fill=TRUE)
p3 <- ggplot(SNP_ALL, aes(concord_prop)) +
  geom_histogram(bins=50) +
  facet_wrap(~region, ncol=1, scales="free_y") +
  labs(x="Within-region concordant proportion (per SNP)",
       title="Distribution of within-region concordance (strict SNP filtering)")
ggsave(file.path(out_dir, "SNP_concordProp_hist.strict.png"),
       p3, width=6.5, height=6.0, dpi=300)

# ========= summary =========
cat("[OK] wrote outputs to:", out_dir, "\n")
cat("[PARAM] eps=", eps, " min_depth=", min_depth,
    " min_n_AK=", min_n_AK, " min_nonzero_AK=", min_nonzero_AK,
    " min_n_BC=", min_n_BC, " min_nonzero_BC=", min_nonzero_BC,
    " SNP_cutoff=", parallel_prop_cutoff,
    " min_snps_per_gene(strict)=", min_snps_per_gene,
    " gene_frac_cutoff(strict)=", gene_parallel_snp_fraction_cutoff,
    "\n", sep="")

cat("[AK] SNPs kept:", nrow(AK_SNP),
    " raw genes:", nrow(AK_GENE_raw),
    " strict genes:", nrow(AK_GENE),
    " strict parallel genes:", length(AK_parallel_genes), "\n")

cat("[BC] SNPs kept:", nrow(BC_SNP),
    " raw genes:", nrow(BC_GENE_raw),
    " strict genes:", nrow(BC_GENE),
    " strict parallel genes:", length(BC_parallel_genes), "\n")

cat("[AK ∩ BC] strict parallel genes:", length(INTERSECT), "\n")









#主图
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ========= 输入 =========
in_gene <- "/mnt/spareHD_2/nu_287/q2_parallelism/within_region_strict_rawGene_min10/within_regions_perGene_merged.raw.tsv"
out_dir <- "/mnt/spareHD_2/nu_287/q2_parallelism/within_region_strict_rawGene_min10"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 阈值：你提的 0.3 / 0.5
thr1 <- 0.30
thr2 <- 0.50

# ========= 读 gene 表（long -> wide）=========
G <- fread(in_gene)
stopifnot(all(c("gene","region","prop_snps_parallel") %in% names(G)))

W <- dcast(G, gene ~ region, value.var = "prop_snps_parallel")
# 只保留同时有 AK & BC 的
W <- W[!is.na(AK) & !is.na(BC)]

# ========= 打标：两地区都超过阈值 =========
W[, pass_03 := (AK >= thr1 & BC >= thr1)]
W[, pass_05 := (AK >= thr2 & BC >= thr2)]

# 输出基因列表（回答“哪些 nuclear subunits 在 AK 和 BC 都 parallel”）
writeLines(W[pass_03 == TRUE, gene], file.path(out_dir, "genes_parallel_AKBC_both_ge0.30.txt"))
writeLines(W[pass_05 == TRUE, gene], file.path(out_dir, "genes_parallel_AKBC_both_ge0.50.txt"))

# ========= 作图 =========
p <- ggplot(W, aes(x = AK, y = BC)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_point(data = W[pass_03 == TRUE], alpha = 0.95, size = 2.5) +
  geom_point(data = W[pass_05 == TRUE], alpha = 1.00, size = 3) +
  geom_vline(xintercept = thr1) + geom_hline(yintercept = thr1) +
  geom_vline(xintercept = thr2, linetype = "dashed") +
  geom_hline(yintercept = thr2, linetype = "dashed") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(
    x = "AK: fraction of parallel SNPs per gene",
    y = "BC: fraction of parallel SNPs per gene",
    title = "Within-region parallelism (gene-level): AK vs BC",
    subtitle = "Solid = 0.3, dashed = 0.5; highlighted points pass both thresholds in both regions"
  )

ggsave(file.path(out_dir, "Fig_Q2_gene_parallel_scatter_AK_vs_BC.png"),
       p, width = 6.5, height = 5.2, dpi = 300)

cat("[OK] wrote:\n")
cat(" - ", file.path(out_dir, "Fig_Q2_gene_parallel_scatter_AK_vs_BC.png"), "\n", sep="")
cat(" - ", file.path(out_dir, "genes_parallel_AKBC_both_ge0.30.txt"), "\n", sep="")
cat(" - ", file.path(out_dir, "genes_parallel_AKBC_both_ge0.50.txt"), "\n", sep="")



#
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ====== 输入输出 ======
in_file <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_parallelism_perSNP.tsv.gz"
out_png <- "/mnt/spareHD_2/nu_287/q2_parallelism/Fig_Q1_crossRegion_perSNP_dAF_trend.png"
out_pdf <- "/mnt/spareHD_2/nu_287/q2_parallelism/Fig_Q1_crossRegion_perSNP_dAF_trend.pdf"

# ====== 读数据 ======
SNP <- fread(in_file)
stopifnot(all(c("dAF_AK","dAF_BC") %in% names(SNP)))

# ====== 相关性（用于图注/正文）=====
pearson <- cor.test(SNP$dAF_AK, SNP$dAF_BC, method="pearson")
spearman <- cor.test(SNP$dAF_AK, SNP$dAF_BC, method="spearman")  # ties 会提示 warning，正常

cat("\n[Q1 trend] n_snps =", nrow(SNP), "\n")
cat("[Q1 trend] Pearson r =", unname(pearson$estimate), "p =", pearson$p.value, "\n")
cat("[Q1 trend] Spearman rho =", unname(spearman$estimate), "p =", spearman$p.value, "\n\n")

# ====== 作图 ======
p <- ggplot(SNP, aes(dAF_AK, dAF_BC)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method="lm", se=FALSE) +
  labs(
    x = expression(Delta~AF[AK]~"(AK freshwater mean - RS)"),
    y = expression(Delta~AF[BC]~"(BC freshwater mean - SAY)"),
    title = "Q1 trend (72 features): per-SNP ΔAF concordance across AK & BC"
  ) +
  theme_classic(base_size = 14)

ggsave(out_png, p, width = 7, height = 5.8, dpi = 300)
ggsave(out_pdf, p, width = 7, height = 5.8)

cat("[OK] wrote:\n", out_png, "\n", out_pdf, "\n")



#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ====== 输入目录 ======
dir_in <- "/mnt/spareHD_2/nu_287/q2_parallelism/within_region_strict"  # ✅ 你现在生成 174 的目录
ak_file <- file.path(dir_in, "within_AK_perSNP.strict.tsv.gz")
bc_file <- file.path(dir_in, "within_BC_perSNP.strict.tsv.gz")
par_file <- file.path(dir_in, "AK_BC_parallel_snps.txt")              # ✅ 174 SNP 清单

# ====== 输出 ======
out_png <- file.path(dir_in, "Fig_Q2_withinRegion_parallel_174_gold15.png")
out_pdf <- file.path(dir_in, "Fig_Q2_withinRegion_parallel_174_gold15.pdf")
out_15_tsv <- file.path(dir_in, "Table_Q2_gold15_highEffect_highConsistency_snps.tsv")

# ====== 参数（与你定义 gold15 一致）=====
cut_dAF <- 0.5
cut_conc <- 0.7

# ====== 读 within SNP 表（AK/BC）=====
AK <- fread(cmd = paste("zcat", shQuote(ak_file)))
BC <- fread(cmd = paste("zcat", shQuote(bc_file)))

# 你的 within 表里列名是 dAF_med / concord_prop（从你贴的输出推断）
stopifnot(all(c("snp","gene","dAF_med","concord_prop") %in% names(AK)))
stopifnot(all(c("snp","gene","dAF_med","concord_prop") %in% names(BC)))

setnames(AK, c("dAF_med","concord_prop"), c("dAF_med_AK","concord_AK"))
setnames(BC, c("dAF_med","concord_prop"), c("dAF_med_BC","concord_BC"))

# 合并到同一张表：只保留 AK 和 BC 都存在的 SNP（你看到的 M rows ~291）
M <- merge(
  AK[, .(snp, chr, pos, gene, dAF_med_AK, concord_AK)],
  BC[, .(snp, dAF_med_BC, concord_BC)],
  by = "snp"
)

cat("Total SNPs after within filtering:", nrow(M), "\n")

# ====== 读 174 SNP 列表 ======
parallel_snps <- fread(par_file, header = FALSE)[[1]]
M[, is_parallel174 := snp %in% parallel_snps]
cat("AK ∩ BC parallel SNPs:", sum(M$is_parallel174), "\n")  # 应该 = 174

# ====== 定义 gold15（在 174 里再筛）=====
M[, is_gold15 :=
    is_parallel174 &
    dAF_med_AK >= cut_dAF & dAF_med_BC >= cut_dAF &
    concord_AK >= cut_conc & concord_BC >= cut_conc
]
cat("High-effect & high-consistency SNPs:", sum(M$is_gold15), "\n")   # 你这里是 15

# ====== 把 15 个 SNP 存表 ======
SNP_15 <- M[is_gold15 == TRUE,
            .(snp, chr, pos, gene, dAF_med_AK, dAF_med_BC, concord_AK, concord_BC)]
fwrite(SNP_15[order(gene, -dAF_med_AK)], out_15_tsv, sep="\t")
cat("[OK] wrote gold15 table:\n", out_15_tsv, "\n")

# ====== 作图（灰底 + 红174 + 金15 + 0线 + 0.5线）=====
p <- ggplot(M, aes(dAF_med_AK, dAF_med_BC)) +
  # 背景：所有 SNP
  geom_point(alpha = 0.25, size = 2) +
  # 红色：174
  geom_point(data = M[is_parallel174 == TRUE], alpha = 0.85, size = 2.3) +
  # 金色：15
  geom_point(data = M[is_gold15 == TRUE], alpha = 0.95, size = 2.8) +
  # 0 线
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # 0.5 阈值线
  geom_hline(yintercept = cut_dAF, linetype = "dotted") +
  geom_vline(xintercept = cut_dAF, linetype = "dotted") +
  labs(
    x = "AK within-region median ΔAF",
    y = "BC within-region median ΔAF",
    title = "Within-region parallel allele-frequency shifts in AK vs BC",
    subtitle = "Grey: all SNPs | Red: 174 AK∩BC parallel SNPs | Gold: high-effect & high-consistency SNPs"
  ) +
  theme_classic(base_size = 16)

ggsave(out_png, p, width = 7, height = 5.8, dpi = 300)
ggsave(out_pdf, p, width = 7, height = 5.8)
cat("[OK] wrote:\n", out_png, "\n", out_pdf, "\n")

