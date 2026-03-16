#removing amo linear model
#去掉amo
Rscript - <<'RS'
library(ape)

# mt tree（和你画 noAMO 图用的是同一个）
tr <- read.tree("/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile")

# drop AMO
tr <- drop.tip(tr, "AMO")

# 距离矩阵
D <- cophenetic(tr)

# PCoA（一定要 eig=TRUE 才能算解释方差）
pco <- cmdscale(D, k = 5, eig = TRUE)

# 解释方差百分比
explained <- 100 * pco$eig / sum(pco$eig)

cat("Explained variance (no AMO):\n")
for(i in 1:5){
  cat(sprintf("  PC%d: %.2f%%\n", i, explained[i]))
}
RS
PC1: 49.54% PC2: 36.83% PC3: 10.06% PC4: 4.36% PC5: 2.81%





# 输入（你现在的 72 subunit AF long）
AF72=/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz

# covariates 源（含 treePC）
COV=/mnt/spareHD_2/nu_287/covariates.treePC.tsv

# mt tree
MT_TREE=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile

# 输出目录
OUTD=/mnt/spareHD_2/nu_287/_assoc72_subunit_AMOdrop
mkdir -p "$OUTD"

Rscript - <<'RS'
suppressPackageStartupMessages({library(data.table); library(dplyr); library(ape)})

cov_in  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
mt_tree <- "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile"
cov_out <- "/mnt/spareHD_2/nu_287/_assoc72_subunit_AMOdrop/covariates.treePC.noAMO.tsv"

COV <- fread(cov_in)

# ---- 1) drop AMO from mt tree then recompute mitoPC ----
tr0 <- read.tree(mt_tree)
if("AMO" %in% tr0$tip.label){
  tr <- drop.tip(tr0, "AMO")
} else {
  tr <- tr0
}

D <- cophenetic(tr)
common <- intersect(COV$pop, rownames(D))
stopifnot(length(common) >= 3)

D2 <- D[common, common, drop=FALSE]
pco <- cmdscale(D2, k=3)
mitoPC <- data.frame(pop=rownames(pco),
                     mitoPC1=pco[,1], mitoPC2=pco[,2], mitoPC3=pco[,3])

# ---- 2) merge back with treePC covariates, and remove AMO row ----
COV2 <- COV %>%
  select(-starts_with("mitoPC")) %>%          # 换成新 mitoPC
  inner_join(mitoPC, by="pop") %>%
  filter(pop != "AMO")

fwrite(COV2, cov_out, sep="\t")
cat("[ok] wrote:", cov_out, " n_pops=", nrow(COV2), "\n")
RS

Rscript - <<'RS'
suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# ======================
# inputs / outputs
# ======================
mt_tree <- "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile"
out_dir <- "/mnt/spareHD_2/nu_287/_assoc72_subunit"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ======================
# helper: region labels
# ======================
add_region <- function(df){
  df %>%
    mutate(region = case_when(
      pop %in% c("FG","LG","SR","SL","TL","WB","WT","WK","LB") ~ "Alaska",
      pop %in% c("SWA","THE","JOE","BEA","MUC","PYE","ROS","BOOT","ECHO","LAW","GOS","ROB") ~ "BC",
      pop %in% c("RS","SAY") ~ "Marine",
      pop %in% c("PACH","FRED","SC","CH") ~ "Recent",
      TRUE ~ "Other"
    )) %>%
    arrange(region, pop)
}

# ======================
# theme: border + grid
# ======================
theme_border_grid <- function(){
  theme_bw(base_size = 12) +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.9),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),  # 如果你想要次网格：改成 element_line(color="grey92", linewidth=0.3)
      panel.background = element_rect(fill = "white"),
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      legend.title     = element_blank(),
      legend.position  = "right"
    )
}

# ======================
# 1) read tree + drop AMO
# ======================
tr <- read.tree(mt_tree)
if(!("AMO" %in% tr$tip.label)){
  stop("AMO not found in mt tree tips. Check tip labels.")
}
tr2 <- drop.tip(tr, "AMO")

# ======================
# 2) PCoA + explained variance
# ======================
D   <- cophenetic(tr2)
pco <- cmdscale(D, k = 3, eig = TRUE)

expl <- 100 * pco$eig / sum(pco$eig)
cat(sprintf("Explained variance (no AMO): PC1=%.2f%% PC2=%.2f%% PC3=%.2f%%\n",
            expl[1], expl[2], expl[3]))

PC <- data.frame(
  pop    = rownames(pco$points),
  mitoPC1 = pco$points[,1],
  mitoPC2 = pco$points[,2],
  mitoPC3 = pco$points[,3]
) %>% add_region()

# ======================
# 3) plots
# ======================
p_base <- ggplot(PC, aes(mitoPC1, mitoPC2, color = region)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    x = sprintf("mitoPC1 (%.1f%%)", expl[1]),
    y = sprintf("mitoPC2 (%.1f%%)", expl[2]),
    title = "Mitochondrial PCoA after removing AMO"
  ) +
  theme_border_grid()

# with labels
p_label <- p_base +
  ggrepel::geom_text_repel(
    aes(label = pop),
    size = 3,
    max.overlaps = Inf,
    seed = 1
  )

# without labels (clean)
p_nolabel <- p_base

# ======================
# 4) save
# ======================
fn1_png <- file.path(out_dir, "mitoPC_noAMO_PC1_PC2.bordergrid.png")
fn1_pdf <- file.path(out_dir, "mitoPC_noAMO_PC1_PC2.bordergrid.pdf")
fn2_png <- file.path(out_dir, "mitoPC_noAMO_PC1_PC2.bordergrid.nolabel.png")
fn2_pdf <- file.path(out_dir, "mitoPC_noAMO_PC1_PC2.bordergrid.nolabel.pdf")

ggsave(fn1_png, p_label,   width = 7.2, height = 5.2, dpi = 300)
ggsave(fn1_pdf, p_label,   width = 7.2, height = 5.2)
ggsave(fn2_png, p_nolabel, width = 7.2, height = 5.2, dpi = 300)
ggsave(fn2_pdf, p_nolabel, width = 7.2, height = 5.2)

cat("[ok] wrote:\n", fn1_png, "\n", fn1_pdf, "\n", fn2_png, "\n", fn2_pdf, "\n")
RS




  Rscript /mnt/spareHD_2/nu_287/run_assoc_perSNP_haplo.R \
  /mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz \
  /mnt/spareHD_2/nu_287/covariates.treePC.noAMO_mitoPC.tsv \
  /mnt/spareHD_2/nu_287/_assoc72_subunit/assoc72_subunit_noAMO







#noamo
# pc percentage  explained    /mnt/spareHD_2/nu_287/_assoc72_subunit/01_plot_mitoPC_noAMO.R

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

# drop AMO
tr <- read.tree(mt_tree)
if("AMO" %in% tr$tip.label) tr <- drop.tip(tr, "AMO")

D <- cophenetic(tr)
pco <- cmdscale(D, k=5, eig=TRUE)

expl <- 100 * pco$eig / sum(pco$eig)
cat(sprintf("[mitoPC explained noAMO] PC1=%.2f%% PC2=%.2f%% PC3=%.2f%% PC4=%.2f%% PC5=%.2f%%\n",
            expl[1], expl[2], expl[3], expl[4], expl[5]))

PC <- data.frame(pop=rownames(pco$points),
                 mitoPC1=pco$points[,1],
                 mitoPC2=pco$points[,2],
                 mitoPC3=pco$points[,3])

PC <- PC %>%
  mutate(region = case_when(
    pop %in% c("FG","LG","SR","SL","TL","WB","WT","WK","LB") ~ "Alaska",
    pop %in% c("SWA","THE","JOE","BEA","MUC","PYE","ROS","BOOT","ECHO","LAW","GOS","ROB") ~ "BC",
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
    title="Mitochondrial PCoA after removing AMO"
  ) +
  theme_bw(base_size=12) +        # ✅ 有边框
  theme(
    panel.grid.major = element_line(linewidth=0.3),
    panel.grid.minor = element_line(linewidth=0.2),
    legend.position="right"
  )

png(file.path(out_dir, "mitoPC_noAMO_PC1_PC2.theme_bw.png"), width=2200, height=1600, res=300)
print(p12)
dev.off()

pdf(file.path(out_dir, "mitoPC_noAMO_PC1_PC2.theme_bw.pdf"), width=7.2, height=5.2)
print(p12)
dev.off()

# 写 explained variance 表，补充材料用
exp_out <- file.path(out_dir, "mitoPC_noAMO_explained.tsv")
fwrite(data.table(PC=paste0("PC",1:5), explained_percent=expl[1:5]), exp_out, sep="\t")
cat("[write] ", exp_out, "\n")




#/mnt/spareHD_2/nu_287/_assoc72_subunit/02_make_covariates_noAMO_treePC2.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ape)
})

# 输入
cov_in  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"  # 你原本 27pop 的 covariates（含 mitoPC、treePC）
mt_tree <- "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/aligned_mt_noDloop.fasta.treefile"

out_dir <- "/mnt/spareHD_2/nu_287/_assoc72_subunit"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
out_cov <- file.path(out_dir, "covariates.noAMO.treePC2_mitoPC2.tsv")

COV <- fread(cov_in)
COV <- COV[pop != "AMO"]
# 重新计算 mitoPC（drop AMO 后重算），避免“轴被旋转但 covariates 没更新”的问题
tr <- read.tree(mt_tree)
if("AMO" %in% tr$tip.label) tr <- drop.tip(tr, "AMO")
D <- cophenetic(tr)
pco <- cmdscale(D, k=3, eig=TRUE)
expl <- 100*pco$eig/sum(pco$eig)

mito <- data.table(pop=rownames(pco$points),
                   mitoPC1=pco$points[,1],
                   mitoPC2=pco$points[,2],
                   mitoPC3=pco$points[,3])

# 只替换 mitoPC*，treePC 仍用你 NJ tree 的 treePC1/2（与你主分析一致）
keep_cols <- setdiff(names(COV), c("mitoPC1","mitoPC2","mitoPC3"))
COV2 <- merge(COV[, ..keep_cols], mito, by="pop", all=FALSE)

fwrite(COV2, out_cov, sep="\t")
cat(sprintf("[write] %s (n=%d)\n", out_cov, nrow(COV2)))
cat(sprintf("[mitoPC explained noAMO] PC1=%.2f%% PC2=%.2f%% PC3=%.2f%%\n", expl[1], expl[2], expl[3]))

# 同时写解释率（给 Supplement）
exp_out <- file.path(out_dir, "covariates.noAMO.mitoPC_explained.tsv")
fwrite(data.table(PC=paste0("PC",1:5), explained_percent=expl[1:5]), exp_out, sep="\t")
cat("[write] ", exp_out, "\n")



Rscript /mnt/spareHD_2/nu_287/run_assoc_perSNP_haplo.R \
  /mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz \
  /mnt/spareHD_2/nu_287/_assoc72_subunit/covariates.noAMO.treePC2_mitoPC2.tsv \
  /mnt/spareHD_2/nu_287/_assoc72_subunit/assoc72_subunit_noAMO


Rscript /mnt/spareHD_2/nu_287/run_assoc_perSNP_haplo_resid.R \
  /mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz \
  /mnt/spareHD_2/nu_287/_assoc72_subunit/covariates.noAMO.treePC2_mitoPC2.tsv \
  /mnt/spareHD_2/nu_287/_assoc72_subunit/assoc72_subunit_noAMO
产出：

assoc72_subunit_noAMO_perSNP.tsv.gz

assoc72_subunit_noAMO_perSNP_resid.tsv.gz

assoc72_subunit_noAMO_covariates_resid.tsv

Rscript /mnt/spareHD_2/nu_287/_assoc72_subunit/check_lambda_fromPerSNP.R \
  /mnt/spareHD_2/nu_287/_assoc72_subunit/assoc72_subunit_noAMO_perSNP.tsv.gz \
  noAMO_raw

Rscript /mnt/spareHD_2/nu_287/_assoc72_subunit/check_lambda_fromPerSNP.R \
  /mnt/spareHD_2/nu_287/_assoc72_subunit/assoc72_subunit_noAMO_perSNP_resid.tsv.gz \
  noAMO_resid


supply 1
/mnt/spareHD_2/nu_287/_assoc72_subunit/QQ_summary_4models/
├── assoc72_QQ_summary_table.tsv        ✅（核心汇总表）
├── QQ_withAMO_raw_PC1.png              ✅
├── QQ_withAMO_resid_PC1.png            ✅
├── QQ_noAMO_raw_PC1.png                ✅
├── QQ_noAMO_resid_PC1.png              ✅
├── QQ_withAMO_raw_PC2.png              ✅
├── QQ_withAMO_resid_PC2.png             ✅
├── QQ_noAMO_raw_PC2.png                 ✅
├── QQ_noAMO_resid_PC2.png               ✅

supply2
/mnt/spareHD_2/nu_287/_assoc72_subunit/
mitoPC_noAMO_PC1_PC2.theme_bw.png        ✅
mitoPC_noAMO_PC1_PC2.theme_bw.pdf        ✅
mitoPC_noAMO_explained.tsv               ✅


/mnt/spareHD_2/nu_287/_assoc72_subunit/gene_stability_4models$ ls
gene_table_all_models.tsv  jaccard_top10_propNominal.tsv  spearman_propNominal_pairwise.tsv  top10_by_propNominal.tsv











#noamo tree rerun
iqtree2 \
  -s /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_noDloop_noAMO.aln.fa \
  -pre /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_noDloop_noAMO.iqtree \
  -m MFP \
  -bb 1000 \
  -nt AUTO \
  -safe


#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

cov_in <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"
pc_in  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"
out_cov <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/covariates.noAMO.rebuiltTree_mitoPC.tsv"

COV <- fread(cov_in)
PC  <- fread(pc_in)

# drop AMO
COV <- COV[pop != "AMO"]

# 替换 mitoPC1-3（用 rebuilt tree 版本）
keep_cols <- setdiff(names(COV), c("mitoPC1","mitoPC2","mitoPC3"))
COV2 <- merge(COV[, ..keep_cols], PC[, .(pop, mitoPC1, mitoPC2, mitoPC3)], by="pop", all=FALSE)

fwrite(COV2, out_cov, sep="\t")
cat(sprintf("[write] %s (n=%d)\n", out_cov, nrow(COV2)))



PC1=49.32% PC2=36.93% PC3=10.20% PC4=4.40% PC5=2.82%


#mtpc1-5 ～ en
Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))

pc_file  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"
env_file <- "/mnt/spareHD_2/nu_287/meta_pop_env.csv"
out_tsv  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mtPC_noAMO_vs_Habitat_LM.tsv"

PC  <- fread(pc_file)
ENV <- fread(env_file)

# ---- find pop column in PC ----
pc_cn_low <- tolower(names(PC))
pc_pop_idx <- which(pc_cn_low %in% c("pop","population","popid","pop_id","site","id"))
if(length(pc_pop_idx)==0) pc_pop_idx <- grep("^pop$|population|pop_id|popid|^id$|site", pc_cn_low)
if(length(pc_pop_idx)==0){
  cat("[ERROR] Cannot find pop column in PC table.\nColumns:\n")
  cat(paste(names(PC), collapse=", "), "\n"); quit("no", status=1)
}
pc_pop_col <- names(PC)[pc_pop_idx[1]]
setnames(PC, pc_pop_col, "pop")

# ---- find pop column in ENV ----
env_cn_low <- tolower(names(ENV))
env_pop_idx <- which(env_cn_low %in% c("pop","population","popid","pop_id","site","id"))
if(length(env_pop_idx)==0) env_pop_idx <- grep("^pop$|population|pop_id|popid|^id$|site", env_cn_low)
if(length(env_pop_idx)==0){
  cat("[ERROR] Cannot find pop column in ENV table.\nColumns:\n")
  cat(paste(names(ENV), collapse=", "), "\n"); quit("no", status=1)
}
env_pop_col <- names(ENV)[env_pop_idx[1]]

# ---- find habitat column in ENV ----
hab_idx <- which(env_cn_low %in% c("habitat","hab","habitat_group","habitatgroup"))
if(length(hab_idx)==0) hab_idx <- grep("habitat", env_cn_low)
if(length(hab_idx)==0){
  cat("[ERROR] Cannot find Habitat column in ENV table.\nColumns:\n")
  cat(paste(names(ENV), collapse=", "), "\n"); quit("no", status=1)
}
hab_col <- names(ENV)[hab_idx[1]]

cat("[info] PC pop col:", pc_pop_col, " -> renamed to pop\n")
cat("[info] ENV pop col:", env_pop_col, "\n")
cat("[info] ENV habitat col:", hab_col, " -> renamed to Habitat\n")

# ---- build df ----
ENV2 <- ENV[, .(pop = get(env_pop_col), Habitat = get(hab_col))]
df <- merge(PC, ENV2, by="pop", all=FALSE)

cat("[info] n populations:", nrow(df), "\n")
cat("[info] Habitat counts:\n")
print(table(df$Habitat))

# baseline Marine if present
df[, Habitat := as.character(Habitat)]
if("Marine" %in% df$Habitat){
  df[, Habitat := factor(Habitat, levels=c("Marine", setdiff(sort(unique(Habitat)), "Marine")))]
} else df[, Habitat := factor(Habitat)]

# PCs present
pc_cols <- grep("^mitoPC[0-9]+$", names(df), value=TRUE)
cat("[info] PCs:", paste(pc_cols, collapse=", "), "\n")

# run LM
res <- rbindlist(lapply(pc_cols, function(pc){
  fit <- lm(as.formula(paste(pc, "~ Habitat")), data=df)
  sm  <- summary(fit)

  co <- as.data.table(coef(sm))
  co[, term := rownames(coef(sm))]
  setnames(co, c("Estimate","Std. Error","t value","Pr(>|t|)"),
              c("beta","se","t","p"))

  out <- co[term!="(Intercept)"]
  out[, `:=`(PC=pc, R2=sm$r.squared, adjR2=sm$adj.r.squared, n=nobs(fit))]
  out[]
}), fill=TRUE)

setcolorder(res, c("PC","term","beta","se","t","p","R2","adjR2","n"))
fwrite(res, out_tsv, sep="\t")
cat("[write]", out_tsv, "\n")

cat("\n=== smallest p per PC ===\n")
print(res[, .(min_p=min(p, na.rm=TRUE)), by=PC][order(min_p)])
RS




#region
Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))

pc_file  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"
env_file <- "/mnt/spareHD_2/nu_287/meta_pop_env.csv"
out_tsv  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mtPC_noAMO_vs_Region_LM.tsv"

PC  <- fread(pc_file)
ENV <- fread(env_file)

# ---- standardize pop column in PC ----
pc_cn_low <- tolower(names(PC))
pc_pop_idx <- grep("^pop$|population|pop_id|popid|^id$|site", pc_cn_low)
stopifnot(length(pc_pop_idx) > 0)
setnames(PC, names(PC)[pc_pop_idx[1]], "pop")

# ---- find pop column in ENV ----
env_cn_low <- tolower(names(ENV))
env_pop_idx <- grep("^pop$|population|pop_id|popid|^id$|site", env_cn_low)
if(length(env_pop_idx)==0){
  cat("[ERROR] Cannot find pop column in ENV.\nColumns:\n")
  cat(paste(names(ENV), collapse=", "), "\n"); quit("no", status=1)
}
env_pop_col <- names(ENV)[env_pop_idx[1]]

# ---- find Region column in ENV ----
reg_idx <- which(env_cn_low %in% c("region","reg","geo_region","georegion"))
if(length(reg_idx)==0) reg_idx <- grep("region", env_cn_low)
if(length(reg_idx)==0){
  cat("[ERROR] Cannot find Region column in ENV.\nColumns:\n")
  cat(paste(names(ENV), collapse=", "), "\n"); quit("no", status=1)
}
reg_col <- names(ENV)[reg_idx[1]]

cat("[info] ENV pop col:", env_pop_col, "\n")
cat("[info] ENV region col:", reg_col, " -> renamed to Region\n")

ENV2 <- ENV[, .(pop = get(env_pop_col), Region = get(reg_col))]
df <- merge(PC, ENV2, by="pop", all=FALSE)

cat("[info] n populations:", nrow(df), "\n")
cat("[info] Region counts:\n"); print(table(df$Region))

df[, Region := factor(as.character(Region))]

pc_cols <- grep("^mitoPC[0-9]+$", names(df), value=TRUE)
cat("[info] PCs:", paste(pc_cols, collapse=", "), "\n")

res <- rbindlist(lapply(pc_cols, function(pc){
  fit <- lm(as.formula(paste(pc, "~ Region")), data=df)
  sm  <- summary(fit)
  co <- as.data.table(coef(sm))
  co[, term := rownames(coef(sm))]
  setnames(co, c("Estimate","Std. Error","t value","Pr(>|t|)"),
              c("beta","se","t","p"))
  out <- co[term!="(Intercept)"]
  out[, `:=`(PC=pc, R2=sm$r.squared, adjR2=sm$adj.r.squared, n=nobs(fit))]
  out[]
}), fill=TRUE)

setcolorder(res, c("PC","term","beta","se","t","p","R2","adjR2","n"))
fwrite(res, out_tsv, sep="\t")
cat("[write]", out_tsv, "\n")

cat("\n=== smallest p per PC ===\n")
print(res[, .(min_p=min(p, na.rm=TRUE)), by=PC][order(min_p)])
RS


#long la
Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))

PC  <- fread("/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv")
ENV <- fread("/mnt/spareHD_2/nu_287/meta_pop_env.csv")

# ---- standardize pop ----
setnames(PC, names(PC)[grep("^pop$|population|pop_id|site", tolower(names(PC)))[1]], "pop")
setnames(ENV, names(ENV)[grep("^pop$|population|pop_id|site", tolower(names(ENV)))[1]], "pop")

# ---- find lat / long ----
cnl <- tolower(names(ENV))
lat_col  <- names(ENV)[grep("lat", cnl)[1]]
lon_col  <- names(ENV)[grep("lon|long", cnl)[1]]

cat("[info] lat =", lat_col, " lon =", lon_col, "\n")

df <- merge(
  PC,
  ENV[, .(pop, lat = get(lat_col), lon = get(lon_col))],
  by="pop", all=FALSE
)

pc_cols <- grep("^mitoPC[0-9]+$", names(df), value=TRUE)

res <- rbindlist(lapply(pc_cols, function(pc){
  fit <- lm(as.formula(paste(pc, "~ lat + lon")), data=df)
  sm  <- summary(fit)

  out <- as.data.table(coef(sm), keep.rownames="term")
  setnames(out, c("Estimate","Std. Error","t value","Pr(>|t|)"),
                c("beta","se","t","p"))

  out <- out[term %in% c("lat","lon")]
  out[, PC := pc]
  out[, R2 := sm$r.squared]
  out[, adjR2 := sm$adj.r.squared]
  out[, n := nobs(fit)]
  out[]
}), fill=TRUE)

setcolorder(res, c("PC","term","beta","se","t","p","R2","adjR2","n"))
fwrite(res, "/mnt/spareHD_2/nu_287/_assoc72_subunit/mtPC_noAMO_vs_lat_lon.tsv", sep="\t")

cat("\n=== smallest p per PC ===\n")
print(res[, .(min_p=min(p)), by=PC][order(min_p)])
RS






#af～ mtpc1+ treepc  run mt12345
Rscript /mnt/spareHD_2/nu_287/_assoc72_subunit/11_AF_vs_mitoPC1to5_plus_treePC12_noAMO.R

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

AF_FILE  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
PC_FILE  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"
COV_FILE <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR   <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_noAMO"
MIN_DEPTH <- 10L

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

OUT_LM  <- file.path(OUTDIR, "AF_LM_perSNP_mitoPC1to5_plus_treePC12.noAMO.tsv.gz")
OUT_SUM <- file.path(OUTDIR, "AF_LM_summary_by_mitoPC.noAMO.tsv")

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
# Read AF (noAMO)
# ----------------------------
cat("[read] AF: ", AF_FILE, "\n", sep="")
AF <- fread(AF_FILE)
need_af <- c("chr","pos","gene","pop","af","depth")
miss_af <- setdiff(need_af, names(AF))
if(length(miss_af)>0) stop("AF missing columns: ", paste(miss_af, collapse=", "))

AF[, pop := normalize_pop(pop)]
AF <- AF[pop != "AMO"]
AF <- AF[depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[info] AF rows after filters: ", nrow(AF), "\n", sep="")
cat("[info] unique pops: ", uniqueN(AF$pop), "\n", sep="")

# ----------------------------
# Read mitoPC (noAMO rebuilt tree)
# ----------------------------
cat("[read] PC: ", PC_FILE, "\n", sep="")
PC <- fread(PC_FILE)
pc_pop_col <- names(PC)[grep("^pop$|population|pop_id|site|id", tolower(names(PC)))[1]]
if(is.na(pc_pop_col)) stop("Cannot find pop column in PC file.")
setnames(PC, pc_pop_col, "pop")
PC[, pop := normalize_pop(pop)]
PC <- PC[pop != "AMO"]

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
COV <- COV[pop != "AMO"]

need_cov <- c("pop","treePC1","treePC2")
miss_cov <- setdiff(need_cov, names(COV))
if(length(miss_cov)>0) stop("COV missing columns: ", paste(miss_cov, collapse=", "))

COV <- unique(COV[, ..need_cov], by="pop")

# ----------------------------
# Merge predictors into AF
# ----------------------------
PRED <- merge(PC, COV, by="pop", all=FALSE)

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
# Per-SNP regressions: 5 separate models per SNP
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
  prop_beta_pos = mean(beta_mitoPC > 0, na.rm=TRUE)
), by=.(mitoPC)]

setorder(SUM, -prop_p05, med_p)
fwrite(SUM, OUT_SUM, sep="\t")
cat("[write] ", OUT_SUM, "\n", sep="")

cat("\n=== summary by mitoPC ===\n")
print(SUM)

cat("\n[OK] done. Output dir:\n  ", OUTDIR, "\n", sep="")

#lamda
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

in_file <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_noAMO/AF_LM_perSNP_mitoPC1to5_plus_treePC12.noAMO.tsv.gz"
out_dir <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mitoPC_noAMO/QQ_lambda"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

DT <- fread(cmd=paste("zcat", shQuote(in_file)))
DT <- DT[is.finite(p_mitoPC) & p_mitoPC>0 & p_mitoPC<=1]

# lambda helper: median(chi^2)/0.456
lambda_gc <- function(p){
  chisq <- qchisq(1 - p, df=1)
  median(chisq, na.rm=TRUE) / qchisq(0.5, df=1)
}

SUM <- DT[, .(
  n = .N,
  lambda = lambda_gc(p_mitoPC),
  min_p = min(p_mitoPC),
  med_p = median(p_mitoPC),
  prop_p05 = mean(p_mitoPC < 0.05)
), by=.(mitoPC)]

fwrite(SUM[order(mitoPC)], file.path(out_dir, "lambda_summary.tsv"), sep="\t")
print(SUM[order(mitoPC)])

# QQ plot per PC
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
       title="QQ plots: AF ~ mitoPCk + treePC1 + treePC2 (noAMO)")

ggsave(file.path(out_dir, "QQ_AF_mitoPC1to5.noAMO.png"), p, width=10.5, height=6.5, dpi=300)

Rscript /mnt/spareHD_2/nu_287/_assoc72_subunit/12_QQ_lambda_AF_mitoPC1to5.noAMO.R
    mitoPC     n      lambda        min_p      med_p   prop_p05
    <char> <int>       <num>        <num>      <num>      <num>
1: mitoPC1 42757  7.88180244 1.110319e-03 0.05827819 0.01683935
2: mitoPC2 42757 10.65457366 3.727465e-04 0.02769179 0.51198634
3: mitoPC3 42757  1.43821869 5.927136e-11 0.41858019 0.01562317
4: mitoPC4 42757  1.39643405 9.131170e-09 0.42542199 0.02493159
5: mitoPC5 42757  0.03513423 8.836412e-19 0.89939376 0.03052132







#af～ mtcluster
#/mnt/spareHD_2/nu_287/_assoc72_subunit/13_AF_vs_mtCluster_plus_treePC12_noAMO.R
#
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# noAMO | per-SNP regression:
#   af ~ mtCluster + treePC1 + treePC2
# Test mtCluster term by nested-model F-test:
#   reduced: af ~ treePC1 + treePC2
#   full:    af ~ mtCluster + treePC1 + treePC2
#
# Inputs:
#   AF long: chr,pos,gene,pop,af,depth
#   mtCluster manual: pop, mtCluster
#   treePC: pop, treePC1, treePC2
#
# Outputs:
#   OUTDIR/
#     AF_LM_perSNP_mtCluster_plus_treePC12.noAMO.tsv.gz
#     AF_LM_summary_mtCluster.noAMO.tsv
# ============================================================

AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR    <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AF_vs_mtCluster_noAMO"
MIN_DEPTH <- 10L

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

OUT_LM  <- file.path(OUTDIR, "AF_LM_perSNP_mtCluster_plus_treePC12.noAMO.tsv.gz")
OUT_SUM <- file.path(OUTDIR, "AF_LM_summary_mtCluster.noAMO.tsv")

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
# Read AF (noAMO)
# ----------------------------
cat("[read] AF: ", AF_FILE, "\n", sep="")
AF <- fread(AF_FILE)
need_af <- c("chr","pos","gene","pop","af","depth")
miss_af <- setdiff(need_af, names(AF))
if(length(miss_af)>0) stop("AF missing columns: ", paste(miss_af, collapse=", "))

AF[, pop := normalize_pop(pop)]
AF <- AF[pop != "AMO"]
AF <- AF[depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[info] AF rows after filters: ", nrow(AF), "\n", sep="")
cat("[info] unique pops: ", uniqueN(AF$pop), "\n", sep="")

# ----------------------------
# Read mtCluster manual
# ----------------------------
cat("[read] mtCluster: ", CL_FILE, "\n", sep="")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
CL <- CL[!(is.na(pop) | pop=="")]
need_cl <- c("pop","mtCluster")
miss_cl <- setdiff(need_cl, names(CL))
if(length(miss_cl)>0) stop("Cluster file missing columns: ", paste(miss_cl, collapse=", "))

CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

cat("[info] pops per cluster:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# ----------------------------
# Read treePC
# ----------------------------
cat("[read] treePC: ", COV_FILE, "\n", sep="")
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
COV <- COV[pop != "AMO"]
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

fwrite(LM, OUT_LM, sep="\t", compress="gzip")
cat("[write] ", OUT_LM, "\n", sep="")

# ----------------------------
# Summary + lambda
# ----------------------------
p <- LM$p_cluster
SUM <- data.table(
  n_snps = nrow(LM),
  min_p = suppressWarnings(min(p, na.rm=TRUE)),
  med_p = suppressWarnings(median(p, na.rm=TRUE)),
  prop_p05 = mean(p < 0.05, na.rm=TRUE),
  prop_p10 = mean(p < 0.10, na.rm=TRUE),
  lambda_gc = lambda_gc(p)
)

fwrite(SUM, OUT_SUM, sep="\t")
cat("[write] ", OUT_SUM, "\n", sep="")
print(SUM)

cat("\n[OK] done. Output dir:\n  ", OUTDIR, "\n", sep="")

cat /mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv
pop	mtCluster
SR	C1_AK
TL	C1_AK
WK	C1_AK
LB	C1_AK
MUC	C1_AK
SWA	C1_AK

PACH	C2_Recent
FRED	C2_Recent
BEA	C2_Recent
THE	C2_Recent

RS	C3_MarineLike
SAY	C3_MarineLike
SC	C3_MarineLike
CH	C3_MarineLike
ROB	C3_MarineLike
WB	C3_MarineLike
LG	C3_MarineLike
SL	C3_MarineLike
LAW	C3_MarineLike
BOOT	C3_MarineLike
JOE	C3_MarineLike
FG	C3_MarineLike

GOS	C4_GOS
PYE	C4_GOS
ECHO	C4_GOS
WT	C4_GOS
(base) cyu@stickleback:/mnt/spareHD_2/nu_287/q2_parallelism$ 

#gene level？ /mnt/spareHD_2/nu_287/_assoc72_subunit/15_geneLevel_AF_vs_mitoPC1to5_plus_treePC12.noAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

AF_FILE  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
PC_FILE  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"
COV_FILE <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/geneAF_vs_mitoPC_noAMO"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

MIN_DEPTH <- 10L

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

cat("[read] AF\n")
AF <- fread(AF_FILE)
stopifnot(all(c("gene","pop","af","depth") %in% names(AF)))
AF[, pop := normalize_pop(pop)]
AF <- AF[pop != "AMO" & depth >= MIN_DEPTH]

cat("[read] mitoPC\n")
PC <- fread(PC_FILE)
PC[, pop := normalize_pop(pop)]
PC <- PC[pop != "AMO"]

pc_cols <- grep("^mitoPC[1-5]$", names(PC), value=TRUE)
stopifnot(length(pc_cols) == 5)

PC <- unique(PC[, c("pop", pc_cols), with=FALSE], by="pop")

cat("[read] treePC\n")
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
COV <- COV[pop != "AMO"]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

PRED <- merge(merge(PC, COV, by="pop", all=FALSE), unique(AF[, .(pop)], by="pop"), by="pop", all=FALSE)

cat("[info] pops used: ", nrow(PRED), "\n", sep="")
cat("[info] pop list: ", paste(PRED$pop, collapse=","), "\n", sep="")

# ---------- gene-level AF per pop ----------
# mean AF across SNPs within each gene for each population
cat("[summarize] gene mean AF per pop\n")
GAF <- AF[, .(
  n_snps = .N,
  geneMeanAF = mean(af, na.rm=TRUE)
), by=.(gene, pop)]

# keep genes with enough SNPs across dataset (avoid tiny genes)
G_tot <- GAF[, .(n_snps_total = sum(n_snps)), by=gene]
GAF <- merge(GAF, G_tot, by="gene", all.x=TRUE)
GAF <- GAF[n_snps_total >= 20]   # 可调：20/50/100

# attach predictors
DT <- merge(GAF, PRED, by="pop", all=FALSE)

cat("[info] genes kept: ", uniqueN(DT$gene), "\n", sep="")
cat("[info] rows gene×pop: ", nrow(DT), "\n", sep="")

# ---------- per gene regressions ----------
run_gene_lm <- function(d, pc){
  # require enough pops to estimate
  if(uniqueN(d$pop) < 10) return(NULL)
  f <- as.formula(paste0("geneMeanAF ~ ", pc, " + treePC1 + treePC2"))
  fit <- try(lm(f, data=d), silent=TRUE)
  if(inherits(fit, "try-error")) return(NULL)
  sm <- summary(fit)
  co <- as.data.table(coef(sm), keep.rownames="term")
  setnames(co, c("Estimate","Std. Error","t value","Pr(>|t|)"),
              c("beta","se","t","p"))
  x <- co[term == pc]
  if(nrow(x)==0) return(NULL)
  x[, `:=`(R2=sm$r.squared, adjR2=sm$adj.r.squared, n=nobs(fit))]
  x[, mitoPC := pc]
  x[, .(beta, se, t, p, R2, adjR2, n, mitoPC)]
}

cat("[lm] per gene × mitoPC\n")
RES_list <- vector("list", length(pc_cols))
names(RES_list) <- pc_cols

for(pc in pc_cols){
  tmp <- DT[, {
    out <- run_gene_lm(.SD, pc)
    if(is.null(out)) NULL else out
  }, by=.(gene)]
  RES_list[[pc]] <- tmp
  cat("  ", pc, " genes: ", nrow(tmp), "\n", sep="")
}

RES <- rbindlist(RES_list, use.names=TRUE, fill=TRUE)

# BH within each PC
RES[, q := p.adjust(p, method="BH"), by=mitoPC]

# write full table
out_full <- file.path(OUTDIR, "geneLevel_AF_LM_byGene_mitoPC1to5_plus_treePC12.noAMO.tsv.gz")
fwrite(RES[order(mitoPC, q, p)], out_full, sep="\t")

# summary
SUM <- RES[, .(
  n_genes = .N,
  min_p = min(p, na.rm=TRUE),
  med_p = median(p, na.rm=TRUE),
  n_q05 = sum(q < 0.05, na.rm=TRUE),
  n_q10 = sum(q < 0.10, na.rm=TRUE)
), by=mitoPC][order(mitoPC)]

out_sum <- file.path(OUTDIR, "geneLevel_AF_summary_by_mitoPC.noAMO.tsv")
fwrite(SUM, out_sum, sep="\t")

cat("[write] ", out_full, "\n", sep="")
cat("[write] ", out_sum, "\n", sep="")
print(SUM)

cat("\n[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/geneAF_vs_mtCluster_noAMO"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

MIN_DEPTH <- 10L
MIN_SNP_GENE <- 20L   # 同你 mitoPC gene-level 脚本一致（可调）

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
# Read AF (noAMO) and compute geneMeanAF per pop
# ----------------------------
cat("[read] AF\n")
AF <- fread(AF_FILE)
stopifnot(all(c("gene","pop","af","depth") %in% names(AF)))

AF[, pop := normalize_pop(pop)]
AF <- AF[pop != "AMO" & depth >= MIN_DEPTH]

cat("[summarize] gene mean AF per pop\n")
GAF <- AF[, .(
  n_snps = .N,
  geneMeanAF = mean(af, na.rm=TRUE)
), by=.(gene, pop)]

# filter genes with enough SNPs (across all pops)
G_tot <- GAF[, .(n_snps_total = sum(n_snps)), by=gene]
GAF <- merge(GAF, G_tot, by="gene", all.x=TRUE)
GAF <- GAF[n_snps_total >= MIN_SNP_GENE]

cat("[info] genes kept: ", uniqueN(GAF$gene), "\n", sep="")
cat("[info] rows gene×pop: ", nrow(GAF), "\n", sep="")

# ----------------------------
# Read mtCluster manual
# ----------------------------
cat("[read] mtCluster\n")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
CL <- CL[!(is.na(pop) | pop=="")]
need_cl <- c("pop","mtCluster")
miss_cl <- setdiff(need_cl, names(CL))
if(length(miss_cl)>0) stop("Cluster file missing columns: ", paste(miss_cl, collapse=", "))

CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

cat("[info] pops per cluster:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# ----------------------------
# Read treePC
# ----------------------------
cat("[read] treePC\n")
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
COV <- COV[pop != "AMO"]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

# ----------------------------
# Merge predictors into gene-level table
# ----------------------------
PRED <- merge(COV, CL, by="pop", all=FALSE)
DT <- merge(GAF, PRED, by="pop", all=FALSE)

cat("[info] pops used: ", uniqueN(DT$pop), "\n", sep="")
cat("[info] genes used: ", uniqueN(DT$gene), "\n", sep="")

# ----------------------------
# Per-gene F-test for mtCluster term
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

# add SNP counts per gene
RES <- merge(RES, unique(GAF[, .(gene, n_snps_total)], by="gene"), by="gene", all.x=TRUE)

# write outputs
out_full <- file.path(OUTDIR, "geneLevel_AF_LM_byGene_mtCluster_plus_treePC12.noAMO.tsv.gz")
out_sum  <- file.path(OUTDIR, "geneLevel_AF_summary_mtCluster.noAMO.tsv")

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

# show top 10
cat("\n=== top 10 genes by p_cluster ===\n")
print(head(RES[order(p_cluster)], 10))

cat("\n[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")


gene level no



#residual
#/mnt/spareHD_2/nu_287/_assoc72_subunit/14_AFresid_vs_mtCluster_noAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/AFresid_vs_mtCluster_noAMO"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

OUT_LM  <- file.path(OUTDIR, "AFresid_LM_perSNP_mtCluster.noAMO.tsv.gz")
OUT_SUM <- file.path(OUTDIR, "AFresid_summary_mtCluster.noAMO.tsv")

MIN_DEPTH <- 10L

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

# ----------------------------
# Read AF
# ----------------------------
cat("[read] AF: ", AF_FILE, "\n", sep="")
AF <- fread(AF_FILE)
AF[, pop := normalize_pop(pop)]
AF <- AF[pop != "AMO" & depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]
cat("[info] AF rows: ", nrow(AF), " | pops: ", uniqueN(AF$pop), "\n", sep="")

# ----------------------------
# Read mtCluster (robust)
# ----------------------------
cat("[read] mtCluster: ", CL_FILE, "\n", sep="")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
CL <- CL[!(is.na(pop) | pop=="")]
CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]
cat("[info] cluster counts:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# ----------------------------
# Read treePC
# ----------------------------
cat("[read] treePC: ", COV_FILE, "\n", sep="")
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
COV <- COV[pop != "AMO"]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

# ----------------------------
# Merge
# ----------------------------
PRED <- merge(COV, CL, by="pop", all=FALSE)
DT <- merge(AF[, .(snp, chr, pos, gene, pop, af)], PRED, by="pop", all=FALSE)

cat("[info] merged rows: ", nrow(DT),
    " | SNPs: ", uniqueN(DT$snp),
    " | pops: ", uniqueN(DT$pop), "\n", sep="")

setkey(DT, snp)

# ----------------------------
# Per SNP: residualize then test mtCluster
# ----------------------------
cat("[lm] per SNP: af_resid = resid(lm(af ~ treePC1 + treePC2)); then af_resid ~ mtCluster\n")

RES <- DT[, {
  # data.table j 里不要 return()，用 NULL
  if(.N < 6) NULL else if(length(unique(mtCluster)) < 2) NULL else {
    fit1 <- try(lm(af ~ treePC1 + treePC2), silent=TRUE)
    if(inherits(fit1, "try-error")) NULL else {
      af_resid <- resid(fit1)
      fit2 <- try(lm(af_resid ~ mtCluster), silent=TRUE)
      if(inherits(fit2, "try-error")) NULL else {
        a <- anova(fit2)
        .(n=.N,
          df1=a$Df[1],
          df2=a$Df[2],
          F_cluster=a$`F value`[1],
          p_cluster=a$`Pr(>F)`[1])
      }
    }
  }
}, by=.(snp, chr, pos, gene)]

fwrite(RES, OUT_LM, sep="\t", compress="gzip")
cat("[write] ", OUT_LM, "\n", sep="")

p <- RES$p_cluster
SUM <- data.table(
  n_snps = nrow(RES),
  min_p = min(p, na.rm=TRUE),
  med_p = median(p, na.rm=TRUE),
  prop_p05 = mean(p < 0.05, na.rm=TRUE),
  prop_p10 = mean(p < 0.10, na.rm=TRUE),
  lambda_gc = lambda_gc(p)
)

fwrite(SUM, OUT_SUM, sep="\t")
cat("[write] ", OUT_SUM, "\n", sep="")
print(SUM)

cat("\n[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")





#/mnt/spareHD_2/nu_287/_assoc72_subunit/20_partialRDA_OXPHOS_AF_noAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

AF_FILE  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
PC_FILE  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"
COV_FILE <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/rda_mitoPC_noAMO"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

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
AF <- AF[pop != "AMO" & depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[read] mitoPC\n")
PC <- fread(PC_FILE)
PC[, pop := normalize_pop(pop)]
PC <- PC[pop != "AMO"]
pc_cols <- paste0("mitoPC", 1:5)
stopifnot(all(c("pop", pc_cols) %in% names(PC)))
PC <- unique(PC[, c("pop", pc_cols), with=FALSE], by="pop")

cat("[read] treePC\n")
COV <- fread(COV_FILE)
COV[, pop := normalize_pop(pop)]
COV <- COV[pop != "AMO"]
stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

# predictors table
PRED <- merge(PC, COV, by="pop", all=FALSE)
cat("[info] pops: ", nrow(PRED), "\n", sep="")

# ---------- build Y matrix: pop × SNP (then transpose) ----------
cat("[build] AF matrix\n")
DT <- merge(AF[, .(pop, snp, af)], PRED[, .(pop)], by="pop", all=FALSE)

# wide: rows=pop, cols=snp
Ypop <- dcast(DT, pop ~ snp, value.var="af")
pop_ids <- Ypop$pop
Y <- as.matrix(Ypop[, -"pop"])

# remove SNP columns with any NA (must be complete for vegan unless you impute)
keep <- colSums(!is.finite(Y)) == 0
cat("[info] SNPs total:", ncol(Y), " | complete:", sum(keep), "\n")
Y <- Y[, keep, drop=FALSE]

# Hellinger transform is standard for allele-frequency like community data
Yh <- decostand(Y, method="hellinger")

# align predictor rows to Y
PRED <- PRED[match(pop_ids, PRED$pop)]
stopifnot(all(PRED$pop == pop_ids))

X <- as.data.frame(PRED[, ..pc_cols])
Z <- as.data.frame(PRED[, .(treePC1, treePC2)])

# ---------- partial RDA ----------
cat("[rda] partial RDA: Y ~ mitoPCs + Condition(treePC1-2)\n")
mod <- rda(Yh ~ mitoPC1 + mitoPC2 + mitoPC3 + mitoPC4 + mitoPC5 + Condition(treePC1 + treePC2), data=cbind(X, Z))

# variance explained
R2  <- RsquareAdj(mod)$r.squared
aR2 <- RsquareAdj(mod)$adj.r.squared

# permutation tests
anova_global <- anova.cca(mod, permutations=9999)
anova_terms  <- anova.cca(mod, by="term", permutations=9999)
anova_axes   <- anova.cca(mod, by="axis", permutations=9999)

# save
saveRDS(mod, file.path(OUTDIR, "partialRDA_model.rds"))
fwrite(data.table(metric=c("R2","adjR2"), value=c(R2,aR2)),
       file.path(OUTDIR, "partialRDA_R2.tsv"), sep="\t")

capture.output(anova_global, file=file.path(OUTDIR, "partialRDA_anova_global.txt"))
capture.output(anova_terms,  file=file.path(OUTDIR, "partialRDA_anova_terms.txt"))
capture.output(anova_axes,   file=file.path(OUTDIR, "partialRDA_anova_axes.txt"))

cat("\n=== partial RDA summary ===\n")
print(data.table(R2=R2, adjR2=aR2))
print(anova_global)
print(anova_terms)

cat("\n[OK] OUTDIR:\n  ", OUTDIR, "\n", sep="")


#/mnt/spareHD_2/nu_287/_assoc72_subunit/21_partialRDA_AF_mtCluster_noAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

# =======================
# Inputs
# =======================
AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR    <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/rda_mtCluster_noAMO"
MIN_DEPTH <- 10L

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# =======================
# Helpers
# =======================
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

# =======================
# Read AF (noAMO) and build SNP id
# =======================
cat("[read] AF\n")
AF <- fread(AF_FILE)
need_af <- c("chr","pos","gene","pop","af","depth")
miss_af <- setdiff(need_af, names(AF))
if(length(miss_af)>0) stop("AF missing columns: ", paste(miss_af, collapse=", "))

AF[, pop := normalize_pop(pop)]
AF <- AF[pop != "AMO" & depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[info] AF rows: ", nrow(AF), " | pops: ", uniqueN(AF$pop), " | SNPs: ", uniqueN(AF$snp), "\n", sep="")

# =======================
# Read mtCluster (noAMO)
# =======================
cat("[read] mtCluster\n")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
CL <- CL[!(is.na(pop) | pop=="")]
if(!all(c("pop","mtCluster") %in% names(CL))) stop("Cluster file must have columns: pop, mtCluster")

CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

cat("[info] cluster counts:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# =======================
# Read treePC (noAMO)
# =======================
cat("[read] treePC\n")
COV <- fread(COV_FILE)
if(!all(c("pop","treePC1","treePC2") %in% names(COV))) stop("COV must have: pop, treePC1, treePC2")

COV[, pop := normalize_pop(pop)]
COV <- COV[pop != "AMO"]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

# =======================
# Build predictors table (meta)
# =======================
meta <- merge(COV, CL, by="pop", all=FALSE)
meta <- meta[order(pop)]

cat("[info] pops after merging predictors: ", nrow(meta), "\n", sep="")
cat("[info] pops: ", paste(meta$pop, collapse=","), "\n", sep="")

# =======================
# Build Y matrix: pop x SNP (complete cases)
# =======================
cat("[build] Y AF matrix\n")

DT <- merge(AF[, .(pop, snp, af)], meta[, .(pop)], by="pop", all=FALSE)

Ywide <- dcast(DT, pop ~ snp, value.var="af")
Ymat  <- as.matrix(Ywide[, -1, with=FALSE])
rownames(Ymat) <- Ywide$pop

# keep complete SNP columns only
ok_col <- colSums(!is.na(Ymat)) == nrow(Ymat)
Ymat2  <- Ymat[, ok_col, drop=FALSE]

cat("[info] SNPs total: ", ncol(Ymat), " | complete: ", ncol(Ymat2), "\n", sep="")

# reorder meta to match Y rows
meta2 <- meta[match(rownames(Ymat2), meta$pop)]
stopifnot(all(meta2$pop == rownames(Ymat2)))

# Hellinger transform is standard for allele-freq / community-like matrix
Y_hel <- decostand(Ymat2, method="hellinger")

# =======================
# partial RDA: Y ~ mtCluster + Condition(treePC1 + treePC2)
# =======================
cat("[rda] partial RDA: Y ~ mtCluster + Condition(treePC1-2)\n")

fit <- rda(
  Y_hel ~ mtCluster + Condition(treePC1 + treePC2),
  data = meta2
)

# global test + term test
a_global <- anova(fit, permutations=999)
a_terms  <- anova(fit, by="term", permutations=999)
r2       <- RsquareAdj(fit)

cat("\n=== anova(global) ===\n"); print(a_global)
cat("\n=== anova(by term) ===\n"); print(a_terms)
cat("\n=== R2 ===\n"); print(r2)

# =======================
# Write outputs
# =======================
out_txt <- file.path(OUTDIR, "partialRDA_mtCluster_noAMO.results.txt")
out_rds <- file.path(OUTDIR, "partialRDA_mtCluster_noAMO.fit.rds")

sink(out_txt)
cat("partial RDA: Y_hel ~ mtCluster + Condition(treePC1 + treePC2)\n\n")
cat("AF_FILE: ", AF_FILE, "\n", sep="")
cat("CL_FILE: ", CL_FILE, "\n", sep="")
cat("COV_FILE:", COV_FILE, "\n", sep="")
cat("MIN_DEPTH: ", MIN_DEPTH, "\n\n", sep="")
cat("Pops: ", nrow(meta2), "\n", sep="")
cat("SNPs complete: ", ncol(Y_hel), "\n\n", sep="")
cat("=== anova(global) ===\n"); print(a_global)
cat("\n=== anova(by term) ===\n"); print(a_terms)
cat("\n=== RsquareAdj ===\n"); print(r2)
sink()

saveRDS(list(
  fit=fit, meta=meta2, Y_hel=Y_hel,
  anova_global=a_global, anova_terms=a_terms, RsquareAdj=r2
), out_rds)

cat("[write] ", out_txt, "\n", sep="")
cat("[write] ", out_rds, "\n", sep="")
cat("[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")




#/mnt/spareHD_2/nu_287/_assoc72_subunit/22_partialdbRDA_AF_mtCluster_noAMO.R
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

# =======================
# Inputs
# =======================
AF_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
CL_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
COV_FILE  <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

OUTDIR    <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/dbrda_mtCluster_noAMO"
MIN_DEPTH <- 10L

# choose distance for dbRDA
DIST_METHOD <- "bray"   # "bray" is common; you can also try "euclidean" on Hellinger

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# =======================
# Helpers
# =======================
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

# =======================
# Read AF (noAMO) and build SNP id
# =======================
cat("[read] AF\n")
AF <- fread(AF_FILE)
need_af <- c("chr","pos","gene","pop","af","depth")
miss_af <- setdiff(need_af, names(AF))
if(length(miss_af)>0) stop("AF missing columns: ", paste(miss_af, collapse=", "))

AF[, pop := normalize_pop(pop)]
AF <- AF[pop != "AMO" & depth >= MIN_DEPTH]
AF[, snp := paste(chr, pos, gene, sep=":")]

cat("[info] AF rows: ", nrow(AF), " | pops: ", uniqueN(AF$pop), " | SNPs: ", uniqueN(AF$snp), "\n", sep="")

# =======================
# Read mtCluster (noAMO)
# =======================
cat("[read] mtCluster\n")
CL <- fread(CL_FILE, sep="\t", header=TRUE, fill=TRUE)
CL <- CL[!(is.na(pop) | pop=="")]
if(!all(c("pop","mtCluster") %in% names(CL))) stop("Cluster file must have columns: pop, mtCluster")

CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL[, .(pop, mtCluster)], by="pop")
CL[, mtCluster := factor(mtCluster)]

cat("[info] cluster counts:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

# =======================
# Read treePC (noAMO)
# =======================
cat("[read] treePC\n")
COV <- fread(COV_FILE)
if(!all(c("pop","treePC1","treePC2") %in% names(COV))) stop("COV must have: pop, treePC1, treePC2")

COV[, pop := normalize_pop(pop)]
COV <- COV[pop != "AMO"]
COV <- unique(COV[, .(pop, treePC1, treePC2)], by="pop")

# =======================
# Build predictors table (meta)
# =======================
meta <- merge(COV, CL, by="pop", all=FALSE)
meta <- meta[order(pop)]

cat("[info] pops after merging predictors: ", nrow(meta), "\n", sep="")

# =======================
# Build Y matrix: pop x SNP (complete cases)
# =======================
cat("[build] Y AF matrix\n")

DT <- merge(AF[, .(pop, snp, af)], meta[, .(pop)], by="pop", all=FALSE)

Ywide <- dcast(DT, pop ~ snp, value.var="af")
Ymat  <- as.matrix(Ywide[, -1, with=FALSE])
rownames(Ymat) <- Ywide$pop

ok_col <- colSums(!is.na(Ymat)) == nrow(Ymat)
Ymat2  <- Ymat[, ok_col, drop=FALSE]

cat("[info] SNPs total: ", ncol(Ymat), " | complete: ", ncol(Ymat2), "\n", sep="")

meta2 <- meta[match(rownames(Ymat2), meta$pop)]
stopifnot(all(meta2$pop == rownames(Ymat2)))

# Hellinger then distance
Y_hel <- decostand(Ymat2, method="hellinger")
D <- vegdist(Y_hel, method=DIST_METHOD)

# =======================
# partial dbRDA: D ~ mtCluster + Condition(treePC1 + treePC2)
# =======================
cat("[dbrda] capscale: D(Y) ~ mtCluster + Condition(treePC1-2)\n")

fit <- capscale(
  D ~ mtCluster + Condition(treePC1 + treePC2),
  data = meta2
)

a_global <- anova(fit, permutations=999)
a_terms  <- anova(fit, by="term", permutations=999)
r2       <- RsquareAdj(fit)

cat("\n=== anova(global) ===\n"); print(a_global)
cat("\n=== anova(by term) ===\n"); print(a_terms)
cat("\n=== R2 ===\n"); print(r2)

# =======================
# Write outputs
# =======================
out_txt <- file.path(OUTDIR, "partialdbRDA_mtCluster_noAMO.results.txt")
out_rds <- file.path(OUTDIR, "partialdbRDA_mtCluster_noAMO.fit.rds")

sink(out_txt)
cat("partial dbRDA (capscale): vegdist(Hellinger(Y)) ~ mtCluster + Condition(treePC1 + treePC2)\n\n")
cat("DIST_METHOD: ", DIST_METHOD, "\n", sep="")
cat("AF_FILE: ", AF_FILE, "\n", sep="")
cat("CL_FILE: ", CL_FILE, "\n", sep="")
cat("COV_FILE:", COV_FILE, "\n", sep="")
cat("MIN_DEPTH: ", MIN_DEPTH, "\n\n", sep="")
cat("Pops: ", nrow(meta2), "\n", sep="")
cat("SNPs complete: ", ncol(Y_hel), "\n\n", sep="")
cat("=== anova(global) ===\n"); print(a_global)
cat("\n=== anova(by term) ===\n"); print(a_terms)
cat("\n=== RsquareAdj ===\n"); print(r2)
sink()

saveRDS(list(
  fit=fit, meta=meta2, Y_hel=Y_hel, D=D,
  anova_global=a_global, anova_terms=a_terms, RsquareAdj=r2
), out_rds)

cat("[write] ", out_txt, "\n", sep="")
cat("[write] ", out_rds, "\n", sep="")
cat("[OK] done. OUTDIR:\n  ", OUTDIR, "\n", sep="")













#delta /mnt/spareHD_2/nu_287/q2_parallelism/
AF long (原料)
  └── af_long_final_72genes_subunit_with_si.tsv.gz
        |
        v
(1) 08_q2_deltaAF_SNPlevel_vs_mtPCk_plus_treePC12.R
        ├── deltaAF_long.noAMO.tsv.gz          [核心中间文件]
        ├── LM_perSNP_mitoPCk_plus_treePC12.tsv.gz
        └── LM_perSNP_treeOnly.tsv.gz
        |
        v
(2) 09_mtCluster_deltaAF_SNPlevel_LM_perm.R    [tree-based cluster + permutation]
        ├── LM_perSNP_mtCluster_plus_treePC12.tsv.gz
        ├── Perm_perSNP_mtCluster_plus_treePC12.tsv.gz
        └── Fig_*.png + Perm_global_summary.tsv
        |
        v
(3) 09_mtClusterManual_deltaAF_SNPlevel_perm0.R [manual cluster + driver_cluster]
        ├── LM_perSNP_mtCluster_manual_plus_treePC12.tsv.gz
        ├── LM_driverCluster_share_sigSNP.tsv
        └── Fig_*.png
        |
        v
(4) 10_gene_enrichment_from_cluster_LM.R       [BH + gene enrichment]
        └── gene_enrichment/*
              ├── LM_with_q_and_driver.tsv.gz
              ├── Gene_enrichment_overall.tsv.gz
              ├── GeneCluster_enrichment.tsv.gz
              └── Fig_topGenes_byCluster_*.png

(5) 06_q2_deltaAF_meanSD_vs_mtPC_noAMO.R       [pop-level meanΔ/sdΔ summary]
        ├── LM_meanDelta_vs_mitoPC1to5_plus_treePC12.tsv
        └── LM_sdDelta_vs_mitoPC1to5_plus_treePC12.tsv
             |
             v
      (你贴的 pop-level 整理脚本) -> LM_summary_all.sorted.tsv / minP / PC1_only

#/mnt/spareHD_2/nu_287/q2_parallelism/08_q2_deltaAF_SNPlevel_vs_mtPCk_plus_treePC12.R
#snp level
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ========= paths =========
af_file  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
pc_file  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"
cov_file <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"

out_dir  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_delta <- file.path(out_dir, "deltaAF_long.noAMO.tsv.gz")
out_mito  <- file.path(out_dir, "LM_perSNP_mitoPCk_plus_treePC12.tsv.gz")
out_tree  <- file.path(out_dir, "LM_perSNP_treeOnly.tsv.gz")
out_log   <- file.path(out_dir, "LOG_summary.txt")

# ========= params =========
min_depth <- 10
AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","AMO","BOOT","ECHO","LAW","GOS","ROB")
AK_marine <- "RS"
BC_marine <- "SAY"
min_n_AK <- 5
min_n_BC <- 8

normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}
infer_region <- function(pop){
  fifelse(pop %in% c(AK_fresh, AK_marine), "AK",
          fifelse(pop %in% c(BC_fresh, BC_marine), "BC", NA_character_))
}

# ========= fast lm helpers =========
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

fast_lm_tree_only <- function(y, t1, t2){
  ok <- is.finite(y) & is.finite(t1) & is.finite(t2)
  y <- y[ok]; t1 <- t1[ok]; t2 <- t2[ok]
  n <- length(y)
  if(n < 4) return(NULL)
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
AF <- fread(cmd = paste("zcat", shQuote(af_file)), sep = ",", header = TRUE)
stopifnot(all(c("chr","pos","gene","pop","af","depth") %in% names(AF)))
AF[, pop := normalize_pop(pop)]
AF <- AF[depth >= min_depth]
AF[, snp := paste(chr, pos, gene, sep=":")]

keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]

# drop AMO everywhere
AF <- AF[pop != "AMO"]

# add region now
AF[, region := infer_region(pop)]
AF <- AF[!is.na(region)]

# ========= read mitoPC + treePC =========
PC <- fread(pc_file)
PC[, pop := toupper(pop)]
PC <- PC[pop != "AMO"]

COV <- fread(cov_file)
COV[, pop := toupper(pop)]
COV <- COV[pop != "AMO"]

stopifnot(all(c("pop","treePC1","treePC2") %in% names(COV)))

pc_cols <- grep("^mitoPC[1-5]$", names(PC), value=TRUE)
if(length(pc_cols) < 1){
  stop("pc_file must contain mitoPC1..mitoPC5 at least one. Found: ",
       paste(names(PC), collapse=","))
}

PC  <- unique(PC[, c("pop", pc_cols), with=FALSE])
COV <- unique(COV[, .(pop, treePC1, treePC2)])

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

# attach SNP annotation
SNP_ANN <- unique(AF[, .(snp, chr, pos, gene)])
DT <- merge(DT, SNP_ANN, by="snp", all.x=TRUE)

# attach predictors
DT <- merge(DT, PRED, by="pop", all=FALSE)

# ---- HARD GUARANTEE region exists ----
if(!("region" %in% names(DT))){
  # if merge created region.x/region.y, rescue it
  if("region.x" %in% names(DT)) setnames(DT, "region.x", "region")
  if("region.y" %in% names(DT) && !("region" %in% names(DT))) setnames(DT, "region.y", "region")
}
if(!("region" %in% names(DT))){
  DT[, region := infer_region(pop)]
}
DT <- DT[!is.na(region)]

# sanity print
cat("[DEBUG] names(DT):\n")
print(names(DT))

stopifnot(!any(duplicated(names(DT))))

# region-specific minimum n per SNP
DT[, n_in_snp := .N, by=.(region, snp)]
DT <- DT[(region=="AK" & n_in_snp >= min_n_AK) | (region=="BC" & n_in_snp >= min_n_BC)]
DT[, n_in_snp := NULL]

# ========= write deltaAF long (SAFE column selection) =========
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

mito_res_list <- vector("list", length(pc_cols))
names(mito_res_list) <- pc_cols
for(pc in pc_cols){
  tmp <- DT[, {
    fit <- fast_lm_1x2t(deltaAF, get(pc), treePC1, treePC2)
    if(is.null(fit)) return(NULL)
    .(chr=chr[1], pos=pos[1], gene=gene[1],
      mitoPC=pc,
      beta_mitoPC=fit$beta,
      se_mitoPC=fit$se,
      t_mitoPC=fit$t,
      p_mitoPC=fit$p,
      R2=fit$R2,
      adjR2=fit$adjR2,
      n=fit$n)
  }, by=.(region, snp)]
  mito_res_list[[pc]] <- tmp
}
mito_res <- rbindlist(mito_res_list, use.names=TRUE, fill=TRUE)
setorder(mito_res, region, mitoPC, p_mitoPC)
fwrite(mito_res, out_mito, sep="\t", compress="gzip")

# ========= log =========
log_txt <- c(
  sprintf("[INFO] AF rows after depth/pop filter: %d", nrow(AF)),
  sprintf("[INFO] deltaAF observations: %d", nrow(DT)),
  sprintf("[INFO] mitoPC cols used: %s", paste(pc_cols, collapse=",")),
  sprintf("[INFO] tree-only rows: %d", nrow(tree_res)),
  sprintf("[INFO] mito rows: %d", nrow(mito_res)),
  "",
  "=== populations used per region ===",
  capture.output(DT[, .N, by=.(region, pop)][order(region, -N)])
)
writeLines(log_txt, out_log)

cat("[OK] wrote:\n")
cat(" - ", out_delta, "\n", sep="")
cat(" - ", out_tree,  "\n", sep="")
cat(" - ", out_mito,  "\n", sep="")
cat(" - ", out_log,   "\n", sep="")

Rscript /mnt/spareHD_2/nu_287/q2_parallelism/08_q2_deltaAF_SNPlevel_vs_mtPCk_plus_treePC12.R





library(data.table)

f <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/LM_perSNP_mitoPCk_plus_treePC12.tsv.gz"
DT <- fread(cmd=paste("zcat", shQuote(f)))

# 只保留有p的行
DT <- DT[p_mitoPC != "" & !is.na(p_mitoPC)]
DT[, p := as.numeric(p_mitoPC)]
DT[, b := as.numeric(beta_mitoPC)]

# summary：每个 region × mitoPC
S <- DT[, .(
  n = .N,
  min_p = min(p, na.rm=TRUE),
  med_p = median(p, na.rm=TRUE),
  prop_p05 = mean(p < 0.05, na.rm=TRUE),
  prop_p10 = mean(p < 0.10, na.rm=TRUE),
  mean_beta = mean(b, na.rm=TRUE),
  prop_beta_pos = mean(b > 0, na.rm=TRUE)
), by = .(region, mitoPC)]

setorder(S, region, -prop_p05, med_p)
print(S)

# 给你每个 region 的“最强PC”（按 p<0.05 富集）
TOP <- S[, .SD[1], by=region]
cat("\n=== top PC per region (by prop_p05) ===\n")
print(TOP)




#pop level
cat 06_q2_deltaAF_meanSD_vs_mtPC_noAMO.R 
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =======================
# Inputs
# =======================
af_file  <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"  # CSV gz: chr,pos,gene,pop,af,depth
pc_file  <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_noAMO_fromRebuiltTree.tsv"          # TSV: pop, mitoPC1..mitoPC5 (noAMO)
cov_file <- "/mnt/spareHD_2/nu_287/covariates.treePC.tsv"                                      # TSV: pop, treePC1, treePC2, ...

out_dir  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_vs_mtPC_noAMO"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_pop_sum <- file.path(out_dir, "deltaAF_pop_summary_mean_sd.tsv")
out_lm_mean <- file.path(out_dir, "LM_meanDelta_vs_mitoPC1to5_plus_treePC12.tsv")
out_lm_sd   <- file.path(out_dir, "LM_sdDelta_vs_mitoPC1to5_plus_treePC12.tsv")
fig_mean_pc1 <- file.path(out_dir, "Fig_meanDelta_vs_mitoPC1.png")
fig_sd_pc1   <- file.path(out_dir, "Fig_sdDelta_vs_mitoPC1.png")

# =======================
# Params
# =======================
min_depth <- 10

AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")  # 9
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS",
               "BOOT","ECHO","LAW","GOS","ROB")               # 12  ✅ AMO removed
AK_marine <- "RS"
BC_marine <- "SAY"

# =======================
# Helpers
# =======================
normalize_pop <- function(x){
  x <- toupper(x)
  x <- gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
  x
}

stopifnot(file.exists(af_file), file.exists(pc_file), file.exists(cov_file))

# =======================
# Read AF (CSV.gz)
# =======================
AF <- fread(cmd = paste("zcat", shQuote(af_file)), sep = ",", header = TRUE)
req_af <- c("chr","pos","gene","pop","af","depth")
if(!all(req_af %in% names(AF))){
  stop("AF file missing required columns. Need: ", paste(req_af, collapse=", "),
       "\nFound: ", paste(names(AF), collapse=", "))
}
AF[, pop := normalize_pop(pop)]
AF[, gene := tolower(gene)]
AF <- AF[depth >= min_depth]

# only keep relevant pops (no AMO here)
keep_pops <- unique(c(AK_fresh, AK_marine, BC_fresh, BC_marine))
AF <- AF[pop %in% keep_pops]

AF[, snp := paste(chr, pos, gene, sep=":")]

# =======================
# Build marine lookup per SNP
# =======================
marine_RS  <- AF[pop == AK_marine, .(snp, af_marine = af)]
marine_SAY <- AF[pop == BC_marine, .(snp, af_marine = af)]

# =======================
# Compute per-pop ΔAF within each region
#   ΔAF(pop) = AF_pop - AF_marine(region)
# =======================
calc_delta <- function(AF_dt, fresh_pops, marine_dt, region_label){
  sub <- AF_dt[pop %in% fresh_pops, .(snp, pop, af)]
  m <- merge(sub, marine_dt, by="snp", all.x=TRUE)
  m[, delta := af - af_marine]
  # summarize per population
  out <- m[is.finite(delta), .(
    region = region_label,
    n_snps = .N,
    mean_delta = mean(delta, na.rm=TRUE),
    sd_delta   = sd(delta,   na.rm=TRUE)
  ), by=.(pop)]
  out[]
}

AK_pop <- calc_delta(AF, AK_fresh, marine_RS,  "AK")
BC_pop <- calc_delta(AF, BC_fresh, marine_SAY, "BC")
POP_SUM <- rbindlist(list(AK_pop, BC_pop), use.names=TRUE, fill=TRUE)

# =======================
# Read mitoPC + treePC
# =======================
PC <- fread(pc_file)
PC[, pop := toupper(pop)]
need_pcs <- paste0("mitoPC", 1:5)
if(!all(need_pcs %in% names(PC))){
  stop("pc_file must contain mitoPC1..mitoPC5. Found: ", paste(names(PC), collapse=", "))
}

COV <- fread(cov_file)
COV[, pop := toupper(pop)]
if(!all(c("pop","treePC1","treePC2") %in% names(COV))){
  stop("cov_file must contain pop, treePC1, treePC2 at minimum. Found: ", paste(names(COV), collapse=", "))
}

# merge
df <- merge(POP_SUM, PC[, c("pop", need_pcs), with=FALSE], by="pop", all=FALSE)
df <- merge(df, COV[, .(pop, treePC1, treePC2)], by="pop", all=FALSE)

# sanity: should be AK=9, BC=12 (unless some pops lack AF/PC/COV)
cat("[info] populations after merge:\n")
print(df[, .N, by=region])

fwrite(df, out_pop_sum, sep="\t")

# =======================
# LM runner
# =======================
fit_lm_grid <- function(dat, response_col){
  out <- list()
  for(reg in sort(unique(dat$region))){
    dd <- dat[region == reg]
    for(k in 1:5){
      pc <- paste0("mitoPC", k)
      fml <- as.formula(paste0(response_col, " ~ ", pc, " + treePC1 + treePC2"))
      m <- lm(fml, data=dd)

      sm <- summary(m)
      co <- coef(sm)

      if(!(pc %in% rownames(co))) next

      out[[length(out)+1]] <- data.table(
        region = reg,
        response = response_col,
        mitoPC = pc,
        beta_mitoPC = unname(co[pc, "Estimate"]),
        se_mitoPC   = unname(co[pc, "Std. Error"]),
        t_mitoPC    = unname(co[pc, "t value"]),
        p_mitoPC    = unname(co[pc, "Pr(>|t|)"]),
        adjR2       = unname(sm$adj.r.squared),
        R2          = unname(sm$r.squared),
        n           = nrow(dd)
      )
    }
  }
  rbindlist(out, use.names=TRUE, fill=TRUE)
}

LM_mean <- fit_lm_grid(df, "mean_delta")
LM_sd   <- fit_lm_grid(df, "sd_delta")

setorder(LM_mean, region, mitoPC)
setorder(LM_sd,   region, mitoPC)

fwrite(LM_mean, out_lm_mean, sep="\t")
fwrite(LM_sd,   out_lm_sd,   sep="\t")

# =======================
# Figures: PC1 only (simple, for quick viewing)
# =======================
p_mean <- ggplot(df, aes(x=mitoPC1, y=mean_delta)) +
  geom_point(size=2, alpha=0.85) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~region, scales="free_y") +
  labs(x="mitoPC1 (noAMO)", y="mean ΔAF (pop vs regional marine)", title="Mean ΔAF vs mitoPC1 (noAMO; +treePC1-2 in LM)") +
  theme_classic(base_size = 14)

p_sd <- ggplot(df, aes(x=mitoPC1, y=sd_delta)) +
  geom_point(size=2, alpha=0.85) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~region, scales="free_y") +
  labs(x="mitoPC1 (noAMO)", y="SD(ΔAF) (pop vs regional marine)", title="SD(ΔAF) vs mitoPC1 (noAMO; +treePC1-2 in LM)") +
  theme_classic(base_size = 14)

ggsave(fig_mean_pc1, p_mean, width=7.2, height=4.5, dpi=300)
ggsave(fig_sd_pc1,   p_sd,   width=7.2, height=4.5, dpi=300)

cat("\n[OK] wrote:\n")
cat(" - ", out_pop_sum, "\n", sep="")
cat(" - ", out_lm_mean, "\n", sep="")
cat(" - ", out_lm_sd, "\n", sep="")
cat(" - ", fig_mean_pc1, "\n", sep="")
cat(" - ", fig_sd_pc1, "\n", sep="")






#ΔAF ~ mtCluster + treePC + permutation（tree-based cluster）
MT_TREE=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_noDloop_noAMO.iqtree.treefile \
N_PERM=200 \
Rscript /mnt/spareHD_2/nu_287/q2_parallelism/09_mtCluster_deltaAF_SNPlevel_LM_perm.R


#delta/af  ～  cluster treepc
#!/usr/bin/env Rscript
# ============================================================
# 09_mtCluster_deltaAF_SNPlevel_LM_perm.R
#
# SNP-level regression (within region AK/BC):
#   deltaAF ~ mtCluster + treePC1 + treePC2
# Test mtCluster term by nested-model F-test (anova(reduced, full)).
#
# Permutation test:
#   permute mtCluster labels among POPs *within each region*,
#   recompute F for each SNP, get empirical p_perm.
#
# Inputs (noAMO):
#   DELTA (long table): deltaAF_long.noAMO.tsv.gz
#     must contain: region, snp, chr, pos, gene, pop, deltaAF, treePC1, treePC2
#   MT_TREE: mt_noDloop_noAMO*.treefile  (non-ultrametric OK)
#
# Outputs:
#   OUTDIR/
#     mtCluster_by_pop.noAMO.tsv
#     LM_perSNP_mtCluster_plus_treePC12.tsv.gz
#     LM_cluster_summary_by_region.tsv
#     Perm_perSNP_mtCluster_plus_treePC12.tsv.gz   (if N_PERM>0)
#     Fig_pHist_mtCluster_LM.png
#     Fig_QQ_mtCluster_LM.png
#     Fig_pHist_mtCluster_perm.png (if N_PERM>0)
#     Perm_global_summary.tsv      (global permutation summary)
#
# Run:
#   Rscript 09_mtCluster_deltaAF_SNPlevel_LM_perm.R
# or:
#   MT_TREE=/path/to/treefile DELTA=/path/to/delta OUTDIR=/path/to/out \
#   K=4 N_PERM=200 SEED=1 Rscript 09_mtCluster_deltaAF_SNPlevel_LM_perm.R
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
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
)

OUT_DIR <- Sys.getenv("OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_noAMO"
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

# Find a treefile if not explicitly provided
find_treefile <- function(){
  # If user provides MT_TREE, use it.
  tf <- Sys.getenv("MT_TREE", unset = "")
  if(nzchar(tf) && file.exists(tf)) return(tf)

  # Otherwise search common locations (no recursive arg; base Sys.glob doesn't support it)
  candidates <- unique(c(
    Sys.glob("/work/cyu/poolseq/PPalign_output/**/mt*noAMO*.treefile"),
    Sys.glob("/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/*noAMO*.treefile"),
    Sys.glob("/mnt/spareHD_2/nu_287/_assoc72_subunit/*noAMO*.treefile"),
    Sys.glob("/mnt/spareHD_2/mt_gene_tree/*noAMO*.treefile")
  ))

  # Sys.glob doesn't support "**" portably; the first pattern may return none.
  # If none found, do shallow glob on a few directories:
  if(length(candidates) == 0){
    candidates <- unique(c(
      Sys.glob("/work/cyu/poolseq/PPalign_output/*/*noAMO*.treefile"),
      Sys.glob("/work/cyu/poolseq/PPalign_output/*/*/*noAMO*.treefile"),
      Sys.glob("/work/cyu/poolseq/PPalign_output/*/*/*/*noAMO*.treefile"),
      Sys.glob("/work/cyu/poolseq/PPalign_output/*/*/*/*/*noAMO*.treefile"),
      Sys.glob("/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/*noAMO*.treefile"),
      Sys.glob("/mnt/spareHD_2/nu_287/_assoc72_subunit/*noAMO*.treefile"),
      Sys.glob("/mnt/spareHD_2/mt_gene_tree/*noAMO*.treefile")
    ))
  }

  candidates <- candidates[file.exists(candidates)]
  if(length(candidates) == 0) return(NA_character_)

  # Prefer iqtree treefile if available
  pref <- candidates[grepl("iqtree", candidates, ignore.case=TRUE)]
  if(length(pref) > 0) return(pref[1])
  candidates[1]
}

stop_msg_tree <- function(){
  cat("\n[ERROR] Cannot find mt noAMO treefile.\n")
  cat("Set it explicitly like:\n")
  cat("  MT_TREE=/full/path/to/mt_noDloop_noAMO.iqtree.treefile Rscript 09_...\n\n")
  stop("Tree file not found.")
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
DEL <- DEL[pop != "AMO"]  # double safety

# Filter NA
DEL <- DEL[is.finite(deltaAF) & is.finite(treePC1) & is.finite(treePC2)]
cat("[info] rows after NA filter: ", nrow(DEL), "\n", sep="")

# ----------------------------
# Load mt tree -> mt clusters (non-ultrametric OK)
# ----------------------------
TREE_FILE <- find_treefile()
if(!is.character(TREE_FILE) || is.na(TREE_FILE) || !file.exists(TREE_FILE)) stop_msg_tree()

cat("[read tree] ", TREE_FILE, "\n", sep="")
tr <- read.tree(TREE_FILE)
cat("[tree] ntip = ", length(tr$tip.label), "\n", sep="")
cat("[tree] first tips: ", paste(head(tr$tip.label, 10), collapse=", "), "\n", sep="")

# Use patristic distances; no need ultrametric
D <- cophenetic(tr)  # matrix tip x tip
hc <- hclust(as.dist(D), method="average")
cl <- cutree(hc, k=K_CLUST)

mt_cluster <- data.table(
  pop = normalize_pop(names(cl)),
  mt_cluster = as.integer(cl)
)
mt_cluster <- unique(mt_cluster, by="pop")

# Keep only pops that appear in DEL (otherwise clusters for other tips unused)
pops_in_DEL <- unique(DEL$pop)
mt_cluster <- mt_cluster[pop %in% pops_in_DEL]

# Merge cluster
DEL <- merge(DEL, mt_cluster, by="pop", all.x=TRUE)

cat("\n[info] rows deltaAF: ", nrow(DEL), "\n", sep="")
cat("[info] unique pops: ", uniqueN(DEL$pop), "\n", sep="")
cat("[info] pops missing mt_cluster: ", sum(is.na(DEL$mt_cluster)), "\n", sep="")
if(sum(is.na(DEL$mt_cluster)) > 0){
  cat("[WARN] pops without mt_cluster (first 30):\n")
  print(head(unique(DEL[is.na(mt_cluster), pop]), 30))
}

# Drop missing cluster rows
DEL <- DEL[!is.na(mt_cluster)]
DEL[, cluster := factor(paste0("C", mt_cluster))]

# Save pop->cluster table (with region)
pop_region <- unique(DEL[, .(pop, region)], by="pop")
pop_cluster_out <- merge(pop_region, unique(DEL[, .(pop, mt_cluster)], by="pop"), by="pop", all.x=TRUE)
fwrite(pop_cluster_out[order(region, mt_cluster, pop)],
       file.path(OUT_DIR, "mtCluster_by_pop.noAMO.tsv"), sep="\t")
cat("[write] mtCluster_by_pop.noAMO.tsv\n")

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
    chr_  <- dt$chr[1]
    pos_  <- dt$pos[1]
    gene_ <- dt$gene[1]

    out <- lm_cluster_Ftest(dt)

    cm <- dt[, .(mu=mean(deltaAF, na.rm=TRUE), n=.N), by=cluster][order(cluster)]
    cm_str <- paste0(cm$cluster, ":", sprintf("%.6g", cm$mu), "(n=", cm$n, ")", collapse=";")

    idx <- idx + 1L
    lm_list[[idx]] <- data.table(
      region=r, snp=s, chr=chr_, pos=pos_, gene=gene_,
      n=out$n, n_clusters=length(unique(dt$cluster)),
      F_cluster=out$F, df1=out$df1, df2=out$df2, p_cluster=out$p,
      cluster_means=cm_str
    )
  }
}

LM <- rbindlist(lm_list[seq_len(idx)], use.names=TRUE, fill=TRUE)

lm_out <- file.path(OUT_DIR, "LM_perSNP_mtCluster_plus_treePC12.tsv.gz")
fwrite(LM, lm_out, sep="\t")
cat("\n[write] ", lm_out, "\n", sep="")

SUM <- LM[, .(
  n_snps = .N,
  prop_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
  prop_p10 = mean(p_cluster < 0.10, na.rm=TRUE),
  min_p = suppressWarnings(min(p_cluster, na.rm=TRUE)),
  med_p = suppressWarnings(median(p_cluster, na.rm=TRUE))
), by=region][order(region)]

fwrite(SUM, file.path(OUT_DIR, "LM_cluster_summary_by_region.tsv"), sep="\t")
cat("[write] LM_cluster_summary_by_region.tsv\n\n")
cat("=== summary (cluster term) ===\n")
print(SUM)

# ----------------------------
# Permutation test
# ----------------------------
perm_out_file <- file.path(OUT_DIR, "Perm_perSNP_mtCluster_plus_treePC12.tsv.gz")
perm_global_file <- file.path(OUT_DIR, "Perm_global_summary.tsv")

if(N_PERM <= 0){
  cat("\n[perm] N_PERM=0 -> skip permutation.\n")
} else {

  cat("\n[perm] running permutations: N_PERM=", N_PERM, "\n", sep="")

  # pop -> region -> cluster table
  pop_cluster_dt <- unique(DEL[, .(pop, region, cluster)], by=c("pop","region"))
  # observed F
  OBS <- LM[, .(region, snp, F_obs=F_cluster)]
  setkey(OBS, region, snp)
  OBS[, ge_count := 0L]
  OBS[, n_perm_used := 0L]

  # Minimal X
  BASE <- DEL[, .(region, snp, pop, deltaAF, treePC1, treePC2)]
  setkey(BASE, region, snp)

  # progress
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

        dt2 <- dt
        dt2[, cluster := cluster_perm]
        out <- lm_cluster_Ftest(dt2)
        if(!is.finite(out$F)) next

        f_obs <- OBS[list(r, s), F_obs]
        if(!is.finite(f_obs)) next

        # update counts
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

  # Global perm summary per region:
  G <- PERM[, .(
    n_snps=.N,
    prop_LM_p05 = mean(p_cluster < 0.05, na.rm=TRUE),
    prop_perm_p05 = mean(p_perm < 0.05, na.rm=TRUE),
    med_perm_p = suppressWarnings(median(p_perm, na.rm=TRUE)),
    n_perm_used_med = suppressWarnings(median(n_perm_used, na.rm=TRUE))
  ), by=region][order(region)]
  fwrite(G, perm_global_file, sep="\t")
  cat("[write] Perm_global_summary.tsv\n")

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
       title="SNP-level association: ΔAF ~ mtCluster + treePC1 + treePC2") +
  theme_classic(base_size=14)
ggsave(file.path(OUT_DIR, "Fig_pHist_mtCluster_LM.png"), p_hist, width=7.2, height=6.5, dpi=300)

LMq <- LM[is.finite(p_cluster) & p_cluster > 0 & p_cluster <= 1]
LMq[, obs := -log10(sort(p_cluster))]
LMq[, exp := -log10(ppoints(.N))]
p_qq <- ggplot(LMq, aes(exp, obs)) +
  geom_point(alpha=0.4, size=0.8) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~region) +
  labs(x="Expected -log10(p)", y="Observed -log10(p)",
       title="QQ plot: mtCluster term p-values") +
  theme_classic(base_size=14)
ggsave(file.path(OUT_DIR, "Fig_QQ_mtCluster_LM.png"), p_qq, width=7.2, height=4.2, dpi=300)

if(N_PERM > 0 && file.exists(perm_out_file)){
  PERM <- safe_zread(perm_out_file)
  if("p_perm" %in% names(PERM)){
    p_hist2 <- ggplot(PERM[is.finite(p_perm)], aes(p_perm)) +
      geom_histogram(bins=50) +
      facet_wrap(~region, ncol=1, scales="free_y") +
      labs(x="Permutation p-value for mtCluster term", y="count",
           title=paste0("Permutation test (N_PERM=", N_PERM, ")")) +
      theme_classic(base_size=14)
    ggsave(file.path(OUT_DIR, "Fig_pHist_mtCluster_perm.png"), p_hist2, width=7.2, height=6.5, dpi=300)
  }
}

# ----------------------------
# Top hits
# ----------------------------
cat("\n=== Top 20 SNPs by LM p_cluster (per region) ===\n")
LM_top <- LM[is.finite(p_cluster)][order(p_cluster)][, head(.SD, 20), by=region]
print(LM_top[, .(region, snp, gene, chr, pos, n, n_clusters, F_cluster, p_cluster)])

cat("\n[OK] done. Outputs in:\n  ", OUT_DIR, "\n", sep="")





#ΔAF ~ manual mtCluster + treePC（perm0 + driver_cluster）
#perm0 delta af～ cluster+treepc
#!/usr/bin/env Rscript
# ============================================================
# /mnt/spareHD_2/nu_287/q2_parallelism/09_mtClusterManual_deltaAF_SNPlevel_perm0.R
# deltaAF ~ mtCluster(manual) + treePC1 + treePC2
# per SNP, within region
# NO permutation (perm=0)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

DELTA_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
OUTDIR       <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

safe_zread <- function(f){
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
cat("[read] deltaAF_long\n")
DEL <- safe_zread(DELTA_FILE)
DEL[, pop := normalize_pop(pop)]
DEL <- DEL[pop != "AMO"]

cat("[read] mtCluster_manual\n")
# fill=TRUE: allow blank lines / ragged lines
CL <- fread(CLUSTER_FILE, sep="\t", header=TRUE, fill=TRUE)

# drop truly empty rows
CL <- CL[!(is.na(pop) | pop=="")]

# enforce two columns
need <- c("pop","mtCluster")
miss <- setdiff(need, names(CL))
if(length(miss) > 0){
  stop("Cluster file missing columns: ", paste(miss, collapse=", "))
}

CL[, pop := normalize_pop(pop)]
CL <- CL[pop != "AMO"]
CL <- unique(CL, by="pop")

# merge
DEL <- merge(DEL, CL, by = "pop", all.x = TRUE)
DEL <- DEL[!is.na(mtCluster)]
DEL[, mtCluster := factor(mtCluster)]

cat("[info] pops per cluster:\n")
print(as.data.table(table(CL$mtCluster))[order(-N)])

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
      mean_deltaAF = mean(deltaAF),
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

fwrite(LM,
       file.path(OUTDIR, "LM_perSNP_mtCluster_manual_plus_treePC12.tsv.gz"),
       sep = "\t")

# ----------------------------
# Summary
# ----------------------------
SUM <- LM[, .(
  n_snps = .N,
  prop_p05 = mean(p_cluster < 0.05, na.rm = TRUE),
  prop_p10 = mean(p_cluster < 0.10, na.rm = TRUE),
  min_p = suppressWarnings(min(p_cluster, na.rm = TRUE)),
  med_p = suppressWarnings(median(p_cluster, na.rm = TRUE))
), by = region][order(region)]

fwrite(SUM,
       file.path(OUTDIR, "LM_cluster_summary_by_region.tsv"),
       sep = "\t")

cat("\n=== summary ===\n")
print(SUM)

# ----------------------------
# Which cluster drives signal?
# ----------------------------
SIG <- LM[is.finite(p_cluster) & p_cluster < 0.05]
DRIVER <- SIG[, .N, by = .(region, driver_cluster)]
DRIVER[, prop := N / sum(N), by = region]
setorder(DRIVER, region, -prop)

fwrite(DRIVER,
       file.path(OUTDIR, "LM_driverCluster_share_sigSNP.tsv"),
       sep = "\t")

cat("\n=== driver cluster share among significant SNPs (p<0.05) ===\n")
print(DRIVER)

# ----------------------------
# Plots
# ----------------------------
p1 <- ggplot(LM[is.finite(p_cluster)], aes(p_cluster)) +
  geom_histogram(bins = 50) +
  facet_wrap(~region, ncol = 1) +
  theme_classic(base_size = 14) +
  labs(x="p-value (mtCluster F-test)", y="count",
       title="deltaAF ~ mtCluster(manual) + treePC1 + treePC2")

ggsave(file.path(OUTDIR, "Fig_pHist_mtCluster_manual.png"),
       p1, width = 7, height = 6, dpi = 300)

p2 <- ggplot(DRIVER, aes(driver_cluster, prop)) +
  geom_col() +
  facet_wrap(~region, ncol = 1, scales="free_x") +
  theme_classic(base_size = 14) +
  labs(x="driver cluster (largest |mean deltaAF|)", y="share among significant SNPs",
       title="Which mtCluster drives deltaAF signal?")

ggsave(file.path(OUTDIR, "Fig_driverCluster_share_sigSNP.png"),
       p2, width = 8, height = 5, dpi = 300)

cat("\n[OK] done. Outputs in:\n", OUTDIR, "\n", sep="")




#bh
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# Gene enrichment from mtCluster-driven SNPs (BH after LM)
#
# Inputs (from your manual-cluster 09 script):
#   OUTDIR contains:
#     - LM_perSNP_mtCluster_plus_treePC12.tsv.gz
#   (optional) If LM file doesn't contain driver_cluster, we will
#   infer it from cluster means if present; but best is: LM already
#   has driver_cluster column from your 09 script output.
#
# Outputs:
#   - LM_with_q_and_driver.tsv.gz
#   - Gene_enrichment_overall.tsv.gz
#   - GeneCluster_enrichment.tsv.gz
#   - Fig_topGenes_byCluster_<region>.png
# ============================================================

# ----------------------------
# Config (env override)
# ----------------------------
IN_LM <- Sys.getenv("IN_LM",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0/LM_perSNP_mtCluster_plus_treePC12.tsv.gz"
)

OUTDIR <- Sys.getenv("OUTDIR",
  unset = "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0/gene_enrichment"
)

# BH threshold for calling SNPs significant
Q_CUT <- as.numeric(Sys.getenv("Q_CUT", unset="0.10"))   # try 0.05 later

# Minimum SNPs per gene in the "universe" to test (avoid tiny genes)
MIN_SNP_GENE <- as.integer(Sys.getenv("MIN_SNP_GENE", unset="20"))

# Only consider these regions (auto if blank)
REGIONS <- Sys.getenv("REGIONS", unset="")  # e.g. "AK,BC" or "" for all

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
# If your LM has driver_cluster column, use it.
# If not, try to infer from cluster_means string if present.
if(!("driver_cluster" %in% names(LM))){
  if("cluster_means" %in% names(LM)){
    # cluster_means like: "C1:0.12(n=3);C2:-0.05(n=2)..."
    # Infer driver as max |mean - overall mean| approx by max |mean|
    infer_driver <- function(x){
      if(is.na(x) || x=="") return(NA_character_)
      parts <- unlist(strsplit(x, ";", fixed=TRUE))
      # extract "Ck:mean"
      cl <- sub("^([^:]+):.*$", "\\1", parts)
      mu <- sub("^[^:]+:([^\\(]+)\\(.*$", "\\1", parts)
      mu <- suppressWarnings(as.numeric(mu))
      if(all(!is.finite(mu))) return(NA_character_)
      cl[ which.max(abs(mu)) ]
    }
    LM[, driver_cluster := vapply(cluster_means, infer_driver, character(1))]
    cat("[warn] driver_cluster not found; inferred from cluster_means.\n")
  } else {
    cat("[ERROR] driver_cluster not found and cannot infer (no cluster_means).\n")
    stop("Please rerun 09 so LM contains driver_cluster (recommended).")
  }
}

# Save LM with q + driver
out_lm2 <- file.path(OUTDIR, "LM_with_q_and_driver.tsv.gz")
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

# precompute region totals
REG_TOT <- LM[, .(
  N_total=.N,
  N_sig=sum(is_sig)
), by=region]

# per gene counts
G <- LM[, .(
  n_gene = .N,
  n_sig_gene = sum(is_sig)
), by=.(region, gene)]

# filter tiny genes
G <- G[n_gene >= MIN_SNP_GENE]

# Fisher per gene
fisher_one <- function(a, b, c, d){
  # matrix [[a,b],[c,d]]; return OR and p
  m <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
  ft <- fisher.test(m)
  list(or=unname(ft$estimate), p=ft$p.value)
}

G[, `:=`(or_overall=NA_real_, p_overall=NA_real_)]

for(r in unique(G$region)){
  tot <- REG_TOT[region==r]
  # a = sig in gene
  # b = non-sig in gene
  # c = sig not in gene
  # d = non-sig not in gene
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

out_g_overall <- file.path(OUTDIR, "Gene_enrichment_overall.tsv.gz")
fwrite(G[order(region, q_overall, p_overall)],
       out_g_overall, sep="\t")
cat("[write] ", out_g_overall, "\n", sep="")

# ----------------------------
# 2) Gene × driver_cluster enrichment among significant SNPs
# Question: within significant SNPs, is driver_cluster overrepresented for a gene?
#
# For each region, gene, cluster:
#   a = # sig SNPs in gene with driver_cluster==cl
#   b = # sig SNPs in gene with driver_cluster!=cl
#   c = # sig SNPs not in gene with driver_cluster==cl
#   d = # sig SNPs not in gene with driver_cluster!=cl
# Fisher exact test.
# ----------------------------
cat("\n[run] gene × driver_cluster enrichment (within significant SNPs)\n")

SIG <- LM[is_sig == TRUE & !is.na(driver_cluster)]

# if no sig SNPs, stop gracefully
if(nrow(SIG)==0){
  cat("[WARN] No significant SNPs at q<", Q_CUT, ". Try Q_CUT=0.10 or 0.20.\n", sep="")
  quit(save="no", status=0)
}

# totals among sig SNPs
SIG_TOT <- SIG[, .N, by=.(region)]        # total sig per region
SIG_CL  <- SIG[, .N, by=.(region, driver_cluster)]  # sig per region×cluster

# per gene×cluster counts among sig
GC <- SIG[, .N, by=.(region, gene, driver_cluster)]
setnames(GC, "N", "n_sig_gene_cl")

# add gene total sig
GS <- SIG[, .N, by=.(region, gene)]
setnames(GS, "N", "n_sig_gene")
GC <- merge(GC, GS, by=c("region","gene"), all.x=TRUE)

# add region totals and region×cluster totals
GC <- merge(GC, SIG_TOT, by="region", all.x=TRUE)
setnames(GC, "N", "n_sig_region")
GC <- merge(GC, SIG_CL, by=c("region","driver_cluster"), all.x=TRUE)
setnames(GC, "N", "n_sig_region_cl")

# add universe size per gene (all SNPs) to allow filtering by MIN_SNP_GENE
GC <- merge(GC, G[, .(region, gene, n_gene)], by=c("region","gene"), all.x=TRUE)

GC <- GC[!is.na(n_gene) & n_gene >= MIN_SNP_GENE]

# fisher for each row
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

out_gc <- file.path(OUTDIR, "GeneCluster_enrichment.tsv.gz")
fwrite(GC[order(region, driver_cluster, q_fisher, p_fisher)],
       out_gc, sep="\t")
cat("[write] ", out_gc, "\n", sep="")

# ----------------------------
# Plot: top enriched genes per cluster per region
# ----------------------------
plot_top <- function(dt, region_name){
  dt <- dt[region == region_name]
  if(nrow(dt)==0) return(NULL)

  # keep only strong signals for plotting (raw p for ranking, but show q)
  # You can tweak these:
  dt2 <- dt[is.finite(p_fisher) & is.finite(or_fisher)]
  if(nrow(dt2)==0) return(NULL)

  # pick top 10 genes per cluster by smallest q (then p)
  dt2 <- dt2[order(q_fisher, p_fisher)]
  top <- dt2[, head(.SD, 10), by=driver_cluster]

  top[, gene := factor(gene, levels=rev(unique(gene[order(driver_cluster, q_fisher, p_fisher)])))]

  p <- ggplot(top, aes(x=gene, y=log2(or_fisher))) +
    geom_col() +
    coord_flip() +
    facet_wrap(~driver_cluster, scales="free_y") +
    labs(
      title=paste0("Top gene enrichments by mt driver_cluster (", region_name,
                   ")  |  significant SNPs: q_cluster<", Q_CUT),
      x="gene",
      y="log2(Fisher OR) for cluster among significant SNPs"
    ) +
    theme_classic(base_size=13)

  p
}

for(r in unique(GC$region)){
  p <- plot_top(GC, r)
  if(!is.null(p)){
    fn <- file.path(OUTDIR, paste0("Fig_topGenes_byCluster_", r, ".png"))
    ggsave(fn, p, width=12, height=7, dpi=300)
    cat("[write] ", fn, "\n", sep="")
  }
}

# ----------------------------
# Quick console summaries
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


IN_LM=/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_mtCluster_manual_perm0/LM_perSNP_mtCluster_manual_plus_treePC12.tsv.gz \
Rscript /mnt/spareHD_2/nu_287/q2_parallelism/10_gene_enrichment_from_cluster_LM.R







#A) AF level（推荐先跑）

问的是：不同 mt cluster 的群体，核基因 AF 是否系统性不同？（控制 treePC）

模型（每个 SNP，每个 region 分开做）

af∼mt_cluster+treePC1+treePC2

04_nuDeltaAF_SNPlevel_vs_mtDeltaAF_popSummary.R
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






#cluster
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