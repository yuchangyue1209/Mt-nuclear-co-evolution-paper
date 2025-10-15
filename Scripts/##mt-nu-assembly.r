##mt-nu-assembly
====== 路径 ======
tsv     <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/codeml_sites_summary.merged.tsv"
out_dir <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/_figs_mt_vs_subunit_R"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ====== 仅用于 dN/dS (omega) 的质量过滤 ======
min_dS_tree <- 0.001   # 可改 0.001 / 0.005 / 0.01

## ====== 选择：只标显著还是全部都标 ======
signif_only <- TRUE      # TRUE=只标显著(FDR<0.05)；FALSE=三条线全部标

## ====== 依赖包 ======
need_pkgs <- c("ggplot2","ggsignif")
for (p in need_pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
library(ggplot2)
library(ggsignif)
have_ragg <- requireNamespace("ragg", quietly = TRUE)

## ====== 读表 ======
df <- read.delim(tsv, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
numify <- function(x) suppressWarnings(as.numeric(x))
df$omega <- numify(df$omega)  # codeml M0 的 dN/dS
df$dN    <- numify(df$dN)
df$dS    <- numify(df$dS)

## ====== 三类角色 ======
keep_roles <- c("mt","subunit","assembly_factor")
df <- df[df$model=="M0" & df$role %in% keep_roles, , drop=FALSE]
df$role <- factor(df$role, levels = keep_roles)

## ====== 标签 & 颜色 ======
label_map <- c(mt="mtOXPHOS", subunit="nuOXPHOS", assembly_factor="nuOXPHOS assembly_factor")
fill_cols <- c(mt="#8CC6EC", subunit="#F5A6A6", assembly_factor="#BFD99B")

## ====== 保存函数 ======
.save_png_pdf <- function(g, file_base, width=7.0, height=4.8, dpi=300){
  if (have_ragg) { ragg::agg_png(paste0(file_base,".png"), width, height, units="in", res=dpi); print(g); dev.off() }
  else           { ggsave(paste0(file_base,".png"), g, width=width, height=height, dpi=dpi) }
  ggsave(paste0(file_base,".pdf"), g, width=width, height=height, dpi=dpi, device=cairo_pdf)
  message("[write] ", file_base, ".png / .pdf")
}

## ====== 成对 Wilcoxon + FDR；返回 ggsignif 所需对象 ======
pairwise_pvals <- function(df_plot, yvar){
  lv  <- levels(df_plot$role)
  prs <- combn(lv, 2, simplify = FALSE)
  raw_p <- sapply(prs, function(pr){
    x <- df_plot[df_plot$role==pr[1], yvar, drop=TRUE]
    y <- df_plot[df_plot$role==pr[2], yvar, drop=TRUE]
    if (length(x)<1 || length(y)<1) return(NA_real_)
    tryCatch(wilcox.test(x, y)$p.value, error=function(e) NA_real_)
  })
  adj_p <- p.adjust(raw_p, method="BH")
  comps <- lapply(prs, function(pr) c(label_map[pr[1]], label_map[pr[2]]))
  ann   <- vapply(adj_p, function(p) ifelse(is.na(p),"p=NA", ifelse(p<0.001,"p<0.001", sprintf("p=%.3f",p))), "")
  data.frame(
    x1 = vapply(comps, `[`, "", 1),
    x2 = vapply(comps, `[`, "", 2),
    annotation = ann,
    p_adj = adj_p,
    stringsAsFactors = FALSE
  )
}

## ====== 画图函数（支持只标显著 / 全部） ======
plot_three_with_brackets <- function(metric, ylab, out_stub, filter_for_omega=FALSE, alpha=0.05){
  dd <- df
  if (filter_for_omega) {
    dd <- dd[!is.na(dd$omega) & dd$omega < 999 & !is.na(dd$dS) & dd$dS >= min_dS_tree, ]
  }
  dd <- dd[!is.na(dd[[metric]]), , drop=FALSE]
  stopifnot(nrow(dd) > 0)

  # p 值表
  pw <- pairwise_pvals(dd, metric)
  if (signif_only) pw <- subset(pw, !is.na(p_adj) & p_adj < alpha)

  # y 轴留白
  rng  <- range(dd[[metric]], na.rm=TRUE)
  spn  <- diff(rng); if (!is.finite(spn) || spn==0) spn <- 1
  top0 <- rng[2] + 0.25*spn
  step <- 0.12*spn
  y_pos <- if (nrow(pw)) top0 + step * seq_len(nrow(pw)) else numeric(0)

  # 样式参数
  BOX_LWD=0.7; FRAME_LWD=0.6; PNT_SIZE=2.0; PNT_STROK=0.35

  g <- ggplot(dd, aes(role, .data[[metric]], fill=role)) +
    geom_boxplot(width=0.55, colour="black", size=BOX_LWD, outlier.shape=NA) +
    geom_point(position=position_jitter(width=0.12, height=0),
               size=PNT_SIZE, shape=21, stroke=PNT_STROK, colour="black") +
    scale_fill_manual(values=fill_cols, guide="none") +
    scale_x_discrete(labels=label_map) +
    labs(x=NULL, y=ylab) +
    theme_bw(base_size=14) +
    theme(
      text = element_text(family="sans"),
      axis.title = element_text(face="plain"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color="#e8e8e8"),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color="black", fill=NA, size=FRAME_LWD),
      axis.text.x = element_text(margin = margin(t=5)),
      plot.margin = unit(c(10,12,10,10), "pt")
    )

  if (nrow(pw)) {
    comps <- Map(function(a,b) c(a,b), pw$x1, pw$x2)
    g <- g + ggsignif::geom_signif(
      comparisons = comps,
      annotations = pw$annotation,
      y_position  = y_pos,
      tip_length  = 0.01,
      textsize    = 4.2,
      vjust       = 0.3,
      size        = 0.5
    ) +
    coord_cartesian(ylim = c(rng[1]-0.15*spn, max(y_pos)+0.2*spn), clip="off")
  } else {
    # 没有需要标的对：兜底显示总体 Kruskal–Wallis p
    kw_p <- tryCatch(kruskal.test(dd[[metric]] ~ dd$role)$p.value, error=function(e) NA_real_)
    txt  <- ifelse(is.na(kw_p),"", ifelse(kw_p<0.001,"p<0.001", sprintf("p=%.3f",kw_p)))
    g <- g +
      annotate("text", x=2, y=rng[2]+0.35*spn, label=txt, size=5) +
      coord_cartesian(ylim = c(rng[1]-0.15*spn, rng[2]+0.55*spn), clip="off")
  }

  .save_png_pdf(g, file.path(out_dir, out_stub))
}

## ====== 出图：dN / dS / dN/dS ======
plot_three_with_brackets("dN",    expression(d[N]),       "M0_dN_mt_nu_nuAF_pairs",     filter_for_omega=FALSE)
plot_three_with_brackets("dS",    expression(d[S]),       "M0_dS_mt_nu_nuAF_pairs",     filter_for_omega=FALSE)
plot_three_with_brackets("omega", expression(d[N]/d[S]),  "M0_dNdS_mt_nu_nuAF_pairs",   filter_for_omega=TRUE)

cat("All figures saved to:\n", out_dir, "\n")





#mt-nu
## ====== 路径 ======
tsv      <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/codeml_sites_summary.merged.tsv"
out_dir  <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/_figs_mt_vs_subunit_R"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ====== 仅用于 dN/dS (omega) 的质量过滤 ======
min_dS_tree <- 0.001

## ====== 只标显著还是总是标 ======
signif_only <- TRUE   # TRUE=只标显著(FDR<0.05)；FALSE=总是标星号

## ====== 依赖 ======
need_pkgs <- c("ggplot2","ggsignif")
for (p in need_pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
library(ggplot2)
library(ggsignif)
have_ragg <- requireNamespace("ragg", quietly = TRUE)

## ====== 读表 ======
df0 <- read.delim(tsv, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
numify <- function(x) suppressWarnings(as.numeric(x))
df0$omega <- numify(df0$omega)  # codeml M0 的 dN/dS
df0$dN    <- numify(df0$dN)
df0$dS    <- numify(df0$dS)

## ====== 只要两组：mt / subunit ======
df0 <- df0[df0$model=="M0" & df0$role %in% c("mt","subunit"), , drop=FALSE]
df0$role <- factor(df0$role, levels=c("mt","subunit"))

## ====== 标签 & 颜色（轴上显示）======
label_map <- c(mt="mtOXPHOS", subunit="nuOXPHOS")
fill_cols <- c(mt="#8CC6EC",  subunit="#F5A6A6")

## ====== 小工具 ======
.save_png_pdf <- function(g, file_base, width=6.6, height=4.6, dpi=300){
  if (have_ragg) { ragg::agg_png(paste0(file_base,".png"), width, height, units="in", res=dpi); print(g); dev.off() }
  else           { ggsave(paste0(file_base,".png"), g, width=width, height=height, dpi=dpi) }
  ggsave(paste0(file_base,".pdf"), g, width=width, height=height, dpi=dpi, device=cairo_pdf)
  message("[write] ", file_base, ".png / .pdf")
}
p_to_stars <- function(p){
  if (is.na(p)) "n.s."
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else "n.s."
}

## ====== 两组图（括号+星号）======
plot_two_with_star <- function(metric, ylab, out_stub, filter_for_omega=FALSE, alpha=0.05){
  dd <- df0
  if (filter_for_omega) {
    dd <- dd[!is.na(dd$omega) & dd$omega < 999 & !is.na(dd$dS) & dd$dS >= min_dS_tree, ]
  }
  dd <- dd[!is.na(dd[[metric]]), , drop=FALSE]
  stopifnot(nrow(dd) > 0)

  # 两组 Wilcoxon（exact=FALSE 避免 ties 告警）；两组时 FDR=原始 p
  p_raw <- tryCatch(wilcox.test(dd[[metric]] ~ dd$role, exact=FALSE)$p.value, error=function(e) NA_real_)
  p_adj <- p_raw
  star  <- p_to_stars(p_adj)

  # y 轴留白
  rng  <- range(dd[[metric]], na.rm=TRUE)
  spn  <- diff(rng); if (!is.finite(spn) || spn==0) spn <- 1
  ybar <- rng[2] + 0.30*spn
  ytxt <- rng[2] + 0.42*spn

  # 样式
  BOX_LWD=0.7; FRAME_LWD=0.6; PNT_SIZE=2.0; PNT_STROK=0.35

  g <- ggplot(dd, aes(role, .data[[metric]], fill=role)) +
    geom_boxplot(width=0.55, colour="black", size=BOX_LWD, outlier.shape=NA) +
    geom_point(position=position_jitter(width=0.12, height=0),
               size=PNT_SIZE, shape=21, stroke=PNT_STROK, colour="black") +
    scale_fill_manual(values=fill_cols, guide="none") +
    scale_x_discrete(labels=label_map) +
    labs(x=NULL, y=ylab) +
    theme_bw(base_size=14) +
    theme(
      text = element_text(family="sans"),
      axis.title = element_text(face="plain"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color="#e8e8e8"),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color="black", fill=NA, linewidth=FRAME_LWD),
      axis.text.x = element_text(margin=margin(t=5)),
      plot.margin = unit(c(10,12,10,10), "pt")
    )

  # 只标显著 or 总是标
  draw_it <- if (signif_only) (!is.na(p_adj) && p_adj < alpha) else TRUE
  if (draw_it) {
    g <- g +
      ggsignif::geom_signif(comparisons=list(c("mt","subunit")),
                            annotations=star, y_position=ytxt,
                            tip_length=0.01, textsize=5, vjust=0.3, size=0.6) +
      coord_cartesian(ylim = c(rng[1]-0.15*spn, ytxt+0.15*spn), clip="off")
  } else {
    g <- g + coord_cartesian(ylim = c(rng[1]-0.15*spn, rng[2]+0.35*spn), clip="off")
  }

  .save_png_pdf(g, file.path(out_dir, out_stub))
}

## ====== 出图：dN / dS / dN/dS（两组，带星号）======
plot_two_with_star("dN",    expression(d[N]),      "M0_dN_mt_vs_nuOXPHOS_star",   filter_for_omega=FALSE)
plot_two_with_star("dS",    expression(d[S]),      "M0_dS_mt_vs_nuOXPHOS_star",   filter_for_omega=FALSE)
plot_two_with_star("omega", expression(d[N]/d[S]), "M0_dNdS_mt_vs_nuOXPHOS_star", filter_for_omega=TRUE)

cat("Saved figures to:\n", out_dir, "\n")





#ars ribo
#!/usr/bin/env Rscript
## ===============================================================
## 6 张独立图（PNG+PDF）：
## cyto-ribo vs Nmt-ribo；cyto-ARS vs Nmt-ARS
## 指标：dN、dS、dN/dS(omega)；顶部括号+星号
## 直接使用 role 列进行分组
## ===============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
})

## ====== 路径 ======
tsv     <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/nu/codeml_sites_summary.with_cat.tsv"
out_dir <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/_figs_cyto_vs_Nmt_individual"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ====== 参数 ======
min_dS_tree <- 0.001     # 仅用于 omega 图的质量过滤
signif_only <- TRUE      # TRUE: 仅 p<0.05 画星号；FALSE: 总是画括号（星号可为 n.s.）

## ====== 读表 & 预处理 ======
df <- read.delim(tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
numify <- function(x) suppressWarnings(as.numeric(x))
if ("omega" %in% names(df)) df$omega <- numify(df$omega)
if ("dN"    %in% names(df)) df$dN    <- numify(df$dN)
if ("dS"    %in% names(df)) df$dS    <- numify(df$dS)

## 只保留 M0
if ("model" %in% names(df)) df <- df[df$model == "M0", , drop = FALSE]

## role 统一小写方便匹配，但保留原始值作为标签
df$role_raw <- df$role
df$role     <- tolower(df$role)

## 兼容命名：ribo/RP；nmt/nu
df$role <- gsub("^nu-", "nmt-", df$role)     # nu- → nmt-
df$role <- sub("rp$", "ribo", df$role)       # RP → ribo（末尾）
df$role <- sub("rp-", "ribo-", df$role)      # RP- → ribo-

## ====== 工具函数 ======
have_ragg <- requireNamespace("ragg", quietly = TRUE)
save_png_pdf <- function(g, file_base, width=6.4, height=4.6, dpi=300){
  if (have_ragg) { ragg::agg_png(paste0(file_base,".png"), width, height, units="in", res=dpi); print(g); dev.off() }
  else           { ggsave(paste0(file_base,".png"), g, width=width, height=height, dpi=dpi) }
  ggsave(paste0(file_base,".pdf"), g, width=width, height=height, dpi=dpi, device=cairo_pdf)
  message("[write] ", file_base, ".png / .pdf")
}
p_to_stars <- function(p){
  if (is.na(p)) "n.s."
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else "n.s."
}

plot_two <- function(dat, level_left, level_right, metric, ylab,
                     label_left, label_right, fill_left, fill_right,
                     out_stub, filter_for_omega=FALSE, alpha=0.05) {

  dd <- dat[dat$role %in% c(level_left, level_right), , drop=FALSE]

  # 质量过滤仅用于 ω
  if (filter_for_omega) {
    dd <- dd[!is.na(dd$omega) & dd$omega < 999 & !is.na(dd$dS) & dd$dS >= min_dS_tree, ]
  }
  # y 变量非缺失
  dd <- dd[!is.na(dd[[metric]]), , drop=FALSE]
  if (nrow(dd) < 2) {
    message("[skip] 有效数据不足：", out_stub,
            " | 已有水平：", paste(unique(dd$role), collapse=", "))
    return(invisible(NULL))
  }

  # 因子顺序（内部名给 ggsignif 用）
  dd$role <- factor(dd$role, levels = c(level_left, level_right))
  # 显示标签与配色
  disp_labels <- setNames(c(label_left, label_right), c(level_left, level_right))
  fills       <- setNames(c(fill_left,  fill_right),  c(level_left, level_right))

  # Wilcoxon（两组；ties -> exact=FALSE）
  p_raw <- tryCatch(wilcox.test(dd[[metric]] ~ dd$role, exact=FALSE)$p.value,
                    error=function(e) NA_real_)
  star  <- p_to_stars(p_raw)

  # y 轴留白
  rng  <- range(dd[[metric]], na.rm=TRUE)
  spn  <- diff(rng); if (!is.finite(spn) || spn==0) spn <- 1
  ytxt <- rng[2] + 0.40 * spn

  # 样式
  BOX_LWD=0.7; FRAME_LWD=0.6; PNT_SIZE=2.0; PNT_STROK=0.35

  g <- ggplot(dd, aes(role, .data[[metric]], fill=role)) +
    geom_boxplot(width=0.55, colour="black", size=BOX_LWD, outlier.shape=NA) +
    geom_point(position=position_jitter(width=0.12, height=0),
               size=PNT_SIZE, shape=21, stroke=PNT_STROK, colour="black") +
    scale_fill_manual(values=fills, guide="none") +
    scale_x_discrete(labels=disp_labels) +
    labs(x=NULL, y=ylab) +
    theme_bw(base_size=14) +
    theme(
      text=element_text(family="sans"),
      axis.title=element_text(face="plain"),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_line(color="#e8e8e8"),
      panel.grid.major.x=element_blank(),
      panel.border=element_rect(color="black", fill=NA, linewidth=FRAME_LWD),
      axis.text.x=element_text(margin=margin(t=5)),
      plot.margin=unit(c(10,12,10,10), "pt")
    )

  draw_it <- if (signif_only) (!is.na(p_raw) && p_raw < alpha) else TRUE
  if (draw_it) {
    g <- g + ggsignif::geom_signif(
      comparisons=list(c(level_left, level_right)),   # 用内部名
      annotations=star, y_position=ytxt,
      tip_length=0.01, textsize=5, vjust=0.3, size=0.6
    ) + coord_cartesian(ylim=c(rng[1]-0.15*spn, ytxt+0.15*spn), clip="off")
  } else {
    g <- g + coord_cartesian(ylim=c(rng[1]-0.15*spn, rng[2]+0.35*spn), clip="off")
  }

  save_png_pdf(g, file.path(out_dir, out_stub))
}

## ====== 颜色（ribo 与 ARS 用不同配色）======
COL_RIBO <- c(cyto="#3BA2D0", nmt="#F29F3D")   # ribo：蓝 / 橙
COL_ARS  <- c(cyto="#4C6FB3", nmt="#E06C9F")   # ARS ：紫蓝 / 洋红

## ====== 出图：共 6 张 ======
# ribo：支持 "cyto-ribo"/"nmt-ribo" 或 "cyto-rp"/"nmt-rp"（前面已做兼容）
plot_two(df, "cyto-ribo", "nmt-ribo",
         "dN", expression(d[N]),
         "cyto-RP", "Nmt-RP",
         COL_RIBO["cyto"], COL_RIBO["nmt"],
         "ribo_dN")

plot_two(df, "cyto-ribo", "nmt-ribo",
         "dS", expression(d[S]),
         "cyto-RP", "Nmt-RP",
         COL_RIBO["cyto"], COL_RIBO["nmt"],
         "ribo_dS")

plot_two(df, "cyto-ribo", "nmt-ribo",
         "omega", expression(d[N]/d[S]),
         "cyto-RP", "Nmt-RP",
         COL_RIBO["cyto"], COL_RIBO["nmt"],
         "ribo_dNdS", filter_for_omega = TRUE)

# ARS
plot_two(df, "cyto-ars", "nmt-ars",
         "dN", expression(d[N]),
         "cyto-ARS", "Nmt-ARS",
         COL_ARS["cyto"], COL_ARS["nmt"],
         "ARS_dN")

plot_two(df, "cyto-ars", "nmt-ars",
         "dS", expression(d[S]),
         "cyto-ARS", "Nmt-ARS",
         COL_ARS["cyto"], COL_ARS["nmt"],
         "ARS_dS")

plot_two(df, "cyto-ars", "nmt-ars",
         "omega", expression(d[N]/d[S]),
         "cyto-ARS", "Nmt-ARS",
         COL_ARS["cyto"], COL_ARS["nmt"],
         "ARS_dNdS", filter_for_omega = TRUE)

cat("Done. Figures -> ", out_dir, "\n")








#complex ass nu mt
#!/usr/bin/env Rscript
## ===============================================================
## 三组：mt vs nu vs nu-assembly_factor
## 每个 Complex 内两两比较（BH），只标显著
## 产出：dN / dS / dN/dS 三张图（PNG+PDF）
## ===============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

## ---------- 路径与参数 ----------
tsv  <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/codeml_sites_summary.merged.tsv"
outd <- "/mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/_figs_mt_vs_subunit_R/complex_mt_nu_ass_stars"
dir.create(outd, recursive = TRUE, showWarnings = FALSE)

min_dS_tree <- 0.001   # 仅用于 omega 的质量过滤
signif_only <- TRUE    # TRUE：只画显著；FALSE：全部都画
alpha       <- 0.05
have_ragg   <- requireNamespace("ragg", quietly = TRUE)

## ---------- 读表 ----------
df <- read.delim(tsv, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
numify <- function(x) suppressWarnings(as.numeric(x))
df$omega <- numify(df$omega); df$dN <- numify(df$dN); df$dS <- numify(df$dS)

df <- df[df$model=="M0" & df$role %in% c("mt","subunit","assembly_factor"), , drop=FALSE]

## Complex 统一 I–V，去 NA
df$complex <- toupper(trimws(df$complex))
df$Complex <- dplyr::recode(df$complex,
  "CI"="I","C I"="I","I"="I",
  "CII"="II","C II"="II","II"="II",
  "CIII"="III","C III"="III","III"="III",
  "CIV"="IV","C IV"="IV","IV"="IV",
  "CV"="V","C V"="V","V"="V",
  .default = NA_character_)
df <- df[!is.na(df$Complex), , drop=FALSE]
df$Complex <- factor(df$Complex, levels=c("I","II","III","IV","V"))

## 角色与配色
df$role <- factor(df$role, levels=c("mt","subunit","assembly_factor"))
role_labels <- c(mt="mtOXPHOS", subunit="nuOXPHOS", assembly_factor="nuOXPHOS assembly_factor")
palette3    <- c(mt="#8CC6EC", subunit="#F5A6A6", assembly_factor="#BFD99B")

## ---------- 保存 ----------
save_png_pdf <- function(g, file_base, w=7.6, h=4.8, dpi=300){
  if (have_ragg) { ragg::agg_png(paste0(file_base,".png"), w, h, units="in", res=dpi); print(g); dev.off() }
  else           { ggsave(paste0(file_base,".png"), g, width=w, height=h, dpi=dpi) }
  ggsave(paste0(file_base,".pdf"), g, width=w, height=h, dpi=dpi, device=cairo_pdf)
  message("[write] ", file_base, ".png / .pdf")
}

## ---------- position_dodge 下每组箱体的 x 坐标 ----------
dodge_x <- function(x_index, role, roles_kept, dodge_width=0.72){
  n <- length(roles_kept)
  j <- match(role, roles_kept)                  # 1..n
  offset <- ((j - (n+1)/2) / n) * dodge_width   # 与 position_dodge 对齐
  x_index + offset
}
p2star <- function(p){
  if (is.na(p)) "n.s."
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else "n.s."
}

## ---------- 生成一张“三组”图 ----------
make_plot3 <- function(dat, metric, ylab, file_stub, filter_omega=FALSE){
  d <- dat
  if (filter_omega) {
    d <- d[!is.na(d$omega) & d$omega < 999 & !is.na(d$dS) & d$dS >= min_dS_tree, ]
  }
  d <- d[!is.na(d[[metric]]), , drop=FALSE]
  if (nrow(d)==0) { warning("No data: ", file_stub); return(invisible(NULL)) }

  roles_keep <- c("mt","subunit","assembly_factor")
  dodge_w <- 0.72

  # 先画箱图
  BOX_LWD=0.7; FRAME_LWD=0.6; PNT_SIZE=2.0; PNT_STROK=0.35
  g <- ggplot(d, aes(x=Complex, y=.data[[metric]], fill=role)) +
    geom_boxplot(position=position_dodge(width=dodge_w), width=0.6,
                 colour="black", size=BOX_LWD, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.10, dodge.width=dodge_w),
               size=PNT_SIZE, shape=21, stroke=PNT_STROK, colour="black") +
    scale_fill_manual(values=palette3, labels=role_labels, name=NULL) +
    labs(x="Complex", y=ylab) +
    theme_bw(base_size=14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color="#e8e8e8"),
      panel.border = element_rect(color="black", fill=NA, linewidth=FRAME_LWD),
      legend.position = "top",
      plot.margin = unit(c(8,10,8,8), "pt")
    )

  # 逐 Complex 做 3 组两两比较（BH）
  yr <- range(d[[metric]], na.rm=TRUE); span <- diff(yr); if (!is.finite(span) || span==0) span <- 1
  base_y <- yr[2] + 0.25*span; step_y <- 0.12*span

  x_lvls  <- levels(d$Complex); x_idx <- seq_along(x_lvls)
  df_lines <- list(); df_text <- list(); row_id <- 1

  for (i in seq_along(x_lvls)) {
    cx <- x_lvls[i]
    di <- d[d$Complex==cx, , drop=FALSE]
    if (nrow(di) < 2) next

    roles_here <- roles_keep[roles_keep %in% unique(as.character(di$role))]
    if (length(roles_here) < 2) next

    prs <- combn(roles_here, 2, simplify=FALSE)
    pvals <- sapply(prs, function(pr){
      x <- di[di$role==pr[1], metric, drop=TRUE]
      y <- di[di$role==pr[2], metric, drop=TRUE]
      if (length(x)<1 || length(y)<1) return(NA_real_)
      tryCatch(wilcox.test(x, y, exact=FALSE)$p.value, error=function(e) NA_real_)
    })
    padj <- if (length(pvals)>1) p.adjust(pvals, "BH") else pvals
    keep <- if (signif_only) which(!is.na(padj) & padj < alpha) else seq_along(padj)
    if (!length(keep)) next

    for (j in seq_along(keep)) {
      k  <- keep[j]; pr <- prs[[k]]
      y0 <- base_y + (j-1)*step_y
      x1 <- dodge_x(x_idx[i], pr[1], roles_keep, dodge_w)
      x2 <- dodge_x(x_idx[i], pr[2], roles_keep, dodge_w)

      df_lines[[row_id]] <- data.frame(x=x1, xend=x2, y=y0, yend=y0)
      df_text [[row_id]] <- data.frame(x=(x1+x2)/2, y=y0 + 0.02*span, lab=p2star(padj[k]))
      row_id <- row_id + 1
    }
  }

  if (length(df_lines)) {
    df_lines <- do.call(rbind, df_lines)
    df_text  <- do.call(rbind, df_text)
    g <- g +
      geom_segment(data=df_lines, aes(x=x, xend=xend, y=y, yend=y), inherit.aes=FALSE, size=0.6) +
      geom_text(   data=df_text,  aes(x=x, y=y, label=lab), inherit.aes=FALSE, size=4.4, vjust=0)
    g <- g + coord_cartesian(ylim=c(yr[1]-0.12*span, max(df_text$y)+0.12*span), clip="off")
  } else {
    g <- g + coord_cartesian(ylim=c(yr[1]-0.12*span, yr[2]+0.35*span), clip="off")
  }

  save_png_pdf(g, file.path(outd, file_stub))
}

## ---------- 出图：三张 ----------
make_plot3(df, "dN",    expression(d[N]),      "mt_nu_ass_dN",    filter_omega=FALSE)
make_plot3(df, "dS",    expression(d[S]),      "mt_nu_ass_dS",    filter_omega=FALSE)
make_plot3(df, "omega", expression(d[N]/d[S]), "mt_nu_ass_dNdS",  filter_omega=TRUE)

cat("Done. Figures -> ", outd, "\n")
