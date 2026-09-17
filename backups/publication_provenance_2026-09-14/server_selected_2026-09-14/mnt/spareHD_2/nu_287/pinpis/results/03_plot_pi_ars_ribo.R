#!/usr/bin/env Rscript
## ===============================================================
## 6 张独立图（PNG+PDF）：
## cyto-ribo vs Nmt-ribo；cyto-ARS vs Nmt-ARS
## 指标：piN / piS / piNpiS；顶部括号+星号
## 输入：Q1b_gene_final_pinpis.rds（gene-level 汇总表）
## ===============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
  library(dplyr)
})

## ====== 路径 ======
in_rds  <- "/mnt/spareHD_2/nu_287/pinpis/results/Q1b_gene_final_pinpis.rds"
out_dir <- "/mnt/spareHD_2/nu_287/pinpis/results/_figs_pi_mt_nu_like_dnds/ARS_Ribo"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ====== 参数 ======
signif_only <- TRUE      # TRUE: 仅 p<0.05 画星号；FALSE: 总是画括号（星号可为 n.s.）
alpha <- 0.05

## ====== 读表 ======
df <- readRDS(in_rds)

## ====== 从 role 推断分组：cyto-ars / nmt-ars / cyto-ribo / nmt-ribo ======
## 兼容你可能出现的命名：cyto-ARS, Nmt-ARS, cyto-RP, Nmt-RP, cyto-ribo, nmt-ribo,
## 也兼容 "cyto-ars"/"nmt-ars" 这种已经规范的。
to_grp <- function(role_chr){
  x <- tolower(trimws(role_chr))
  x <- gsub("_","-",x)

  is_ars  <- grepl("ars", x)
  is_ribo <- grepl("ribo|\\brp\\b|rp-", x)  # ribo / RP

  is_cyto <- grepl("^cyto|cytoplas", x)
  is_nmt  <- grepl("^nmt|mito|nuclear|nu-", x)  # 把 nu- 也归到 nmt 这边（核编码线粒体相关）

  if (is_ars && is_cyto) return("cyto-ars")
  if (is_ars && is_nmt)  return("nmt-ars")
  if (is_ribo && is_cyto) return("cyto-ribo")
  if (is_ribo && is_nmt)  return("nmt-ribo")
  return(NA_character_)
}

df <- df %>%
  mutate(
    grp = vapply(role, to_grp, character(1))
  ) %>%
  filter(!is.na(grp))

df$grp <- factor(df$grp, levels=c("cyto-ribo","nmt-ribo","cyto-ars","nmt-ars"))

message("[grp counts]"); print(table(df$grp))

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
                     out_stub) {

  dd <- dat[dat$grp %in% c(level_left, level_right), , drop=FALSE]
  dd <- dd[!is.na(dd[[metric]]), , drop=FALSE]
  if (nrow(dd) < 2) {
    message("[skip] 数据不足：", out_stub,
            " | levels present: ", paste(unique(dd$grp), collapse=", "))
    return(invisible(NULL))
  }

  dd$grp <- factor(dd$grp, levels=c(level_left, level_right))

  disp_labels <- setNames(c(label_left, label_right), c(level_left, level_right))
  fills       <- setNames(c(fill_left,  fill_right),  c(level_left, level_right))

  p_raw <- tryCatch(wilcox.test(dd[[metric]] ~ dd$grp, exact=FALSE)$p.value,
                    error=function(e) NA_real_)
  star <- p_to_stars(p_raw)

  rng <- range(dd[[metric]], na.rm=TRUE)
  spn <- diff(rng); if (!is.finite(spn) || spn==0) spn <- 1
  ytxt <- rng[2] + 0.40*spn

  BOX_LWD=0.7; FRAME_LWD=0.6; PNT_SIZE=2.0; PNT_STROK=0.35

  g <- ggplot(dd, aes(grp, .data[[metric]], fill=grp)) +
    geom_boxplot(width=0.55, colour="black", linewidth=BOX_LWD, outlier.shape=NA) +
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
    g <- g +
      ggsignif::geom_signif(
        comparisons=list(c(level_left, level_right)),
        annotations=star, y_position=ytxt,
        tip_length=0.01, textsize=5, vjust=0.3, linewidth=0.6
      ) +
      coord_cartesian(ylim=c(rng[1]-0.15*spn, ytxt+0.15*spn), clip="off")
  } else {
    g <- g + coord_cartesian(ylim=c(rng[1]-0.15*spn, rng[2]+0.35*spn), clip="off")
  }

  save_png_pdf(g, file.path(out_dir, out_stub))
}

## ====== 颜色（完全沿用你 dN/dS 的 ribo/ARS 配色）======
COL_RIBO <- c(cyto="#3BA2D0", nmt="#F29F3D")   # ribo：蓝 / 橙
COL_ARS  <- c(cyto="#4C6FB3", nmt="#E06C9F")   # ARS ：紫蓝 / 洋红

## ====== 出图：共 6 张 ======
# ribo
plot_two(df, "cyto-ribo", "nmt-ribo",
         "piN_mean", expression(pi[N]),
         "cyto-RP", "Nmt-RP",
         COL_RIBO["cyto"], COL_RIBO["nmt"],
         "PI_ribo_piN_cyto_vs_Nmt")

plot_two(df, "cyto-ribo", "nmt-ribo",
         "piS_mean", expression(pi[S]),
         "cyto-RP", "Nmt-RP",
         COL_RIBO["cyto"], COL_RIBO["nmt"],
         "PI_ribo_piS_cyto_vs_Nmt")

plot_two(df, "cyto-ribo", "nmt-ribo",
         "piN_piS_mean", expression(pi[N]/pi[S]),
         "cyto-RP", "Nmt-RP",
         COL_RIBO["cyto"], COL_RIBO["nmt"],
         "PI_ribo_piNpiS_cyto_vs_Nmt")

# ARS
plot_two(df, "cyto-ars", "nmt-ars",
         "piN_mean", expression(pi[N]),
         "cyto-ARS", "Nmt-ARS",
         COL_ARS["cyto"], COL_ARS["nmt"],
         "PI_ARS_piN_cyto_vs_Nmt")

plot_two(df, "cyto-ars", "nmt-ars",
         "piS_mean", expression(pi[S]),
         "cyto-ARS", "Nmt-ARS",
         COL_ARS["cyto"], COL_ARS["nmt"],
         "PI_ARS_piS_cyto_vs_Nmt")

plot_two(df, "cyto-ars", "nmt-ars",
         "piN_piS_mean", expression(pi[N]/pi[S]),
         "cyto-ARS", "Nmt-ARS",
         COL_ARS["cyto"], COL_ARS["nmt"],
         "PI_ARS_piNpiS_cyto_vs_Nmt")

cat("Done. Figures -> ", out_dir, "\n")
