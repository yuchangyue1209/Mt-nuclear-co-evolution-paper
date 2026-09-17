#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
  library(dplyr)
})

## ===== 输入 =====
nu_rds <- "/mnt/spareHD_2/nu_287/pinpis/results/Q1b_gene_final_pinpis.rds"
out_dir <- "/mnt/spareHD_2/nu_287/pinpis/results/_figs_pi_mt_nu_like_dnds/core_noncore"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

signif_only <- TRUE
alpha <- 0.05
have_ragg <- requireNamespace("ragg", quietly = TRUE)

## ===== 颜色（跟 dN/dS 一致的 nu 配色）=====
label_map <- c(nu_core="nu OXPHOS (core)", nu_noncore="nu OXPHOS (noncore)")
fill_cols <- c(nu_core="#F5A6A6", nu_noncore="#E5C16C")

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

## ===== 读 nu gene_final =====
df <- readRDS(nu_rds)

## 只保留 nu subunit（core/noncore 只对 subunit 有意义）
df <- df %>%
  filter(role == "subunit") %>%
  mutate(
    GENE_UP = toupper(gene),
    COMPLEX_UP = toupper(gsub("\\s+","", complex))
  )

has_suffix <- function(x, pat) grepl(paste0("^(", pat, ")[A-Z]?$"), x)

## CI 核心（NDUFS1-3/7/8；NDUFV1-2）
is_core_CI <- with(df, COMPLEX_UP=="CI" &
  ( has_suffix(GENE_UP, "NDUFS1|NDUFS2|NDUFS3|NDUFS7|NDUFS8") |
    has_suffix(GENE_UP, "NDUFV1|NDUFV2") ))

## CII 核心（SDHA-D）
is_core_CII <- with(df, COMPLEX_UP=="CII" &
  has_suffix(GENE_UP, "SDHA|SDHB|SDHC|SDHD"))

## CV 核心（兼容旧/新命名）
core_CV_F1 <- with(df, COMPLEX_UP=="CV" &
  ( has_suffix(GENE_UP, "ATP5A1|ATP5B|ATP5C1|ATP5D") |
    has_suffix(GENE_UP, "ATP5F1A|ATP5F1B|ATP5F1C|ATP5F1D") ))

core_CV_c <- with(df, COMPLEX_UP=="CV" &
  has_suffix(GENE_UP, "ATP5G1|ATP5G2|ATP5G3|ATP5MC1|ATP5MC2|ATP5MC3"))

is_nu_core <- is_core_CI | is_core_CII | core_CV_F1 | core_CV_c

df$grp <- ifelse(is_nu_core, "nu_core", "nu_noncore")
df$grp <- factor(df$grp, levels=c("nu_core","nu_noncore"))

message("[counts]"); print(table(df$grp))

## ===== 两组图函数（括号+星号）=====
plot_two <- function(metric, ylab, out_stub){
  dd <- df[!is.na(df[[metric]]), , drop=FALSE]
  stopifnot(nrow(dd) > 0)

  p_raw <- tryCatch(wilcox.test(dd[[metric]] ~ dd$grp, exact=FALSE)$p.value,
                    error=function(e) NA_real_)
  star <- p_to_stars(p_raw)

  rng <- range(dd[[metric]], na.rm=TRUE)
  spn <- diff(rng); if (!is.finite(spn) || spn==0) spn <- 1
  ytxt <- rng[2] + 0.42*spn

  BOX_LWD=0.7; FRAME_LWD=0.6; PNT_SIZE=2.0; PNT_STROK=0.35

  g <- ggplot(dd, aes(grp, .data[[metric]], fill=grp)) +
    geom_boxplot(width=0.55, colour="black", linewidth=BOX_LWD, outlier.shape=NA) +
    geom_point(position=position_jitter(width=0.12, height=0),
               size=PNT_SIZE, shape=21, stroke=PNT_STROK, colour="black") +
    scale_fill_manual(values=fill_cols, guide="none") +
    scale_x_discrete(labels=label_map) +
    labs(x=NULL, y=ylab) +
    theme_bw(base_size=14) +
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(color="#e8e8e8"),
      panel.border=element_rect(color="black", fill=NA, linewidth=FRAME_LWD),
      axis.text.x=element_text(margin=margin(t=5)),
      plot.margin=unit(c(10,12,10,10), "pt")
    )

  draw_it <- if (signif_only) (!is.na(p_raw) && p_raw < alpha) else TRUE
  if (draw_it) {
    g <- g +
      ggsignif::geom_signif(comparisons=list(c("nu_core","nu_noncore")),
                            annotations=star, y_position=ytxt,
                            tip_length=0.01, textsize=5, vjust=0.3, linewidth=0.6) +
      coord_cartesian(ylim=c(rng[1]-0.15*spn, ytxt+0.15*spn), clip="off")
  } else {
    g <- g + coord_cartesian(ylim=c(rng[1]-0.15*spn, rng[2]+0.35*spn), clip="off")
  }

  .save_png_pdf(g, file.path(out_dir, out_stub))
}

## ===== 出图：πN / πS / πNπS =====
plot_two("piN_mean",     expression(pi[N]),       "PI_nu_core_vs_noncore_piN")
plot_two("piS_mean",     expression(pi[S]),       "PI_nu_core_vs_noncore_piS")
plot_two("piN_piS_mean", expression(pi[N]/pi[S]), "PI_nu_core_vs_noncore_piNpiS")

cat("Done. Figures -> ", out_dir, "\n")
