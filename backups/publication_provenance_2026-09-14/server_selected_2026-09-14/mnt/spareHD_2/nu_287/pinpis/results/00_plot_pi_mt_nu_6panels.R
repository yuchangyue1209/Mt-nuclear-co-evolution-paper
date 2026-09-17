#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
  library(dplyr)
})

## ========= 输入 =========
mt_rds <- "/mnt/spareHD_2/nu_287/pinpis/results/Q1b_gene_final_pinpis_mt13.rds"
nu_rds <- "/mnt/spareHD_2/nu_287/pinpis/results/Q1b_gene_final_pinpis.rds"

out_dir <- "/mnt/spareHD_2/nu_287/pinpis/results/_figs_pi_mt_nu_like_dnds"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========= 参数 =========
alpha <- 0.05
signif_only <- TRUE  # TRUE: 只画显著星号；FALSE: 总是画（n.s.）
have_ragg <- requireNamespace("ragg", quietly = TRUE)

## ========= 颜色与标签（完全沿用你 dN/dS）=========
fill_cols_3 <- c(mt="#8CC6EC", subunit="#F5A6A6", assembly_factor="#BFD99B")
label_map_3 <- c(mt="mtOXPHOS", subunit="nuOXPHOS", assembly_factor="nuOXPHOS assembly_factor")

fill_cols_2 <- c(mt="#8CC6EC", subunit="#F5A6A6")
label_map_2 <- c(mt="mtOXPHOS", subunit="nuOXPHOS")

## ========= 读表 =========
mt <- readRDS(mt_rds)
nu <- readRDS(nu_rds)

## 只保留需要的 role
mt <- mt %>% mutate(role = "mt")
nu2 <- nu %>% filter(role %in% c("subunit","assembly_factor"))

## 合并成一个用于画图的数据框
df <- bind_rows(
  mt %>% select(gene, role, piN_mean, piS_mean, piN_piS_mean),
  nu2 %>% select(gene, role, piN_mean, piS_mean, piN_piS_mean)
)

## ========= 小工具 =========
.save_png_pdf <- function(g, file_base, width=7.0, height=4.8, dpi=300){
  if (have_ragg) {
    ragg::agg_png(paste0(file_base,".png"), width, height, units="in", res=dpi); print(g); dev.off()
  } else {
    ggsave(paste0(file_base,".png"), g, width=width, height=height, dpi=dpi)
  }
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

## ========= 三组：两两 Wilcoxon + BH（画星号）=========
pairwise_stars3 <- function(dd, yvar){
  lv <- levels(dd$role)
  prs <- combn(lv, 2, simplify=FALSE)

  raw_p <- sapply(prs, function(pr){
    x <- dd[dd$role==pr[1], yvar, drop=TRUE]
    y <- dd[dd$role==pr[2], yvar, drop=TRUE]
    if (length(x)<2 || length(y)<2) return(NA_real_)
    tryCatch(wilcox.test(x, y, exact=FALSE)$p.value, error=function(e) NA_real_)
  })
  adj_p <- p.adjust(raw_p, method="BH")

  out <- data.frame(
    a = vapply(prs, `[`, "", 1),
    b = vapply(prs, `[`, "", 2),
    p = adj_p,
    star = vapply(adj_p, p_to_stars, ""),
    stringsAsFactors = FALSE
  )

  if (signif_only) out <- out %>% filter(!is.na(p) & p < alpha)
  out
}

plot_three_with_stars <- function(metric, ylab, out_stub){
  dd <- df %>%
    filter(!is.na(.data[[metric]])) %>%
    mutate(role = factor(role, levels=c("mt","subunit","assembly_factor")))

  # 只保留三组都在的情况也行，但这里不强制
  dd <- dd %>% filter(role %in% c("mt","subunit","assembly_factor"))

  stopifnot(nrow(dd) > 0)

  pw <- pairwise_stars3(dd, metric)

  rng <- range(dd[[metric]], na.rm=TRUE)
  spn <- diff(rng); if (!is.finite(spn) || spn==0) spn <- 1
  top0 <- rng[2] + 0.25*spn
  step <- 0.12*spn
  y_pos <- if (nrow(pw)) top0 + step*seq_len(nrow(pw)) else numeric(0)

  g <- ggplot(dd, aes(role, .data[[metric]], fill=role)) +
    geom_boxplot(width=0.55, colour="black", linewidth=0.7, outlier.shape=NA) +
    geom_point(position=position_jitter(width=0.12, height=0),
               size=2.0, shape=21, stroke=0.35, colour="black") +
    scale_fill_manual(values=fill_cols_3, guide="none") +
    scale_x_discrete(labels=label_map_3) +
    labs(x=NULL, y=ylab) +
    theme_bw(base_size=14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color="#e8e8e8"),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color="black", fill=NA, linewidth=0.6),
      axis.text.x = element_text(margin=margin(t=5)),
      plot.margin = unit(c(10,12,10,10), "pt")
    )

  if (nrow(pw)) {
    comps <- Map(function(a,b) c(a,b), pw$a, pw$b)
    g <- g + ggsignif::geom_signif(
      comparisons = comps,
      annotations = pw$star,
      y_position  = y_pos,
      tip_length  = 0.01,
      textsize    = 5,
      vjust       = 0.25,
      size        = 0.6
    ) +
      coord_cartesian(ylim=c(rng[1]-0.15*spn, max(y_pos)+0.15*spn), clip="off")
  } else {
    g <- g + coord_cartesian(ylim=c(rng[1]-0.15*spn, rng[2]+0.35*spn), clip="off")
  }

  .save_png_pdf(g, file.path(out_dir, out_stub))
}

## ========= 两组：Wilcoxon（画星号）=========
plot_two_with_star <- function(metric, ylab, out_stub){
  dd <- df %>%
    filter(role %in% c("mt","subunit")) %>%
    filter(!is.na(.data[[metric]])) %>%
    mutate(role = factor(role, levels=c("mt","subunit")))

  stopifnot(nrow(dd) > 0)

  p_raw <- tryCatch(wilcox.test(dd[[metric]] ~ dd$role, exact=FALSE)$p.value, error=function(e) NA_real_)
  star <- p_to_stars(p_raw)
  draw_it <- if (signif_only) (!is.na(p_raw) && p_raw < alpha) else TRUE

  rng <- range(dd[[metric]], na.rm=TRUE)
  spn <- diff(rng); if (!is.finite(spn) || spn==0) spn <- 1
  ytxt <- rng[2] + 0.42*spn

  g <- ggplot(dd, aes(role, .data[[metric]], fill=role)) +
    geom_boxplot(width=0.55, colour="black", linewidth=0.7, outlier.shape=NA) +
    geom_point(position=position_jitter(width=0.12, height=0),
               size=2.0, shape=21, stroke=0.35, colour="black") +
    scale_fill_manual(values=fill_cols_2, guide="none") +
    scale_x_discrete(labels=label_map_2) +
    labs(x=NULL, y=ylab) +
    theme_bw(base_size=14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color="#e8e8e8"),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color="black", fill=NA, linewidth=0.6),
      axis.text.x = element_text(margin=margin(t=5)),
      plot.margin = unit(c(10,12,10,10), "pt")
    )

  if (draw_it) {
    g <- g + ggsignif::geom_signif(
      comparisons=list(c("mt","subunit")),
      annotations=star, y_position=ytxt,
      tip_length=0.01, textsize=6, vjust=0.25, size=0.7
    ) + coord_cartesian(ylim=c(rng[1]-0.15*spn, ytxt+0.15*spn), clip="off")
  } else {
    g <- g + coord_cartesian(ylim=c(rng[1]-0.15*spn, rng[2]+0.35*spn), clip="off")
  }

  .save_png_pdf(g, file.path(out_dir, out_stub))
}

## ========= 6 张图 =========
plot_three_with_stars("piN_mean",     expression(pi[N]),      "PI_gene_piN_mt_nu_ass_star")
plot_three_with_stars("piS_mean",     expression(pi[S]),      "PI_gene_piS_mt_nu_ass_star")
plot_three_with_stars("piN_piS_mean", expression(pi[N]/pi[S]),"PI_gene_piNpiS_mt_nu_ass_star")

plot_two_with_star("piN_mean",     expression(pi[N]),      "PI_gene_piN_mt_vs_nuOXPHOS_star")
plot_two_with_star("piS_mean",     expression(pi[S]),      "PI_gene_piS_mt_vs_nuOXPHOS_star")
plot_two_with_star("piN_piS_mean", expression(pi[N]/pi[S]),"PI_gene_piNpiS_mt_vs_nuOXPHOS_star")

cat("Done. 6 π figures -> ", out_dir, "\n")
