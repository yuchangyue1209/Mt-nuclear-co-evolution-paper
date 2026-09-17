#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

## ---------- 输入 ----------
mt_rds <- "/mnt/spareHD_2/nu_287/pinpis/results/Q1b_gene_final_pinpis_mt13.rds"
nu_rds <- "/mnt/spareHD_2/nu_287/pinpis/results/Q1b_gene_final_pinpis.rds"

outd <- "/mnt/spareHD_2/nu_287/pinpis/results/_figs_pi_mt_nu_like_dnds/complex_mt_nu_ass_stars"
dir.create(outd, recursive = TRUE, showWarnings = FALSE)

## ---------- 参数 ----------
signif_only <- TRUE
alpha <- 0.05
have_ragg <- requireNamespace("ragg", quietly = TRUE)

## ---------- 颜色/标签（与 dN/dS 一致） ----------
role_labels <- c(mt="mtOXPHOS", subunit="nuOXPHOS", assembly_factor="nuOXPHOS assembly_factor")
palette3    <- c(mt="#8CC6EC",  subunit="#F5A6A6", assembly_factor="#BFD99B")

## ---------- 读 RDS ----------
mt <- readRDS(mt_rds) %>% mutate(role="mt")
nu <- readRDS(nu_rds) %>% filter(role %in% c("subunit","assembly_factor"))

## 需要 complex：mt13 自带 complex；nu 的 complex 来自你 codeml annot join 后的 gene_final（你现在的 nu rds 里应该有 complex）
stopifnot("complex" %in% names(mt))
stopifnot("complex" %in% names(nu))

df <- bind_rows(
  mt %>% select(gene, role, complex, piN_mean, piS_mean, piN_piS_mean),
  nu %>% select(gene, role, complex, piN_mean, piS_mean, piN_piS_mean)
)

## ---------- Complex 统一 I–V ----------
df$complex <- toupper(trimws(df$complex))
df$Complex <- dplyr::recode(df$complex,
  "CI"="I","C I"="I","I"="I",
  "CII"="II","C II"="II","II"="II",
  "CIII"="III","C III"="III","III"="III",
  "CIV"="IV","C IV"="IV","IV"="IV",
  "CV"="V","C V"="V","V"="V",
  .default = NA_character_
)
df <- df[!is.na(df$Complex), , drop=FALSE]
df$Complex <- factor(df$Complex, levels=c("I","II","III","IV","V"))

## role 顺序固定
df$role <- factor(df$role, levels=c("mt","subunit","assembly_factor"))

## ---------- 保存 ----------
save_png_pdf <- function(g, file_base, w=7.6, h=4.8, dpi=300){
  if (have_ragg) { ragg::agg_png(paste0(file_base,".png"), w, h, units="in", res=dpi); print(g); dev.off() }
  else           { ggsave(paste0(file_base,".png"), g, width=w, height=h, dpi=dpi) }
  ggsave(paste0(file_base,".pdf"), g, width=w, height=h, dpi=dpi, device=cairo_pdf)
  message("[write] ", file_base, ".png / .pdf")
}

## ---------- position_dodge 下每组箱体 x 坐标 ----------
dodge_x <- function(x_index, role, roles_kept, dodge_width=0.72){
  n <- length(roles_kept)
  j <- match(role, roles_kept)
  offset <- ((j - (n+1)/2) / n) * dodge_width
  x_index + offset
}
p2star <- function(p){
  if (is.na(p)) "n.s."
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else "n.s."
}

## ---------- 生成一张（三组，每个 Complex 内 BH 后只标显著） ----------
make_plot3 <- function(dat, metric, ylab, file_stub){
  d <- dat[!is.na(dat[[metric]]), , drop=FALSE]
  if (nrow(d)==0) { warning("No data: ", file_stub); return(invisible(NULL)) }

  roles_keep <- c("mt","subunit","assembly_factor")
  dodge_w <- 0.72

  BOX_LWD=0.7; FRAME_LWD=0.6; PNT_SIZE=2.0; PNT_STROK=0.35

  g <- ggplot(d, aes(x=Complex, y=.data[[metric]], fill=role)) +
    geom_boxplot(position=position_dodge(width=dodge_w), width=0.6,
                 colour="black", linewidth=BOX_LWD, outlier.shape=NA) +
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

  ## ---- 每个 Complex 内两两 Wilcoxon + BH ----
  yr <- range(d[[metric]], na.rm=TRUE)
  span <- diff(yr); if (!is.finite(span) || span==0) span <- 1
  base_y <- yr[2] + 0.25*span
  step_y <- 0.12*span

  x_lvls <- levels(d$Complex); x_idx <- seq_along(x_lvls)
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
      if (length(x)<2 || length(y)<2) return(NA_real_)
      tryCatch(wilcox.test(x, y, exact=FALSE)$p.value, error=function(e) NA_real_)
    })
    padj <- if (length(pvals)>1) p.adjust(pvals, "BH") else pvals

    keep <- if (signif_only) which(!is.na(padj) & padj < alpha) else seq_along(padj)
    if (!length(keep)) next

    for (j in seq_along(keep)) {
      k <- keep[j]; pr <- prs[[k]]
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
      geom_segment(data=df_lines, aes(x=x, xend=xend, y=y, yend=y),
                   inherit.aes=FALSE, linewidth=0.6) +
      geom_text(data=df_text, aes(x=x, y=y, label=lab),
                inherit.aes=FALSE, size=4.4, vjust=0)
    g <- g + coord_cartesian(ylim=c(yr[1]-0.12*span, max(df_text$y)+0.12*span), clip="off")
  } else {
    g <- g + coord_cartesian(ylim=c(yr[1]-0.12*span, yr[2]+0.35*span), clip="off")
  }

  save_png_pdf(g, file.path(outd, file_stub))
}

## ---------- 出图：三张（πN / πS / πNπS） ----------
make_plot3(df, "piN_mean",     expression(pi[N]),       "PI_complex_mt_nu_ass_piN")
make_plot3(df, "piS_mean",     expression(pi[S]),       "PI_complex_mt_nu_ass_piS")
make_plot3(df, "piN_piS_mean", expression(pi[N]/pi[S]), "PI_complex_mt_nu_ass_piNpiS")

cat("Done. π complex figures -> ", outd, "\n")
