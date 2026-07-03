#figure s1
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
  library(dplyr)
  library(patchwork)
})

## =========================================================
## INPUT PATHS
## =========================================================

mt_rds <- "/mnt/spareHD_2/nu_287/pinpis/results/Q1b_gene_final_pinpis_mt13.rds"
nu_rds <- "/mnt/spareHD_2/nu_287/pinpis/results/Q1b_gene_final_pinpis.rds"

out_dir <- "/mnt/spareHD_2/nu_287/pinpis/results/_figs_pnps_combined_big"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## PARAMETERS
## =========================================================

alpha <- 0.05
signif_only <- TRUE   # FALSE = 强制显示 n.s.

## =========================================================
## COLORS
## =========================================================

fill_cols_3 <- c(
  mt              = "#8CC6EC",
  subunit         = "#F5A6A6",
  assembly_factor = "#BFD99B"
)

label_map_3 <- c(
  mt              = "mtOXPHOS",
  subunit         = "nuOXPHOS",
  assembly_factor = "nuOXPHOS assembly factor"
)

fill_cols_2 <- c(
  mt      = "#8CC6EC",
  subunit = "#F5A6A6"
)

label_map_2 <- c(
  mt      = "mtOXPHOS",
  subunit = "nuOXPHOS"
)

fill_core <- c(
  nu_core    = "#F5A6A6",
  nu_noncore = "#E5C16C"
)

label_core <- c(
  nu_core    = "nuOXPHOS core",
  nu_noncore = "nuOXPHOS noncore"
)

COL_RIBO <- c(
  "cyto-ribo" = "#3BA2D0",
  "nmt-ribo"  = "#F29F3D"
)

COL_ARS <- c(
  "cyto-ars" = "#4C6FB3",
  "nmt-ars"  = "#E06C9F"
)

label_ribo <- c(
  "cyto-ribo" = "cyto-RP",
  "nmt-ribo"  = "Nmt-RP"
)

label_ars <- c(
  "cyto-ars" = "cyto-ARS",
  "nmt-ars"  = "Nmt-ARS"
)

## =========================================================
## READ DATA
## =========================================================

mt <- readRDS(mt_rds) %>%
  mutate(role = "mt")

nu <- readRDS(nu_rds)

cat("mt rows:", nrow(mt), "\n")
cat("nu rows:", nrow(nu), "\n")
cat("nu role counts:\n")
print(table(nu$role))

## =========================================================
## HELPERS
## =========================================================

p_to_stars <- function(p){
  if (is.na(p)) "n.s."
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else "n.s."
}

save_plot <- function(p, outfile, width = 6.5, height = 4.6){
  print(p)

  ggsave(paste0(outfile, ".pdf"), p, width = width, height = height)
  ggsave(paste0(outfile, ".png"), p, width = width, height = height, dpi = 300)

  message("[write] ", outfile, ".pdf / .png")
  return(p)
}

add_sig_two <- function(p, dd, xvar, yvar, comp){

  p_raw <- tryCatch(
    wilcox.test(dd[[yvar]] ~ dd[[xvar]], exact = FALSE)$p.value,
    error = function(e) NA_real_
  )

  star <- p_to_stars(p_raw)

  rng <- range(dd[[yvar]], na.rm = TRUE)
  spn <- diff(rng)
  if (!is.finite(spn) || spn == 0) spn <- 1

  ytxt <- rng[2] + 0.40 * spn

  draw_it <- if (signif_only) (!is.na(p_raw) && p_raw < alpha) else TRUE

  if (draw_it) {
    p <- p +
      geom_signif(
        comparisons = list(comp),
        annotations = star,
        y_position = ytxt,
        tip_length = 0.01,
        textsize = 5,
        vjust = 0.25,
        linewidth = 0.6
      ) +
      coord_cartesian(
        ylim = c(rng[1] - 0.15 * spn, ytxt + 0.15 * spn),
        clip = "off"
      )
  } else {
    p <- p +
      coord_cartesian(
        ylim = c(rng[1] - 0.15 * spn, rng[2] + 0.35 * spn),
        clip = "off"
      )
  }

  p
}

plot_box_two <- function(dd, xvar, yvar, ylab, cols, labels, comp, outfile){

  dd <- dd %>% filter(!is.na(.data[[yvar]]))

  p <- ggplot(dd, aes(.data[[xvar]], .data[[yvar]], fill = .data[[xvar]])) +
    geom_boxplot(width = 0.55, colour = "black", linewidth = 0.7, outlier.shape = NA) +
    geom_point(
      position = position_jitter(width = 0.12, height = 0),
      shape = 21, size = 2, stroke = 0.35, colour = "black"
    ) +
    scale_fill_manual(values = cols, labels = labels, guide = "none") +
    scale_x_discrete(labels = labels) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#e8e8e8"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.text.x = element_text(margin = margin(t = 5))
    )

  p <- add_sig_two(p, dd, xvar, yvar, comp)
  return(save_plot(p, file.path(out_dir, outfile)))
}

plot_box_three <- function(dd, yvar, ylab, outfile){

  dd <- dd %>%
    filter(role %in% c("mt", "subunit", "assembly_factor")) %>%
    filter(!is.na(.data[[yvar]])) %>%
    mutate(role = factor(role, levels = c("mt", "subunit", "assembly_factor")))

  prs <- combn(levels(dd$role), 2, simplify = FALSE)

  pvals <- sapply(prs, function(pr){
    x <- dd[dd$role == pr[1], yvar, drop = TRUE]
    y <- dd[dd$role == pr[2], yvar, drop = TRUE]
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
  })

  padj <- p.adjust(pvals, method = "BH")
  stars <- vapply(padj, p_to_stars, character(1))
  keep <- if (signif_only) which(!is.na(padj) & padj < alpha) else seq_along(padj)

  rng <- range(dd[[yvar]], na.rm = TRUE)
  spn <- diff(rng)
  if (!is.finite(spn) || spn == 0) spn <- 1

  p <- ggplot(dd, aes(role, .data[[yvar]], fill = role)) +
    geom_boxplot(width = 0.55, colour = "black", linewidth = 0.7, outlier.shape = NA) +
    geom_point(
      position = position_jitter(width = 0.12, height = 0),
      shape = 21, size = 2, stroke = 0.35, colour = "black"
    ) +
    scale_fill_manual(values = fill_cols_3, guide = "none") +
    scale_x_discrete(labels = label_map_3) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#e8e8e8"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.text.x = element_text(margin = margin(t = 5))
    )

  if (length(keep) > 0) {
    comps <- prs[keep]
    ypos <- rng[2] + 0.25 * spn + seq_along(keep) * 0.12 * spn

    p <- p +
      geom_signif(
        comparisons = comps,
        annotations = stars[keep],
        y_position = ypos,
        tip_length = 0.01,
        textsize = 5,
        vjust = 0.25,
        linewidth = 0.6
      ) +
      coord_cartesian(
        ylim = c(rng[1] - 0.15 * spn, max(ypos) + 0.15 * spn),
        clip = "off"
      )
  } else {
    p <- p +
      coord_cartesian(
        ylim = c(rng[1] - 0.15 * spn, rng[2] + 0.35 * spn),
        clip = "off"
      )
  }

  return(save_plot(p, file.path(out_dir, outfile)))
}

## =========================================================
## OXPHOS DATA
## =========================================================

df_ox <- bind_rows(
  mt %>% select(gene, role, complex, piN_mean, piS_mean, piN_piS_mean),
  nu %>%
    filter(role %in% c("subunit", "assembly_factor")) %>%
    select(gene, role, complex, piN_mean, piS_mean, piN_piS_mean)
)

## =========================================================
## A: mtOXPHOS / nuOXPHOS / assembly
## =========================================================

pA1 <- plot_box_three(df_ox, "piN_mean", expression(pi[N]), "A_OXPHOS_mt_nu_ass_piN")
pA2 <- plot_box_three(df_ox, "piS_mean", expression(pi[S]), "A_OXPHOS_mt_nu_ass_piS")
pA3 <- plot_box_three(df_ox, "piN_piS_mean", expression(pi[N]/pi[S]), "A_OXPHOS_mt_nu_ass_piNpiS")

## =========================================================
## B: mtOXPHOS vs nuOXPHOS
## =========================================================

df_mt_nu <- df_ox %>%
  filter(role %in% c("mt", "subunit")) %>%
  mutate(role = factor(role, levels = c("mt", "subunit")))

pB1 <- plot_box_two(df_mt_nu, "role", "piN_mean", expression(pi[N]),
                    fill_cols_2, label_map_2, c("mt", "subunit"),
                    "B_OXPHOS_mt_vs_nu_piN")

pB2 <- plot_box_two(df_mt_nu, "role", "piS_mean", expression(pi[S]),
                    fill_cols_2, label_map_2, c("mt", "subunit"),
                    "B_OXPHOS_mt_vs_nu_piS")

pB3 <- plot_box_two(df_mt_nu, "role", "piN_piS_mean", expression(pi[N]/pi[S]),
                    fill_cols_2, label_map_2, c("mt", "subunit"),
                    "B_OXPHOS_mt_vs_nu_piNpiS")

## =========================================================
## C: Complex-level mt / nu / assembly WITH STARS
## =========================================================

df_complex <- df_ox %>%
  mutate(
    complex = toupper(trimws(complex)),
    Complex = recode(
      complex,
      "CI" = "I", "C I" = "I", "I" = "I",
      "CII" = "II", "C II" = "II", "II" = "II",
      "CIII" = "III", "C III" = "III", "III" = "III",
      "CIV" = "IV", "C IV" = "IV", "IV" = "IV",
      "CV" = "V", "C V" = "V", "V" = "V",
      .default = NA_character_
    )
  ) %>%
  filter(!is.na(Complex)) %>%
  mutate(
    Complex = factor(Complex, levels = c("I", "II", "III", "IV", "V")),
    role = factor(role, levels = c("mt", "subunit", "assembly_factor"))
  )

dodge_x <- function(x_index, role, roles_kept, dodge_width = 0.72){
  n <- length(roles_kept)
  j <- match(role, roles_kept)
  offset <- ((j - (n + 1) / 2) / n) * dodge_width
  x_index + offset
}

plot_complex <- function(yvar, ylab, outfile){

  dd <- df_complex %>% filter(!is.na(.data[[yvar]]))

  roles_keep <- c("mt", "subunit", "assembly_factor")
  dodge_w <- 0.72

  p <- ggplot(dd, aes(Complex, .data[[yvar]], fill = role)) +
    geom_boxplot(
      position = position_dodge(width = dodge_w),
      width = 0.6,
      colour = "black",
      linewidth = 0.7,
      outlier.shape = NA
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.10, dodge.width = dodge_w),
      shape = 21,
      size = 2,
      stroke = 0.35,
      colour = "black"
    ) +
    scale_fill_manual(values = fill_cols_3, labels = label_map_3, name = NULL) +
    labs(x = "Complex", y = ylab) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#e8e8e8"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )

  yr <- range(dd[[yvar]], na.rm = TRUE)
  span <- diff(yr)
  if (!is.finite(span) || span == 0) span <- 1

  base_y <- yr[2] + 0.25 * span
  step_y <- 0.12 * span

  x_lvls <- levels(dd$Complex)
  x_idx <- seq_along(x_lvls)

  df_lines <- list()
  df_text <- list()
  row_id <- 1

  for (i in seq_along(x_lvls)) {

    cx <- x_lvls[i]
    di <- dd[dd$Complex == cx, , drop = FALSE]

    roles_here <- roles_keep[roles_keep %in% unique(as.character(di$role))]
    if (length(roles_here) < 2) next

    prs <- combn(roles_here, 2, simplify = FALSE)

    pvals <- sapply(prs, function(pr){
      x <- di[di$role == pr[1], yvar, drop = TRUE]
      y <- di[di$role == pr[2], yvar, drop = TRUE]
      if (length(x) < 2 || length(y) < 2) return(NA_real_)
      tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    })

    padj <- p.adjust(pvals, method = "BH")
    keep <- if (signif_only) which(!is.na(padj) & padj < alpha) else seq_along(padj)

    if (!length(keep)) next

    for (j in seq_along(keep)) {

      k <- keep[j]
      pr <- prs[[k]]

      y0 <- base_y + (j - 1) * step_y
      x1 <- dodge_x(x_idx[i], pr[1], roles_keep, dodge_w)
      x2 <- dodge_x(x_idx[i], pr[2], roles_keep, dodge_w)

      df_lines[[row_id]] <- data.frame(x = x1, xend = x2, y = y0, yend = y0)
      df_text[[row_id]]  <- data.frame(
        x = (x1 + x2) / 2,
        y = y0 + 0.02 * span,
        lab = p_to_stars(padj[k])
      )

      row_id <- row_id + 1
    }
  }

  if (length(df_lines)) {

    df_lines <- do.call(rbind, df_lines)
    df_text <- do.call(rbind, df_text)

    p <- p +
      geom_segment(
        data = df_lines,
        aes(x = x, xend = xend, y = y, yend = yend),
        inherit.aes = FALSE,
        linewidth = 0.6
      ) +
      geom_text(
        data = df_text,
        aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        size = 4.4,
        vjust = 0
      ) +
      coord_cartesian(
        ylim = c(yr[1] - 0.12 * span, max(df_text$y) + 0.12 * span),
        clip = "off"
      )

  } else {

    p <- p +
      coord_cartesian(
        ylim = c(yr[1] - 0.12 * span, yr[2] + 0.35 * span),
        clip = "off"
      )
  }

  return(save_plot(p, file.path(out_dir, outfile), width = 7.6, height = 4.8))
}

pC1 <- plot_complex("piN_mean", expression(pi[N]), "C_complex_mt_nu_ass_piN")
pC2 <- plot_complex("piS_mean", expression(pi[S]), "C_complex_mt_nu_ass_piS")
pC3 <- plot_complex("piN_piS_mean", expression(pi[N]/pi[S]), "C_complex_mt_nu_ass_piNpiS")

## =========================================================
## D: nuOXPHOS core vs noncore
## =========================================================

df_core <- nu %>%
  filter(role == "subunit") %>%
  mutate(
    GENE_UP = toupper(gene),
    COMPLEX_UP = toupper(gsub("\\s+", "", complex))
  )

has_suffix <- function(x, pat){
  grepl(paste0("^(", pat, ")[A-Z]?$"), x)
}

is_core_CI <- with(df_core, COMPLEX_UP == "CI" &
  (
    has_suffix(GENE_UP, "NDUFS1|NDUFS2|NDUFS3|NDUFS7|NDUFS8") |
      has_suffix(GENE_UP, "NDUFV1|NDUFV2")
  ))

is_core_CII <- with(df_core, COMPLEX_UP == "CII" &
  has_suffix(GENE_UP, "SDHA|SDHB|SDHC|SDHD"))

core_CV_F1 <- with(df_core, COMPLEX_UP == "CV" &
  (
    has_suffix(GENE_UP, "ATP5A1|ATP5B|ATP5C1|ATP5D") |
      has_suffix(GENE_UP, "ATP5F1A|ATP5F1B|ATP5F1C|ATP5F1D")
  ))

core_CV_c <- with(df_core, COMPLEX_UP == "CV" &
  has_suffix(GENE_UP, "ATP5G1|ATP5G2|ATP5G3|ATP5MC1|ATP5MC2|ATP5MC3"))

df_core$grp <- ifelse(
  is_core_CI | is_core_CII | core_CV_F1 | core_CV_c,
  "nu_core",
  "nu_noncore"
)

df_core$grp <- factor(df_core$grp, levels = c("nu_core", "nu_noncore"))

cat("core/noncore counts:\n")
print(table(df_core$grp))

pD1 <- plot_box_two(df_core, "grp", "piN_mean", expression(pi[N]),
                    fill_core, label_core, c("nu_core", "nu_noncore"),
                    "D_core_vs_noncore_piN")

pD2 <- plot_box_two(df_core, "grp", "piS_mean", expression(pi[S]),
                    fill_core, label_core, c("nu_core", "nu_noncore"),
                    "D_core_vs_noncore_piS")

pD3 <- plot_box_two(df_core, "grp", "piN_piS_mean", expression(pi[N]/pi[S]),
                    fill_core, label_core, c("nu_core", "nu_noncore"),
                    "D_core_vs_noncore_piNpiS")

## =========================================================
## E/F: ribo and ARS
## =========================================================

to_grp <- function(role_chr){
  x <- tolower(trimws(role_chr))
  x <- gsub("_", "-", x)

  is_ars  <- grepl("ars", x)
  is_ribo <- grepl("ribo|\\brp\\b|rp-", x)

  is_cyto <- grepl("^cyto|cytoplas", x)
  is_nmt  <- grepl("^nmt|mito|nuclear|nu-", x)

  if (is_ars && is_cyto) return("cyto-ars")
  if (is_ars && is_nmt)  return("nmt-ars")
  if (is_ribo && is_cyto) return("cyto-ribo")
  if (is_ribo && is_nmt)  return("nmt-ribo")

  return(NA_character_)
}

df_extra <- nu %>%
  mutate(grp = vapply(role, to_grp, character(1))) %>%
  filter(!is.na(grp))

cat("ARS/ribo counts:\n")
print(table(df_extra$grp))

df_ribo <- df_extra %>%
  filter(grp %in% c("cyto-ribo", "nmt-ribo")) %>%
  mutate(grp = factor(grp, levels = c("cyto-ribo", "nmt-ribo")))

df_ars <- df_extra %>%
  filter(grp %in% c("cyto-ars", "nmt-ars")) %>%
  mutate(grp = factor(grp, levels = c("cyto-ars", "nmt-ars")))

pE1 <- plot_box_two(df_ribo, "grp", "piN_mean", expression(pi[N]),
                    COL_RIBO, label_ribo, c("cyto-ribo", "nmt-ribo"),
                    "E_ribo_cyto_vs_Nmt_piN")

pE2 <- plot_box_two(df_ribo, "grp", "piS_mean", expression(pi[S]),
                    COL_RIBO, label_ribo, c("cyto-ribo", "nmt-ribo"),
                    "E_ribo_cyto_vs_Nmt_piS")

pE3 <- plot_box_two(df_ribo, "grp", "piN_piS_mean", expression(pi[N]/pi[S]),
                    COL_RIBO, label_ribo, c("cyto-ribo", "nmt-ribo"),
                    "E_ribo_cyto_vs_Nmt_piNpiS")

pF1 <- plot_box_two(df_ars, "grp", "piN_mean", expression(pi[N]),
                    COL_ARS, label_ars, c("cyto-ars", "nmt-ars"),
                    "F_ARS_cyto_vs_Nmt_piN")

pF2 <- plot_box_two(df_ars, "grp", "piS_mean", expression(pi[S]),
                    COL_ARS, label_ars, c("cyto-ars", "nmt-ars"),
                    "F_ARS_cyto_vs_Nmt_piS")

pF3 <- plot_box_two(df_ars, "grp", "piN_piS_mean", expression(pi[N]/pi[S]),
                    COL_ARS, label_ars, c("cyto-ars", "nmt-ars"),
                    "F_ARS_cyto_vs_Nmt_piNpiS")

## =========================================================
## COMBINE BIG FIGURE
## one letter per category, three metrics per row
## =========================================================

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 20),
  plot.tag.position = c(0.01, 0.98)
)

rowA <- (pA1 + labs(tag = "A") + tag_theme) | pA2 | pA3
rowB <- (pB1 + labs(tag = "B") + tag_theme) | pB2 | pB3
rowC <- (pC1 + labs(tag = "C") + tag_theme) | pC2 | pC3
rowD <- (pD1 + labs(tag = "D") + tag_theme) | pD2 | pD3
rowE <- (pE1 + labs(tag = "E") + tag_theme) | pE2 | pE3
rowF <- (pF1 + labs(tag = "F") + tag_theme) | pF2 | pF3

big_fig <- rowA / rowB / rowC / rowD / rowE / rowF +
  plot_layout(heights = c(1, 1, 1.15, 1, 1, 1))

print(big_fig)

ggsave(
  file.path(out_dir, "PNPS_all_panels_combined_big.pdf"),
  big_fig,
  width = 18,
  height = 26
)

ggsave(
  file.path(out_dir, "PNPS_all_panels_combined_big.png"),
  big_fig,
  width = 18,
  height = 26,
  dpi = 300
)

cat("\nDone. Big combined figure written to:\n")
cat(file.path(out_dir, "PNPS_all_panels_combined_big.pdf"), "\n")
cat(file.path(out_dir, "PNPS_all_panels_combined_big.png"), "\n")

