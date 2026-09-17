#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

## ============================================================
## 1. Paths
## ============================================================

nuclear_file <- paste0(
    "/mnt/spareHD_2/genomewide_codeml_kuster/",
    "07_codeml_genomewide/",
    "codeml_master_analysis.tsv"
)

old_mt_file <- paste0(
    "/mnt/spareHD_2/oxphos_codeml_ready/",
    "09_codeml_sites_models/",
    "codeml_sites_summary.merged.tsv"
)

outdir <- paste0(
    "/mnt/spareHD_2/genomewide_codeml_kuster/",
    "09_figures/Figure2_complete_updated_large_fonts"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

statistics_file <- file.path(
    outdir,
    "Figure2_complete_statistics.tsv"
)

## ============================================================
## 2. Figure settings
## ============================================================

BORDER_LWD <- 0.6
BOX_LWD <- 0.8
TICK_LWD <- 0.6

AXIS_TITLE_SIZE <- 30
AXIS_TEXT_SIZE <- 21
STAR_SIZE <- 8
ROW_TAG_SIZE <- 32

POINT_SIZE <- 1.8
POINT_ALPHA <- 0.65

ALPHA_LEVEL <- 0.05
SIGNIF_ONLY <- TRUE

# Only affects displayed dots in the very large Kuster panel.
# Boxplots and tests always use all genes.
MAX_DISPLAY_POINTS_PER_GROUP <- 1500

set.seed(20260810)

## ============================================================
## 3. Morandi palettes
## ============================================================

palette_overall <- c(
    "mt"  = "#7FA4B2",
    "nu"  = "#BC8984",
    "ass" = "#9EAD8D"
)

labels_overall <- c(
    "mt"  = "mtOXPHOS",
    "nu"  = "nuOXPHOS",
    "ass" = "nu Assembly Factor"
)

palette_rp <- c(
    "cyto-ribo" = "#7795A8",
    "Nmt-ribo"  = "#C28C73"
)

labels_rp <- c(
    "cyto-ribo" = "cyto-RP",
    "Nmt-ribo"  = "Nmt-RP"
)

palette_ars <- c(
    "cyto-ARS" = "#7589A9",
    "Nmt-ARS"  = "#B77D9A"
)

labels_ars <- c(
    "cyto-ARS" = "cyto-ARS",
    "Nmt-ARS"  = "Nmt-ARS"
)

palette_core <- c(
    "mt"         = "#7FA4B2",
    "nu_core"    = "#BC8984",
    "nu_noncore" = "#B6A174"
)

labels_core <- c(
    "mt"         = "mtOXPHOS",
    "nu_core"    = "nu-core",
    "nu_noncore" = "nu-noncore"
)

palette_kuster <- c(
    "direct_n-mt"   = "#BC8984",
    "indirect_n-mt" = "#948C9D",
    "non-n-mt"      = "#AAA8A1"
)

labels_kuster <- c(
    "direct_n-mt"   = "direct n-mt",
    "indirect_n-mt" = "indirect n-mt",
    "non-n-mt"      = "non-n-mt"
)

## ============================================================
## 4. Read new nuclear data
## ============================================================

nu <- read.delim(
    nuclear_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

nu$dN <- suppressWarnings(
    as.numeric(nu$tree_length_dN)
)

nu$dS <- suppressWarnings(
    as.numeric(nu$tree_length_dS)
)

nu$omega_plot <- suppressWarnings(
    as.numeric(nu$omega_ES1)
)

if (nrow(nu) != 17965) {
    stop(
        "Expected 17,965 eligible nuclear genes; found ",
        nrow(nu)
    )
}

cat("===== Updated nuclear input =====\n")
cat("Eligible genes:", nrow(nu), "\n")
cat(
    "Invariant genes:",
    sum(nu$rate_status=="dN0_dS0"),
    "\n"
)
cat(
    "Reliable omega:",
    sum(is.finite(nu$omega_plot)),
    "\n"
)

## ============================================================
## 5. Read old mtOXPHOS M0 data
## ============================================================

mt_raw <- read.delim(
    old_mt_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

mt_raw$role <- tolower(
    trimws(mt_raw$role)
)

mt <- mt_raw[
    mt_raw$model == "M0" &
    mt_raw$role == "mt",
    ,
    drop = FALSE
]

mt$dN <- suppressWarnings(
    as.numeric(mt$dN)
)

mt$dS <- suppressWarnings(
    as.numeric(mt$dS)
)

mt$omega_original <- suppressWarnings(
    as.numeric(mt$omega)
)

# Preserve old mt filtering rule for identifiable omega.
mt$omega_plot <- ifelse(
    is.finite(mt$omega_original) &
    mt$omega_original < 999 &
    is.finite(mt$dS) &
    mt$dS >= 0.001,
    mt$omega_original,
    NA_real_
)

cat("\n===== Old mt input =====\n")
cat("mt genes:", nrow(mt), "\n")
cat(
    "mt genes with plotted omega:",
    sum(is.finite(mt$omega_plot)),
    "\n"
)

## ============================================================
## 6. Genome-wide dotted-line backgrounds
## ============================================================

# dN and dS use all 17,965 genes, including invariant zeros.
# omega uses only identifiable omega_ES1.
genome_medians <- c(
    dN = median(
        nu$dN,
        na.rm = TRUE
    ),
    dS = median(
        nu$dS,
        na.rm = TRUE
    ),
    omega_plot = median(
        nu$omega_plot,
        na.rm = TRUE
    )
)

cat("\n===== Genome-wide dotted lines =====\n")
print(genome_medians)

## ============================================================
## 7. Complex normalization
## ============================================================

normalize_complex <- function(value) {
    value <- toupper(
        gsub(
            "[[:space:]]+",
            "",
            trimws(as.character(value))
        )
    )

    output <- rep(
        NA_character_,
        length(value)
    )

    output[value %in% c("CI", "I")] <- "I"
    output[value %in% c("CII", "II")] <- "II"
    output[value %in% c("CIII", "III")] <- "III"
    output[value %in% c("CIV", "IV")] <- "IV"
    output[value %in% c("CV", "V")] <- "V"

    output
}

mt$Complex <- normalize_complex(
    mt$complex
)

nu$Complex <- normalize_complex(
    nu$own_complex
)

## ============================================================
## 8. Prepare datasets for rows A-F
## ============================================================

common_columns <- c(
    "group",
    "Complex",
    "dN",
    "dS",
    "omega_plot"
)

# ---------- Row A: mt / nu subunit / assembly ----------

A_mt <- data.frame(
    group = "mt",
    Complex = mt$Complex,
    dN = mt$dN,
    dS = mt$dS,
    omega_plot = mt$omega_plot
)

A_nu <- data.frame(
    group = ifelse(
        nu$own_role == "subunit",
        "nu",
        ifelse(
            nu$own_role == "assembly_factor",
            "ass",
            NA_character_
        )
    ),
    Complex = nu$Complex,
    dN = nu$dN,
    dS = nu$dS,
    omega_plot = nu$omega_plot
)

A_data <- rbind(
    A_mt[, common_columns],
    A_nu[
        !is.na(A_nu$group),
        common_columns
    ]
)

# ---------- Row B: complexes ----------

B_data <- A_data[
    !is.na(A_data$Complex),
    ,
    drop = FALSE
]

B_data$Complex <- factor(
    B_data$Complex,
    levels = c("I", "II", "III", "IV", "V")
)

# ---------- Row C: RP ----------

C_data <- data.frame(
    group = nu$own_role,
    Complex = NA_character_,
    dN = nu$dN,
    dS = nu$dS,
    omega_plot = nu$omega_plot
)

C_data <- C_data[
    !is.na(C_data$group) &
    C_data$group %in% c(
        "cyto-ribo",
        "Nmt-ribo"
    ),
    common_columns
]

# ---------- Row D: ARS ----------

D_data <- data.frame(
    group = nu$own_role,
    Complex = NA_character_,
    dN = nu$dN,
    dS = nu$dS,
    omega_plot = nu$omega_plot
)

D_data <- D_data[
    !is.na(D_data$group) &
    D_data$group %in% c(
        "cyto-ARS",
        "Nmt-ARS"
    ),
    common_columns
]

# ---------- Row E: mt / nu-core / nu-noncore ----------

E_mt <- data.frame(
    group = "mt",
    Complex = mt$Complex,
    dN = mt$dN,
    dS = mt$dS,
    omega_plot = mt$omega_plot
)

E_nu <- data.frame(
    group = ifelse(
        nu$own_role == "subunit" &
        nu$core_status == "nu_core",
        "nu_core",
        ifelse(
            nu$own_role == "subunit" &
            nu$core_status == "nu_noncore",
            "nu_noncore",
            NA_character_
        )
    ),
    Complex = nu$Complex,
    dN = nu$dN,
    dS = nu$dS,
    omega_plot = nu$omega_plot
)

E_data <- rbind(
    E_mt[, common_columns],
    E_nu[
        !is.na(E_nu$group),
        common_columns
    ]
)

# ---------- Row F: Kuster three classes ----------

F_data <- data.frame(
    group = nu$Kuster_class,
    Complex = NA_character_,
    dN = nu$dN,
    dS = nu$dS,
    omega_plot = nu$omega_plot
)

F_data <- F_data[
    F_data$group %in% c(
        "direct_n-mt",
        "indirect_n-mt",
        "non-n-mt"
    ),
    common_columns
]

## ============================================================
## 9. Statistical helpers
## ============================================================

p_to_stars <- function(p) {
    if (is.na(p)) {
        return("n.s.")
    }

    if (p < 0.001) {
        return("***")
    }

    if (p < 0.01) {
        return("**")
    }

    if (p < 0.05) {
        return("*")
    }

    "n.s."
}

sample_display_points <- function(
    data,
    max_per_group
) {
    pieces <- split(
        data,
        data$group
    )

    sampled <- lapply(
        pieces,
        function(piece) {
            if (nrow(piece) <= max_per_group) {
                return(piece)
            }

            piece[
                sample(
                    seq_len(nrow(piece)),
                    max_per_group
                ),
                ,
                drop = FALSE
            ]
        }
    )

    do.call(
        rbind,
        sampled
    )
}

calculate_simple_statistics <- function(
    data,
    metric,
    group_levels,
    row_name
) {
    pairs <- combn(
        group_levels,
        2,
        simplify = FALSE
    )

    result <- lapply(
        pairs,
        function(pair) {
            a <- data[[metric]][
                data$group == pair[1]
            ]

            b <- data[[metric]][
                data$group == pair[2]
            ]

            a <- a[is.finite(a)]
            b <- b[is.finite(b)]

            if (
                length(a) < 2 ||
                length(b) < 2
            ) {
                p <- NA_real_
                U <- NA_real_
                effect <- NA_real_
            } else {
                test <- suppressWarnings(
                    wilcox.test(
                        a,
                        b,
                        exact = FALSE
                    )
                )

                p <- test$p.value
                U <- unname(test$statistic)

                effect <- (
                    2 * U /
                    (length(a) * length(b))
                ) - 1
            }

            data.frame(
                row = row_name,
                metric = metric,
                Complex = NA_character_,
                group1 = pair[1],
                group2 = pair[2],
                n1 = length(a),
                n2 = length(b),
                median1 = median(a, na.rm=TRUE),
                median2 = median(b, na.rm=TRUE),
                mean1 = mean(a, na.rm=TRUE),
                mean2 = mean(b, na.rm=TRUE),
                rank_biserial = effect,
                p_raw = p,
                stringsAsFactors = FALSE
            )
        }
    )

    result <- do.call(
        rbind,
        result
    )

    result$p_adj <- p.adjust(
        result$p_raw,
        method = "BH"
    )

    result$stars <- vapply(
        result$p_adj,
        p_to_stars,
        character(1)
    )

    result
}

calculate_complex_statistics <- function(
    data,
    metric,
    group_levels
) {
    results <- list()
    index <- 1

    for (complex_name in levels(data$Complex)) {
        subset_data <- data[
            data$Complex == complex_name,
            ,
            drop = FALSE
        ]

        present <- group_levels[
            group_levels %in%
            unique(subset_data$group)
        ]

        if (length(present) < 2) {
            next
        }

        pairs <- combn(
            present,
            2,
            simplify = FALSE
        )

        complex_results <- lapply(
            pairs,
            function(pair) {
                a <- subset_data[[metric]][
                    subset_data$group == pair[1]
                ]

                b <- subset_data[[metric]][
                    subset_data$group == pair[2]
                ]

                a <- a[is.finite(a)]
                b <- b[is.finite(b)]

                if (
                    length(a) < 2 ||
                    length(b) < 2
                ) {
                    p <- NA_real_
                    U <- NA_real_
                    effect <- NA_real_
                } else {
                    test <- suppressWarnings(
                        wilcox.test(
                            a,
                            b,
                            exact = FALSE
                        )
                    )

                    p <- test$p.value
                    U <- unname(test$statistic)

                    effect <- (
                        2 * U /
                        (length(a) * length(b))
                    ) - 1
                }

                data.frame(
                    row = "B_complex",
                    metric = metric,
                    Complex = complex_name,
                    group1 = pair[1],
                    group2 = pair[2],
                    n1 = length(a),
                    n2 = length(b),
                    median1 = median(a, na.rm=TRUE),
                    median2 = median(b, na.rm=TRUE),
                    mean1 = mean(a, na.rm=TRUE),
                    mean2 = mean(b, na.rm=TRUE),
                    rank_biserial = effect,
                    p_raw = p,
                    stringsAsFactors = FALSE
                )
            }
        )

        complex_results <- do.call(
            rbind,
            complex_results
        )

        # BH correction within each complex and metric
        complex_results$p_adj <- p.adjust(
            complex_results$p_raw,
            method = "BH"
        )

        complex_results$stars <- vapply(
            complex_results$p_adj,
            p_to_stars,
            character(1)
        )

        results[[index]] <- complex_results
        index <- index + 1
    }

    do.call(
        rbind,
        results
    )
}

## ============================================================
## 10. Plot theme
## ============================================================

figure_theme <- theme_bw(base_size = AXIS_TEXT_SIZE) +
    theme(
        text = element_text(
            family = "sans",
            face = "bold",
            colour = "black"
        ),
        axis.title.y = element_text(
            size = AXIS_TITLE_SIZE,
            margin = margin(r = 8)
        ),
        axis.text.x = element_text(
            size = AXIS_TEXT_SIZE,
            colour = "black",
            angle = 45,
            hjust = 1
        ),
        axis.text.y = element_text(
            size = AXIS_TEXT_SIZE,
            colour = "black"
        ),
        axis.ticks = element_line(
            colour = "black",
            linewidth = TICK_LWD
        ),
        panel.border = element_rect(
            colour = "black",
            fill = NA,
            linewidth = BORDER_LWD
        ),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        aspect.ratio = 1,
        plot.margin = margin(
            8, 9, 8, 8
        )
    )

## ============================================================
## 11. Simple-row plotting function
## ============================================================

make_simple_plot <- function(
    data,
    metric,
    ylab,
    group_levels,
    colors,
    labels,
    row_name,
    point_cap = Inf
) {
    dd <- data[
        is.finite(data[[metric]]),
        ,
        drop = FALSE
    ]

    dd$group <- factor(
        dd$group,
        levels = group_levels
    )

    statistics <- calculate_simple_statistics(
        dd,
        metric,
        group_levels,
        row_name
    )

    display_statistics <- statistics

    if (SIGNIF_ONLY) {
        display_statistics <- display_statistics[
            !is.na(display_statistics$p_adj) &
            display_statistics$p_adj < ALPHA_LEVEL,
            ,
            drop = FALSE
        ]
    }

    point_data <- if (
        is.finite(point_cap)
    ) {
        sample_display_points(
            dd,
            point_cap
        )
    } else {
        dd
    }

    value_range <- range(
        dd[[metric]],
        na.rm = TRUE
    )

    span <- diff(value_range)

    if (
        !is.finite(span) ||
        span == 0
    ) {
        span <- 1
    }

    g <- ggplot(
        dd,
        aes(
            x = group,
            y = .data[[metric]],
            fill = group
        )
    ) +
        geom_hline(
            yintercept = genome_medians[[metric]],
            linetype = "dotted",
            linewidth = 0.9,
            colour = "black"
        ) +
        geom_boxplot(
            width = 0.52,
            colour = "black",
            linewidth = BOX_LWD,
            outlier.shape = NA
        ) +
        geom_point(
            data = point_data,
            position = position_jitter(
                width = 0.11,
                height = 0,
                seed = 20260810
            ),
            shape = 21,
            size = POINT_SIZE,
            stroke = 0.3,
            colour = "black",
            alpha = POINT_ALPHA
        ) +
        scale_fill_manual(
            values = colors,
            guide = "none"
        ) +
        scale_x_discrete(
            labels = labels,
            drop = FALSE
        ) +
        labs(
            x = NULL,
            y = ylab
        ) +
        figure_theme

    if (nrow(display_statistics) > 0) {
        for (
            index in seq_len(
                nrow(display_statistics)
            )
        ) {
            current <- display_statistics[index, ]

            x1 <- match(
                current$group1,
                group_levels
            )

            x2 <- match(
                current$group2,
                group_levels
            )

            y <- (
                value_range[2] +
                index * 0.13 * span
            )

            tip <- 0.025 * span

            g <- g +
                annotate(
                    "segment",
                    x = x1,
                    xend = x2,
                    y = y,
                    yend = y,
                    linewidth = 0.6
                ) +
                annotate(
                    "segment",
                    x = x1,
                    xend = x1,
                    y = y,
                    yend = y-tip,
                    linewidth = 0.6
                ) +
                annotate(
                    "segment",
                    x = x2,
                    xend = x2,
                    y = y,
                    yend = y-tip,
                    linewidth = 0.6
                ) +
                annotate(
                    "text",
                    x = (x1+x2)/2,
                    y = y,
                    label = current$stars,
                    size = STAR_SIZE,
                    fontface = "bold",
                    vjust = -0.3
                )
        }

        y_max <- (
            value_range[2] +
            (nrow(display_statistics)+1) *
            0.15 * span
        )

        g <- g +
            coord_cartesian(
                ylim = c(
                    min(0, value_range[1]),
                    y_max
                ),
                clip = "off"
            )
    }

    list(
        plot = g,
        statistics = statistics
    )
}

## ============================================================
## 12. Complex-row plotting function
## ============================================================

make_complex_plot <- function(
    data,
    metric,
    ylab
) {
    dd <- data[
        is.finite(data[[metric]]),
        ,
        drop = FALSE
    ]

    group_levels <- c(
        "mt",
        "nu",
        "ass"
    )

    dd$group <- factor(
        dd$group,
        levels = group_levels
    )

    statistics <- calculate_complex_statistics(
        dd,
        metric,
        group_levels
    )

    point_data <- dd

    dodge_width <- 0.76

    g <- ggplot(
        dd,
        aes(
            x = Complex,
            y = .data[[metric]],
            fill = group
        )
    ) +
        geom_hline(
            yintercept = genome_medians[[metric]],
            linetype = "dotted",
            linewidth = 0.9,
            colour = "black"
        ) +
        geom_boxplot(
            position = position_dodge(
                width = dodge_width
            ),
            width = 0.62,
            colour = "black",
            linewidth = BOX_LWD,
            outlier.shape = NA
        ) +
        geom_point(
            data = point_data,
            position = position_jitterdodge(
                jitter.width = 0.10,
                dodge.width = dodge_width,
                seed = 20260810
            ),
            shape = 21,
            size = 1.2,
            stroke = 0.25,
            colour = "black",
            alpha = POINT_ALPHA
        ) +
        scale_fill_manual(
            values = palette_overall,
            labels = c(
                mt = "mt",
                nu = "nu",
                ass = "ass"
            )
        ) +
        labs(
            x = NULL,
            y = ylab
        ) +
        figure_theme +
        theme(
            axis.text.x = element_text(
                size = AXIS_TEXT_SIZE,
                colour = "black",
                angle = 0,
                hjust = 0.5
            ),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(
                size = 20,
                face = "bold"
            )
        )

    value_range <- range(
        dd[[metric]],
        na.rm = TRUE
    )

    span <- diff(value_range)

    if (
        !is.finite(span) ||
        span == 0
    ) {
        span <- 1
    }

    offsets <- c(
        mt = -0.25,
        nu = 0,
        ass = 0.25
    )

    highest_annotation <- value_range[2]

    for (
        complex_name in levels(dd$Complex)
    ) {
        current_statistics <- statistics[
            statistics$Complex == complex_name,
            ,
            drop = FALSE
        ]

        if (SIGNIF_ONLY) {
            current_statistics <- current_statistics[
                !is.na(current_statistics$p_adj) &
                current_statistics$p_adj < ALPHA_LEVEL,
                ,
                drop = FALSE
            ]
        }

        if (nrow(current_statistics) == 0) {
            next
        }

        complex_index <- match(
            complex_name,
            levels(dd$Complex)
        )

        for (
            row_index in seq_len(
                nrow(current_statistics)
            )
        ) {
            current <- current_statistics[
                row_index,
            ]

            x1 <- (
                complex_index +
                offsets[[current$group1]]
            )

            x2 <- (
                complex_index +
                offsets[[current$group2]]
            )

            y <- (
                value_range[2] +
                row_index * 0.12 * span
            )

            tip <- 0.025 * span

            g <- g +
                annotate(
                    "segment",
                    x = x1,
                    xend = x2,
                    y = y,
                    yend = y,
                    linewidth = 0.55
                ) +
                annotate(
                    "segment",
                    x = x1,
                    xend = x1,
                    y = y,
                    yend = y-tip,
                    linewidth = 0.55
                ) +
                annotate(
                    "segment",
                    x = x2,
                    xend = x2,
                    y = y,
                    yend = y-tip,
                    linewidth = 0.55
                ) +
                annotate(
                    "text",
                    x = (x1+x2)/2,
                    y = y,
                    label = current$stars,
                    size = STAR_SIZE,
                    fontface = "bold",
                    vjust = -0.3
                )

            highest_annotation <- max(
                highest_annotation,
                y
            )
        }
    }

    if (
        highest_annotation >
        value_range[2]
    ) {
        g <- g +
            coord_cartesian(
                ylim = c(
                    min(0, value_range[1]),
                    highest_annotation +
                    0.12 * span
                ),
                clip = "off"
            )
    }

    list(
        plot = g,
        statistics = statistics
    )
}

## ============================================================
## 13. Build each row
## ============================================================

metric_names <- c(
    "dN",
    "dS",
    "omega_plot"
)

metric_labels <- list(
    dN = expression(d[N]),
    dS = expression(d[S]),
    omega_plot = expression(d[N]/d[S])
)

make_three_simple_panels <- function(
    data,
    group_levels,
    colors,
    labels,
    row_name,
    point_cap = Inf
) {
    results <- lapply(
        metric_names,
        function(metric) {
            make_simple_plot(
                data = data,
                metric = metric,
                ylab = metric_labels[[metric]],
                group_levels = group_levels,
                colors = colors,
                labels = labels,
                row_name = row_name,
                point_cap = point_cap
            )
        }
    )

    names(results) <- metric_names

    row_plot <- (
        results$dN$plot |
        results$dS$plot |
        results$omega_plot$plot
    )

    list(
        plot = row_plot,
        results = results
    )
}

row_A <- make_three_simple_panels(
    A_data,
    c("mt", "nu", "ass"),
    palette_overall,
    labels_overall,
    "A_overall"
)

B_results <- lapply(
    metric_names,
    function(metric) {
        make_complex_plot(
            B_data,
            metric,
            metric_labels[[metric]]
        )
    }
)

names(B_results) <- metric_names

row_B <- list(
    plot = (
        B_results$dN$plot |
        B_results$dS$plot |
        B_results$omega_plot$plot
    ),
    results = B_results
)

row_C <- make_three_simple_panels(
    C_data,
    c("cyto-ribo", "Nmt-ribo"),
    palette_rp,
    labels_rp,
    "C_RP"
)

row_D <- make_three_simple_panels(
    D_data,
    c("cyto-ARS", "Nmt-ARS"),
    palette_ars,
    labels_ars,
    "D_ARS"
)

row_E <- make_three_simple_panels(
    E_data,
    c("mt", "nu_core", "nu_noncore"),
    palette_core,
    labels_core,
    "E_core"
)

row_F <- make_three_simple_panels(
    F_data,
    c(
        "direct_n-mt",
        "indirect_n-mt",
        "non-n-mt"
    ),
    palette_kuster,
    labels_kuster,
    "F_Kuster",
    point_cap = MAX_DISPLAY_POINTS_PER_GROUP
)

## ============================================================
## 14. Add row labels A-F
## ============================================================

add_row_tag <- function(
    row_object,
    tag
) {
    row_object +
        plot_annotation(
            tag_levels = list(
                c(tag, "", "")
            )
        ) &
        theme(
            plot.tag = element_text(
                size = ROW_TAG_SIZE,
                face = "bold",
                colour = "black"
            )
        )
}

plot_A <- add_row_tag(
    row_A$plot,
    "A"
)

plot_B <- add_row_tag(
    row_B$plot,
    "B"
)

plot_C <- add_row_tag(
    row_C$plot,
    "C"
)

plot_D <- add_row_tag(
    row_D$plot,
    "D"
)

plot_E <- add_row_tag(
    row_E$plot,
    "E"
)

plot_F <- add_row_tag(
    row_F$plot,
    "F"
)

## ============================================================
## 15. Save each row separately
## ============================================================

save_row <- function(
    plot,
    filename_stub,
    height = 8.5
) {
    ggsave(
        paste0(filename_stub, ".png"),
        plot,
        width = 21,
        height = height,
        dpi = 300,
        bg = "white"
    )

    ggsave(
        paste0(filename_stub, ".pdf"),
        plot,
        width = 21,
        height = height,
        device = cairo_pdf,
        bg = "white"
    )
}

save_row(
    plot_A,
    file.path(
        outdir,
        "Figure2A_mt_nu_assembly"
    )
)

save_row(
    plot_B,
    file.path(
        outdir,
        "Figure2B_complex"
    ),
    height = 9
)

save_row(
    plot_C,
    file.path(
        outdir,
        "Figure2C_RP"
    )
)

save_row(
    plot_D,
    file.path(
        outdir,
        "Figure2D_ARS"
    )
)

save_row(
    plot_E,
    file.path(
        outdir,
        "Figure2E_core_noncore"
    )
)

save_row(
    plot_F,
    file.path(
        outdir,
        "Figure2F_Kuster"
    )
)

## ============================================================
## 16. Combine large Figure 2
## ============================================================

combined <- wrap_plots(
    plot_A,
    plot_B,
    plot_C,
    plot_D,
    plot_E,
    plot_F,
    ncol = 1,
    heights = c(
        1,
        1.08,
        1,
        1,
        1,
        1
    )
) +
    plot_annotation(
        caption = paste0(
            "Dotted lines indicate genome-wide medians. ",
            "All 17,965 quality-eligible genes, including ",
            "invariant genes with dN=dS=0, contributed to ",
            "the dN and dS backgrounds. Omega was retained ",
            "only when identifiable and supported by at least ",
            "one expected synonymous change."
        )
    ) &
    theme(
        plot.caption = element_text(
            size = 15,
            face = "plain",
            hjust = 0
        )
    )

ggsave(
    file.path(
        outdir,
        "Figure2_combined_A_to_F.png"
    ),
    combined,
    width = 21,
    height = 52,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
)

ggsave(
    file.path(
        outdir,
        "Figure2_combined_A_to_F.pdf"
    ),
    combined,
    width = 21,
    height = 52,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
)

## ============================================================
## 17. Combine statistics
## ============================================================

collect_statistics <- function(
    row_object
) {
    do.call(
        rbind,
        lapply(
            row_object$results,
            function(result) {
                result$statistics
            }
        )
    )
}

all_statistics <- rbind(
    collect_statistics(row_A),
    collect_statistics(row_B),
    collect_statistics(row_C),
    collect_statistics(row_D),
    collect_statistics(row_E),
    collect_statistics(row_F)
)

write.table(
    all_statistics,
    statistics_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

## ============================================================
## 18. Final report
## ============================================================

cat("\n===== Figure 2 completed =====\n")
cat("Output directory:", outdir, "\n")
cat("Statistics:", statistics_file, "\n")

cat("\nIndividual PDFs:\n")
cat("Figure2A_mt_nu_assembly.pdf\n")
cat("Figure2B_complex.pdf\n")
cat("Figure2C_RP.pdf\n")
cat("Figure2D_ARS.pdf\n")
cat("Figure2E_core_noncore.pdf\n")
cat("Figure2F_Kuster.pdf\n")

cat("\nCombined PDF:\n")
cat("Figure2_combined_A_to_F.pdf\n")

cat("\nGenome-wide dotted lines:\n")
print(genome_medians)
