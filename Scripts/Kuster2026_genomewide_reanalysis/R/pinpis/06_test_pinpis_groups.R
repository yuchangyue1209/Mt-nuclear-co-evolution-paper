#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

root <- "/work/cyu/Kuster2026_genomewide_reanalysis/genomewide_pinpis"

input_file <- file.path(
    root,
    "results",
    "genomewide_nuclear_mt13_pinpis_combined.tsv"
)

test_output <- file.path(
    root,
    "results",
    "pinpis_group_tests.tsv"
)

summary_output <- file.path(
    root,
    "results",
    "pinpis_group_summaries.tsv"
)

x <- fread(
    input_file,
    na.strings=c("NA","NaN","")
)

metrics <- c(
    piN="piN_mean",
    piS="piS_mean",
    piN_piS="piN_piS_ratio_of_means"
)

results <- list()
summaries <- list()

add_summary <- function(
    data,
    panel,
    group_variable
) {
    for (metric_name in names(metrics)) {
        metric <- metrics[[metric_name]]

        z <- data[
            is.finite(get(metric)) &
            !is.na(get(group_variable))
        ]

        if (nrow(z)==0L)
            next

        s <- z[
            ,
            .(
                n_genes=.N,
                median=median(
                    get(metric)
                ),
                mean=mean(
                    get(metric)
                ),
                q25=quantile(
                    get(metric),
                    0.25
                ),
                q75=quantile(
                    get(metric),
                    0.75
                )
            ),
            by=group_variable
        ]

        setnames(
            s,
            group_variable,
            "group"
        )

        s[
            ,
            `:=`(
                panel=panel,
                metric=metric_name
            )
        ]

        summaries[[
            length(summaries)+1L
        ]] <<- s
    }
}

run_comparison <- function(
    data,
    panel,
    group_variable,
    group1,
    group2,
    metric_name,
    comparison_family=NULL
) {
    metric <- metrics[[metric_name]]

    a <- data[
        get(group_variable)==group1,
        get(metric)
    ]

    b <- data[
        get(group_variable)==group2,
        get(metric)
    ]

    a <- a[is.finite(a)]
    b <- b[is.finite(b)]

    if (length(a)<2L || length(b)<2L)
        return(NULL)

    test <- wilcox.test(
        a,
        b,
        exact=FALSE,
        correct=FALSE
    )

    mann_whitney_u <- unname(
        test$statistic
    )

    rank_biserial <- (
        2 * mann_whitney_u /
        (length(a)*length(b))
    ) - 1

    data.table(
        panel=panel,
        comparison_family=ifelse(
            is.null(comparison_family),
            panel,
            comparison_family
        ),
        metric=metric_name,
        group1=group1,
        group2=group2,
        n1=length(a),
        n2=length(b),
        median1=median(a),
        median2=median(b),
        mean1=mean(a),
        mean2=mean(b),
        rank_biserial=rank_biserial,
        p=unname(test$p.value)
    )
}

add_pair <- function(
    data,
    panel,
    group_variable,
    group1,
    group2,
    comparison_family=NULL
) {
    for (metric_name in names(metrics)) {
        result <- run_comparison(
            data=data,
            panel=panel,
            group_variable=group_variable,
            group1=group1,
            group2=group2,
            metric_name=metric_name,
            comparison_family=comparison_family
        )

        if (!is.null(result)) {
            results[[
                length(results)+1L
            ]] <<- result
        }
    }
}

# Panel A: mtOXPHOS, nuOXPHOS, assembly factor
panel_a <- x[
    plot_role %in% c(
        "mtOXPHOS",
        "nuOXPHOS",
        "assembly factor"
    )
]

add_summary(
    panel_a,
    "A_OXPHOS_roles",
    "plot_role"
)

add_pair(
    panel_a,
    "A_OXPHOS_roles",
    "plot_role",
    "mtOXPHOS",
    "nuOXPHOS"
)

add_pair(
    panel_a,
    "A_OXPHOS_roles",
    "plot_role",
    "mtOXPHOS",
    "assembly factor"
)

add_pair(
    panel_a,
    "A_OXPHOS_roles",
    "plot_role",
    "nuOXPHOS",
    "assembly factor"
)

# Panel B: within-complex comparisons
panel_b <- x[
    plot_role %in% c(
        "mtOXPHOS",
        "nuOXPHOS",
        "assembly factor"
    ) &
    !is.na(plot_complex)
]

add_summary(
    panel_b,
    "B_OXPHOS_complex",
    "plot_role"
)

for (complex_name in c(
    "CI",
    "CII",
    "CIII",
    "CIV",
    "CV"
)) {
    z <- panel_b[
        plot_complex==complex_name
    ]

    present <- unique(
        as.character(z$plot_role)
    )

    candidate_pairs <- list(
        c("mtOXPHOS","nuOXPHOS"),
        c("mtOXPHOS","assembly factor"),
        c("nuOXPHOS","assembly factor")
    )

    for (pair in candidate_pairs) {
        if (all(pair %in% present)) {
            add_pair(
                z,
                paste0(
                    "B_",
                    complex_name
                ),
                "plot_role",
                pair[1],
                pair[2],
                comparison_family=paste0(
                    "B_",
                    complex_name
                )
            )
        }
    }
}

# Panel C: ribosomal proteins
panel_c <- x[
    plot_role %in% c(
        "cyto-RP",
        "Nmt-RP"
    )
]

add_summary(
    panel_c,
    "C_RP",
    "plot_role"
)

add_pair(
    panel_c,
    "C_RP",
    "plot_role",
    "Nmt-RP",
    "cyto-RP"
)

# Panel D: aminoacyl-tRNA synthetases
panel_d <- x[
    plot_role %in% c(
        "cyto-ARS",
        "Nmt-ARS"
    )
]

add_summary(
    panel_d,
    "D_ARS",
    "plot_role"
)

add_pair(
    panel_d,
    "D_ARS",
    "plot_role",
    "Nmt-ARS",
    "cyto-ARS"
)

# Panel E: mtOXPHOS, nuclear core and noncore
panel_e <- x[
    !is.na(core_plot_group)
]

add_summary(
    panel_e,
    "E_core_status",
    "core_plot_group"
)

add_pair(
    panel_e,
    "E_core_status",
    "core_plot_group",
    "mtOXPHOS",
    "nu core"
)

add_pair(
    panel_e,
    "E_core_status",
    "core_plot_group",
    "mtOXPHOS",
    "nu noncore"
)

add_pair(
    panel_e,
    "E_core_status",
    "core_plot_group",
    "nu core",
    "nu noncore"
)

# Panel F: Kuster genome-wide classes
panel_f <- x[
    gene_source=="nuclear" &
    Kuster_class %in% c(
        "direct_n-mt",
        "indirect_n-mt",
        "non-n-mt"
    )
]

add_summary(
    panel_f,
    "F_Kuster",
    "Kuster_class"
)

add_pair(
    panel_f,
    "F_Kuster",
    "Kuster_class",
    "direct_n-mt",
    "indirect_n-mt"
)

add_pair(
    panel_f,
    "F_Kuster",
    "Kuster_class",
    "direct_n-mt",
    "non-n-mt"
)

add_pair(
    panel_f,
    "F_Kuster",
    "Kuster_class",
    "indirect_n-mt",
    "non-n-mt"
)

tests <- rbindlist(
    results,
    use.names=TRUE,
    fill=TRUE
)

summary_table <- rbindlist(
    summaries,
    use.names=TRUE,
    fill=TRUE
)

# Correct separately within each panel and metric.
tests[
    ,
    p_bh := p.adjust(
        p,
        method="BH"
    ),
    by=.(comparison_family,metric)
]

tests[
    ,
    significance := fifelse(
        p_bh < 0.001,
        "***",
        fifelse(
            p_bh < 0.01,
            "**",
            fifelse(
                p_bh < 0.05,
                "*",
                "ns"
            )
        )
    )
]

setorder(
    tests,
    panel,
    metric,
    p_bh
)

setcolorder(
    tests,
    c(
        "panel",
        "comparison_family",
        "metric",
        "group1",
        "group2",
        "n1",
        "n2",
        "median1",
        "median2",
        "mean1",
        "mean2",
        "rank_biserial",
        "p",
        "p_bh",
        "significance"
    )
)

fwrite(
    tests,
    test_output,
    sep="\t",
    quote=FALSE,
    na="NA"
)

fwrite(
    summary_table,
    summary_output,
    sep="\t",
    quote=FALSE,
    na="NA"
)

cat("===== Significant comparisons after BH =====\n")

print(
    tests[
        p_bh < 0.05
    ]
)

cat("\n===== Kuster comparisons =====\n")

print(
    tests[
        grepl(
            "^F_",
            panel
        )
    ]
)

cat("\n===== RP comparisons =====\n")

print(
    tests[
        panel=="C_RP"
    ]
)

cat("\n===== ARS comparisons =====\n")

print(
    tests[
        panel=="D_ARS"
    ]
)

cat("\n===== Core comparisons =====\n")

print(
    tests[
        panel=="E_core_status"
    ]
)

cat("\n[tests]",test_output,"\n")
cat("[summaries]",summary_output,"\n")
