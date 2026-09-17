#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

root <- "/work/cyu/Kuster2026_genomewide_reanalysis/genomewide_pinpis"

nuclear_file <- file.path(
    root,
    "results",
    "genomewide_20347_pinpis_gene_level.tsv"
)

mt_file <- file.path(
    root,
    "results",
    "mt13_pinpis_population.corrected.tsv"
)

output_file <- file.path(
    root,
    "results",
    "genomewide_nuclear_mt13_pinpis_combined.tsv"
)

message("[read] Nuclear gene-level table")

nu <- fread(
    nuclear_file,
    na.strings=c("NA","NaN","")
)

nu <- nu[
    pinpis_main_eligible==TRUE
]

nu[
    ,
    gene_source := "nuclear"
]

nu[
    ,
    plot_role := fifelse(
        own_role=="subunit",
        "nuOXPHOS",
        fifelse(
            own_role=="assembly_factor",
            "assembly factor",
            fifelse(
                own_role=="Nmt-ribo",
                "Nmt-RP",
                fifelse(
                    own_role=="cyto-ribo",
                    "cyto-RP",
                    fifelse(
                        own_role=="Nmt-ARS",
                        "Nmt-ARS",
                        fifelse(
                            own_role=="cyto-ARS",
                            "cyto-ARS",
                            NA_character_
                        )
                    )
                )
            )
        )
    )
]

nu[
    ,
    plot_complex := fifelse(
        own_complex %in% c("CI","CII","CIII","CIV","CV"),
        own_complex,
        NA_character_
    )
]

nu[
    ,
    core_plot_group := fifelse(
        own_role=="subunit" &
        core_status=="nu_core",
        "nu core",
        fifelse(
            own_role=="subunit" &
            core_status=="nu_noncore",
            "nu noncore",
            NA_character_
        )
    )
]

message("[read] mt population-level table")

mt_population <- fread(
    mt_file,
    na.strings=c("NA","NaN","")
)

mt <- mt_population[
    ,
    {
        mean_piN <- mean(
            piN,
            na.rm=TRUE
        )

        mean_piS <- mean(
            piS,
            na.rm=TRUE
        )

        ratio <- if (
            is.finite(mean_piN) &&
            is.finite(mean_piS) &&
            mean_piS > 0
        ) {
            mean_piN/mean_piS
        } else {
            NA_real_
        }

        list(
            symbol=first(symbol),
            Kuster_class="mtDNA",
            eligible_populations=sum(
                is.finite(piN) &
                is.finite(piS)
            ),
            piN_mean=mean_piN,
            piS_mean=mean_piS,
            piN_piS_ratio_of_means=ratio,
            chromosome="mtDNA",
            autosome_only=FALSE,
            dnds_primary_gene=FALSE,
            old_name=first(symbol),
            new_symbol=first(symbol),
            own_role="mt_subunit",
            own_complex=first(complex),
            core_status="mt_core",
            old_target_gene=TRUE,
            gene_source="mitochondrial",
            plot_role="mtOXPHOS",
            plot_complex=first(complex),
            core_plot_group="mtOXPHOS"
        )
    },
    by=gene_id
]

required_columns <- c(
    "gene_id",
    "symbol",
    "Kuster_class",
    "eligible_populations",
    "piN_mean",
    "piS_mean",
    "piN_piS_ratio_of_means",
    "chromosome",
    "autosome_only",
    "dnds_primary_gene",
    "old_name",
    "new_symbol",
    "own_role",
    "own_complex",
    "core_status",
    "old_target_gene",
    "gene_source",
    "plot_role",
    "plot_complex",
    "core_plot_group"
)

for (column in setdiff(
    required_columns,
    names(nu)
)) {
    nu[
        ,
        (column) := NA
    ]
}

combined <- rbindlist(
    list(
        nu[
            ,
            ..required_columns
        ],
        mt[
            ,
            ..required_columns
        ]
    ),
    use.names=TRUE,
    fill=TRUE
)

combined[
    ,
    plot_complex := factor(
        plot_complex,
        levels=c(
            "CI",
            "CII",
            "CIII",
            "CIV",
            "CV"
        )
    )
]

combined[
    ,
    plot_role := factor(
        plot_role,
        levels=c(
            "mtOXPHOS",
            "nuOXPHOS",
            "assembly factor",
            "cyto-RP",
            "Nmt-RP",
            "cyto-ARS",
            "Nmt-ARS"
        )
    )
]

combined[
    ,
    core_plot_group := factor(
        core_plot_group,
        levels=c(
            "mtOXPHOS",
            "nu core",
            "nu noncore"
        )
    )
]

setorder(
    combined,
    gene_source,
    gene_id
)

fwrite(
    combined,
    output_file,
    sep="\t",
    quote=FALSE,
    na="NA"
)

cat("===== Combined dimensions =====\n")
cat("Rows:",nrow(combined),"\n")
cat("Nuclear:",sum(combined$gene_source=="nuclear"),"\n")
cat("Mitochondrial:",sum(combined$gene_source=="mitochondrial"),"\n")

cat("\n===== Plot roles =====\n")
print(
    combined[
        !is.na(plot_role),
        .N,
        by=plot_role
    ]
)

cat("\n===== OXPHOS complex x role =====\n")
print(
    combined[
        !is.na(plot_complex) &
        plot_role %in% c(
            "mtOXPHOS",
            "nuOXPHOS",
            "assembly factor"
        ),
        .N,
        by=.(plot_complex,plot_role)
    ][
        order(plot_complex,plot_role)
    ]
)

cat("\n===== Core groups =====\n")
print(
    combined[
        !is.na(core_plot_group),
        .N,
        by=core_plot_group
    ]
)

cat("\n===== Genome-wide medians =====\n")
print(
    combined[
        gene_source=="nuclear",
        .(
            median_piN=median(
                piN_mean,
                na.rm=TRUE
            ),
            median_piS=median(
                piS_mean,
                na.rm=TRUE
            ),
            median_piN_piS=median(
                piN_piS_ratio_of_means,
                na.rm=TRUE
            )
        )
    ]
)

cat("\n===== mtOXPHOS medians =====\n")
print(
    combined[
        gene_source=="mitochondrial",
        .(
            median_piN=median(
                piN_mean,
                na.rm=TRUE
            ),
            median_piS=median(
                piS_mean,
                na.rm=TRUE
            ),
            median_piN_piS=median(
                piN_piS_ratio_of_means,
                na.rm=TRUE
            )
        )
    ]
)

cat("\n[output]",output_file,"\n")
