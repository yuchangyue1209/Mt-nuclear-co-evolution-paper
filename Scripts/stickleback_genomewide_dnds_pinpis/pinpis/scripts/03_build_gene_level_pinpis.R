#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

root <- "/work/cyu/Kuster2026_genomewide_reanalysis/genomewide_pinpis"

population_file <- file.path(
    root,
    "results",
    "genomewide_20347_pinpis_population.tsv"
)

output_file <- file.path(
    root,
    "results",
    "genomewide_20347_pinpis_gene_level.tsv"
)

chromosome_file <- paste0(
    "/work/cyu/Kuster2026_genomewide_reanalysis/",
    "genomewide_pbs/stickleback_20426_gene_spans.bed"
)

dnds_master_file <- paste0(
    "/mnt/spareHD_2/genomewide_codeml_kuster/",
    "07_codeml_genomewide/codeml_master_analysis.tsv"
)

old_label_file <- paste0(
    "/work/cyu/Kuster2026_genomewide_reanalysis/",
    "old287_vs_genomewide/old287_three_labels_codeml20347.tsv"
)

minimum_eligible_populations <- 22L
total_populations <- 27L

message("[read] Population-level piN/piS")

x <- fread(
    population_file,
    na.strings=c("NA","NaN","")
)

numeric_columns <- c(
    "piN",
    "piS",
    "piN_piS",
    "piN_sum",
    "piS_sum",
    "N_opportunities",
    "S_opportunities",
    "callable_sites",
    "CDS_length",
    "callable_fraction"
)

for (column in numeric_columns) {
    set(
        x,
        j=column,
        value=as.numeric(x[[column]])
    )
}

if (nrow(x) != 20347L * 27L) {
    stop(
        "Unexpected population-table dimensions: ",
        nrow(x)
    )
}

if (uniqueN(x$gene_id) != 20347L) {
    stop("Unexpected number of genes")
}

if (uniqueN(x$population) != 27L) {
    stop("Unexpected number of populations")
}

if (any(x$population=="Norway")) {
    stop("Norway unexpectedly present")
}

x[
    ,
    row_eligible :=
        pinpis_eligible=="yes" &
        is.finite(piN) &
        is.finite(piS)
]

message("[summarize] Gene-level estimates")

gene <- x[
    ,
    {
        keep <- row_eligible

        n_eligible <- sum(keep)
        n_finite_population_ratios <- sum(
            keep &
            is.finite(piN_piS)
        )

        if (n_eligible > 0L) {
            mean_piN <- mean(piN[keep])
            mean_piS <- mean(piS[keep])

            median_piN <- median(piN[keep])
            median_piS <- median(piS[keep])

            pooled_piN <- sum(
                piN_sum[keep],
                na.rm=TRUE
            ) / sum(
                N_opportunities[keep],
                na.rm=TRUE
            )

            pooled_piS <- sum(
                piS_sum[keep],
                na.rm=TRUE
            ) / sum(
                S_opportunities[keep],
                na.rm=TRUE
            )

            mean_callable_fraction <- mean(
                callable_fraction[keep]
            )

            minimum_callable_fraction <- min(
                callable_fraction[keep]
            )
        } else {
            mean_piN <- NA_real_
            mean_piS <- NA_real_

            median_piN <- NA_real_
            median_piS <- NA_real_

            pooled_piN <- NA_real_
            pooled_piS <- NA_real_

            mean_callable_fraction <- NA_real_
            minimum_callable_fraction <- NA_real_
        }

        if (
            is.finite(mean_piN) &&
            is.finite(mean_piS) &&
            mean_piS > 0
        ) {
            mean_ratio <- mean_piN / mean_piS
            ratio_status <- "finite"
        } else if (
            is.finite(mean_piN) &&
            is.finite(mean_piS) &&
            mean_piN==0 &&
            mean_piS==0
        ) {
            mean_ratio <- NA_real_
            ratio_status <- "piN0_piS0"
        } else if (
            is.finite(mean_piN) &&
            is.finite(mean_piS) &&
            mean_piN > 0 &&
            mean_piS==0
        ) {
            mean_ratio <- NA_real_
            ratio_status <- "piNpositive_piS0"
        } else {
            mean_ratio <- NA_real_
            ratio_status <- "insufficient_data"
        }

        if (
            is.finite(pooled_piN) &&
            is.finite(pooled_piS) &&
            pooled_piS > 0
        ) {
            pooled_ratio <- pooled_piN / pooled_piS
        } else {
            pooled_ratio <- NA_real_
        }

        list(
            symbol=first(symbol),
            Kuster_class=first(Kuster_class),

            total_populations=.N,
            eligible_populations=n_eligible,
            eligible_population_fraction=
                n_eligible/total_populations,

            finite_population_ratios=
                n_finite_population_ratios,

            piN_mean=mean_piN,
            piS_mean=mean_piS,
            piN_piS_ratio_of_means=mean_ratio,

            piN_median=median_piN,
            piS_median=median_piS,

            piN_opportunity_weighted=pooled_piN,
            piS_opportunity_weighted=pooled_piS,
            piN_piS_opportunity_weighted=
                pooled_ratio,

            mean_callable_fraction=
                mean_callable_fraction,

            minimum_callable_fraction=
                minimum_callable_fraction,

            ratio_status=ratio_status
        )
    },
    by=gene_id
]

gene[
    ,
    pinpis_main_eligible :=
        eligible_populations >=
        minimum_eligible_populations
]

message("[join] Chromosomes")

bed <- fread(
    chromosome_file,
    header=FALSE,
    select=c(1,4)
)

setnames(
    bed,
    c("chromosome","gene_id")
)

bed <- unique(
    bed,
    by="gene_id"
)

gene <- merge(
    gene,
    bed,
    by="gene_id",
    all.x=TRUE,
    sort=FALSE
)

gene[
    ,
    chromosome_set := fifelse(
        chromosome %in% c("chrUn","chrXIX"),
        "excluded_chrUn_or_sex",
        fifelse(
            is.na(chromosome),
            "missing_chromosome",
            "autosome"
        )
    )
]

gene[
    ,
    autosome_only :=
        chromosome_set=="autosome"
]

message("[join] dN/dS primary set")

dnds <- fread(
    dnds_master_file,
    select="gene_id"
)

dnds <- unique(dnds)

dnds[
    ,
    dnds_primary_gene := TRUE
]

gene <- merge(
    gene,
    dnds,
    by="gene_id",
    all.x=TRUE,
    sort=FALSE
)

gene[
    is.na(dnds_primary_gene),
    dnds_primary_gene := FALSE
]

gene[
    ,
    analysis_set := fifelse(
        pinpis_main_eligible &
        dnds_primary_gene,
        "pi_main_and_dnds_matched",
        fifelse(
            pinpis_main_eligible,
            "pi_main_only",
            "pi_not_eligible"
        )
    )
]

message("[join] Original functional labels")

if (file.exists(old_label_file)) {
    old <- fread(
        old_label_file,
        na.strings=c("NA","")
    )

    wanted <- intersect(
        c(
            "gene_id",
            "old_name",
            "new_symbol",
            "own_role",
            "own_complex",
            "core_status",
            "old_target_gene"
        ),
        names(old)
    )

    old <- unique(
        old[
            ,
            ..wanted
        ],
        by="gene_id"
    )

    gene <- merge(
        gene,
        old,
        by="gene_id",
        all.x=TRUE,
        sort=FALSE
    )
}

setorder(
    gene,
    Kuster_class,
    gene_id
)

fwrite(
    gene,
    output_file,
    sep="\t",
    quote=FALSE,
    na="NA"
)

cat("\n===== Gene-level dimensions =====\n")
cat("Genes:",nrow(gene),"\n")

cat("\n===== πN/πS eligibility =====\n")
print(
    gene[
        ,
        .N,
        by=pinpis_main_eligible
    ]
)

cat("\n===== Eligible-population counts =====\n")
print(
    gene[
        ,
        .N,
        by=eligible_populations
    ][
        order(eligible_populations)
    ]
)

cat("\n===== Main eligible genes by Kuster class =====\n")
print(
    gene[
        pinpis_main_eligible==TRUE,
        .N,
        by=Kuster_class
    ]
)

cat("\n===== Main/matched/autosome counts =====\n")
print(
    gene[
        ,
        .(
            genes=.N,
            autosomal=sum(autosome_only),
            dnds_matched=sum(dnds_primary_gene)
        ),
        by=analysis_set
    ]
)

cat("\n===== Gene-level ratio status =====\n")
print(
    gene[
        ,
        .N,
        by=ratio_status
    ]
)

cat("\n===== Main eligible medians =====\n")

print(
    gene[
        pinpis_main_eligible==TRUE,
        .(
            n_genes=.N,
            median_piN=median(
                piN_mean,
                na.rm=TRUE
            ),
            median_piS=median(
                piS_mean,
                na.rm=TRUE
            ),
            ratio_identifiable_n=sum(
                is.finite(
                    piN_piS_ratio_of_means
                )
            ),
            median_piN_piS=median(
                piN_piS_ratio_of_means,
                na.rm=TRUE
            )
        ),
        by=Kuster_class
    ]
)

cat("\n===== Completely missing genes =====\n")

print(
    gene[
        eligible_populations==0,
        .(
            gene_id,
            symbol,
            Kuster_class,
            chromosome
        )
    ]
)

cat("\n[output]",output_file,"\n")
