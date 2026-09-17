#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/mnt/spareHD_2/nu_287/q2_parallelism"

OUTDIR <- file.path(
  BASE,
  "q2_deltaAF_withAMO_noLB_Jesse"
)

IN_DELTA <- file.path(
  OUTDIR,
  "deltaAF_long.withAMO_noLB_Jesse.tsv.gz"
)

GENE_LIST <- file.path(
  BASE,
  "oxphos72_genes.withAMO_noLB_Jesse.list"
)

OUT_ALL <- file.path(
  OUTDIR,
  "allSNP_sharedScores_Fig3B.withAMO_noLB_Jesse.tsv.gz"
)

OUT_TOP <- file.path(
  OUTDIR,
  "Fig3B_sharedRepresentativeSNP.withAMO_noLB_Jesse.tsv"
)

OUT_LONG <- file.path(
  OUTDIR,
  "Fig3B_sharedSNP_populationData.withAMO_noLB_Jesse.tsv.gz"
)

OUT_SUM <- file.path(
  OUTDIR,
  "Fig3B_sharedSNP_regionSummary.withAMO_noLB_Jesse.tsv"
)

OUT_COVERAGE <- file.path(
  OUTDIR,
  "Fig3B_sharedSNP_coverage.withAMO_noLB_Jesse.tsv"
)

EXPECTED_N <- c(
  AK = 8L,
  BC = 12L
)

# ============================================================
# Read data
# ============================================================

if (!file.exists(IN_DELTA)) {
  stop("Missing input: ", IN_DELTA)
}

if (!file.exists(GENE_LIST)) {
  stop("Missing gene list: ", GENE_LIST)
}

DT <- fread(IN_DELTA)

genes72 <- readLines(
  GENE_LIST,
  warn = FALSE
)

genes72 <- unique(
  tolower(
    trimws(genes72)
  )
)

genes72 <- genes72[
  nzchar(genes72)
]

if (length(genes72) != 72L) {
  stop(
    "Expected 72 genes, found ",
    length(genes72)
  )
}

DT[, gene := tolower(trimws(gene))]
DT[, region := toupper(trimws(region))]

DT <- DT[
  gene %in% genes72 &
  region %in% c("AK", "BC") &
  is.finite(deltaAF) &
  is.finite(af) &
  is.finite(marine_af)
]

# ============================================================
# Regional score for every SNP
# ============================================================

REGIONAL <- DT[
  ,
  {
    med <- median(
      deltaAF,
      na.rm = TRUE
    )

    n_positive <- sum(
      deltaAF > 0,
      na.rm = TRUE
    )

    n_negative <- sum(
      deltaAF < 0,
      na.rm = TRUE
    )

    expected_n <- unname(
      EXPECTED_N[
        region[1]
      ]
    )

    list(
      n_pop = uniqueN(pop),

      expected_n_pop = expected_n,

      complete_coverage = (
        uniqueN(pop) == expected_n
      ),

      median_deltaAF = med,

      abs_median_deltaAF = abs(med),

      mean_deltaAF = mean(
        deltaAF,
        na.rm = TRUE
      ),

      abs_mean_deltaAF = abs(
        mean(
          deltaAF,
          na.rm = TRUE
        )
      ),

      directional_concordance = max(
        n_positive,
        n_negative
      ) / .N,

      mean_depth = mean(
        depth,
        na.rm = TRUE
      ),

      min_depth = min(
        depth,
        na.rm = TRUE
      ),

      marine_af = marine_af[1],

      marine_pop = marine_pop[1],

      marine_depth = marine_depth[1],

      focal_allele = focal_allele[1],

      chr = chr[1],

      pos = pos[1]
    )
  },
  by = .(
    region,
    gene,
    snp
  )
]

# ============================================================
# Put AK and BC values for the same SNP on one row
# ============================================================

WIDE <- dcast(
  REGIONAL,
  gene + snp + chr + pos + focal_allele ~ region,
  value.var = c(
    "n_pop",
    "complete_coverage",
    "median_deltaAF",
    "abs_median_deltaAF",
    "mean_deltaAF",
    "abs_mean_deltaAF",
    "directional_concordance",
    "mean_depth",
    "min_depth",
    "marine_af",
    "marine_pop",
    "marine_depth"
  )
)

# Require exact same SNP with complete AK and BC coverage

WIDE <- WIDE[
  complete_coverage_AK == TRUE &
  complete_coverage_BC == TRUE &
  is.finite(median_deltaAF_AK) &
  is.finite(median_deltaAF_BC)
]

# ============================================================
# Shared-SNP score
#
# shared score =
# min(
#   |AK median signed deltaAF|,
#   |BC median signed deltaAF|
# )
# ============================================================

WIDE[, same_nonzero_direction := (
  sign(median_deltaAF_AK) ==
    sign(median_deltaAF_BC) &
  sign(median_deltaAF_AK) != 0
)]

WIDE[, shared_score := pmin(
  abs_median_deltaAF_AK,
  abs_median_deltaAF_BC
)]

WIDE[, joint_mean_abs_median := (
  abs_median_deltaAF_AK +
  abs_median_deltaAF_BC
) / 2]

WIDE[, minimum_directional_concordance := pmin(
  directional_concordance_AK,
  directional_concordance_BC
)]

WIDE[, minimum_mean_depth := pmin(
  mean_depth_AK,
  mean_depth_BC
)]

# ============================================================
# Prefer same-direction SNPs
#
# If a gene contains one or more same-direction SNPs:
#   only those SNPs are eligible.
#
# If no same-direction SNP exists:
#   allow the highest shared-score opposite/zero SNP.
# ============================================================

WIDE[, has_concordant_candidate := any(
  same_nonzero_direction
), by = gene]

WIDE[, eligible_for_selection := fifelse(
  has_concordant_candidate,
  same_nonzero_direction,
  TRUE
)]

WIDE[, selection_class := fifelse(
  has_concordant_candidate,
  "concordant_preferred",
  "opposite_or_zero_fallback"
)]

fwrite(
  WIDE,
  OUT_ALL,
  sep = "\t",
  compress = "gzip"
)

# ============================================================
# Coverage report
# ============================================================

COVERAGE <- WIDE[
  ,
  .(
    n_shared_complete_snps = .N,

    n_concordant_shared_snps = sum(
      same_nonzero_direction
    )
  ),
  by = gene
]

COVERAGE <- merge(
  data.table(
    gene = genes72
  ),
  COVERAGE,
  by = "gene",
  all.x = TRUE
)

COVERAGE[
  is.na(n_shared_complete_snps),
  `:=`(
    n_shared_complete_snps = 0L,
    n_concordant_shared_snps = 0L
  )
]

setorder(
  COVERAGE,
  gene
)

fwrite(
  COVERAGE,
  OUT_COVERAGE,
  sep = "\t"
)

missing_genes <- COVERAGE[
  n_shared_complete_snps == 0L
]

if (nrow(missing_genes) > 0L) {
  print(missing_genes)

  stop(
    "Some genes have no exact SNP with complete coverage ",
    "in both AK and BC. See: ",
    OUT_COVERAGE
  )
}

# ============================================================
# Select shared SNP
#
# Primary:
#   maximum shared score
#
# Tie breakers:
#   mean regional effect
#   population directional concordance
#   sequencing depth
#   genomic position
# ============================================================

ELIGIBLE <- WIDE[
  eligible_for_selection == TRUE
]

setorder(
  ELIGIBLE,
  gene,
  -shared_score,
  -joint_mean_abs_median,
  -minimum_directional_concordance,
  -minimum_mean_depth,
  chr,
  pos
)

TOP <- ELIGIBLE[
  ,
  .SD[1],
  by = gene
]

TOP[, selection_method := paste0(
  "prefer_same_direction_then_",
  "maximize_minimum_absolute_regional_median_deltaAF"
)]

setorder(
  TOP,
  gene
)

# ============================================================
# Extract population-level rows for selected shared SNPs
# ============================================================

KEY <- TOP[
  ,
  .(
    gene,
    snp
  )
]

LONG <- merge(
  DT,
  KEY,
  by = c(
    "gene",
    "snp"
  ),
  all = FALSE
)

setorder(
  LONG,
  gene,
  region,
  pop
)

# ============================================================
# Regional summary for selected shared SNP
# ============================================================

SUMMARY <- LONG[
  ,
  {
    med <- median(
      deltaAF,
      na.rm = TRUE
    )

    n_positive <- sum(
      deltaAF > 0,
      na.rm = TRUE
    )

    n_negative <- sum(
      deltaAF < 0,
      na.rm = TRUE
    )

    list(
      n_pop = uniqueN(pop),

      median_deltaAF = med,

      abs_median_deltaAF = abs(med),

      mean_deltaAF = mean(
        deltaAF,
        na.rm = TRUE
      ),

      directional_concordance = max(
        n_positive,
        n_negative
      ) / .N,

      median_freshwater_af = median(
        af,
        na.rm = TRUE
      ),

      marine_af = marine_af[1],

      marine_pop = marine_pop[1],

      focal_allele = focal_allele[1],

      shared_snp = snp[1],

      chr = chr[1],

      pos = pos[1]
    )
  },
  by = .(
    gene,
    region
  )
]

SUMMARY <- merge(
  SUMMARY,
  TOP[
    ,
    .(
      gene,
      shared_score,
      same_nonzero_direction,
      selection_class
    )
  ],
  by = "gene",
  all.x = TRUE
)

setorder(
  SUMMARY,
  gene,
  region
)

# ============================================================
# Save
# ============================================================

fwrite(
  TOP,
  OUT_TOP,
  sep = "\t"
)

fwrite(
  LONG,
  OUT_LONG,
  sep = "\t",
  compress = "gzip"
)

fwrite(
  SUMMARY,
  OUT_SUM,
  sep = "\t"
)

# ============================================================
# Checks
# ============================================================

if (nrow(TOP) != 72L) {
  stop(
    "Expected 72 shared representative SNPs, found ",
    nrow(TOP)
  )
}

if (
  uniqueN(
    LONG,
    by = c(
      "gene",
      "snp"
    )
  ) != 72L
) {
  stop(
    "Population table does not contain exactly ",
    "one shared SNP per gene."
  )
}

cat(
  "\n[shared-SNP selection classes]\n"
)

print(
  TOP[
    ,
    .N,
    by = selection_class
  ]
)

cat(
  "\n[opposite/zero fallback genes]\n"
)

print(
  TOP[
    selection_class ==
      "opposite_or_zero_fallback",
    .(
      gene,
      snp,
      median_deltaAF_AK,
      median_deltaAF_BC,
      shared_score
    )
  ]
)

cat(
  "\n[OK] wrote:\n"
)

cat(
  "  ",
  OUT_ALL,
  "\n",
  sep = ""
)

cat(
  "  ",
  OUT_TOP,
  "\n",
  sep = ""
)

cat(
  "  ",
  OUT_LONG,
  "\n",
  sep = ""
)

cat(
  "  ",
  OUT_SUM,
  "\n",
  sep = ""
)

cat(
  "  ",
  OUT_COVERAGE,
  "\n",
  sep = ""
)
