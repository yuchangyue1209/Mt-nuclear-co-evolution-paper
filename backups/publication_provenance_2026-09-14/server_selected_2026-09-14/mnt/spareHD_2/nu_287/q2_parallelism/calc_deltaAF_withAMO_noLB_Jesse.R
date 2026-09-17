#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# Inputs and outputs
# ============================================================

AF_FILE <- paste0(
  "/mnt/spareHD_2/nu_287/q2_parallelism/",
  "af_long.withAMO_noLB_Jesse.tsv.gz"
)

OUTDIR <- paste0(
  "/mnt/spareHD_2/nu_287/q2_parallelism/",
  "q2_deltaAF_withAMO_noLB_Jesse"
)

dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)

OUT_DELTA <- file.path(
  OUTDIR,
  "deltaAF_long.withAMO_noLB_Jesse.tsv.gz"
)

OUT_SUM <- file.path(
  OUTDIR,
  "deltaAF_pop_summary.withAMO_noLB_Jesse.tsv"
)


# ============================================================
# Parameters
# ============================================================

MIN_DEPTH <- 20

DROP_POPS <- c(
  "LB"
)

AK_fresh <- c(
  "FG",
  "LG",
  "SR",
  "SL",
  "TL",
  "WB",
  "WT",
  "WK"
)

BC_fresh <- c(
  "SWA",
  "THE",
  "JOE",
  "BEA",
  "MUC",
  "PYE",
  "AMO",
  "BOOT",
  "ECHO",
  "LAW",
  "GOS",
  "ROB"
)

AK_marine <- "RS"
BC_marine <- "SAY"


# ============================================================
# Helper
# ============================================================

normalize_pop <- function(x) {

  x <- toupper(x)

  gsub(
    "^(\\d+_)?([A-Z]+)(?:_S\\d+)?$",
    "\\2",
    x,
    perl = TRUE
  )
}


# ============================================================
# Read allele-frequency table
# ============================================================

AF <- fread(AF_FILE)

required_columns <- c(
  "chr",
  "pos",
  "gene",
  "pop",
  "af",
  "depth",
  "focal_allele"
)

missing_columns <- setdiff(
  required_columns,
  names(AF)
)

if (length(missing_columns) > 0) {

  stop(
    "Missing required columns in AF file: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ============================================================
# Clean and filter
# ============================================================

AF[, pop := normalize_pop(pop)]

AF[, gene := tolower(
  trimws(gene)
)]

AF <- AF[
  depth >= MIN_DEPTH
]

AF <- AF[
  !pop %in% DROP_POPS
]

AF[, snp := paste(
  chr,
  pos,
  gene,
  sep = ":"
)]


# ============================================================
# Keep only populations used in the analysis
# ============================================================

keep_pops <- unique(
  c(
    AK_fresh,
    BC_fresh,
    AK_marine,
    BC_marine
  )
)

AF <- AF[
  pop %in% keep_pops
]


# ============================================================
# Assign region
# ============================================================

AF[, region := fifelse(
  pop %in% c(AK_fresh, AK_marine),
  "AK",
  fifelse(
    pop %in% c(BC_fresh, BC_marine),
    "BC",
    NA_character_
  )
)]

AF <- AF[
  !is.na(region)
]


# ============================================================
# Marine reference AF
# ============================================================

marine_AK <- AF[
  pop == AK_marine,
  .(
    region = "AK",
    snp,
    marine_pop = AK_marine,
    marine_af = af,
    marine_depth = depth
  )
]

marine_BC <- AF[
  pop == BC_marine,
  .(
    region = "BC",
    snp,
    marine_pop = BC_marine,
    marine_af = af,
    marine_depth = depth
  )
]

MAR <- rbindlist(
  list(
    marine_AK,
    marine_BC
  ),
  use.names = TRUE
)


# ============================================================
# Freshwater data
# ============================================================

FRESH <- AF[
  (
    region == "AK" &
      pop %in% AK_fresh
  ) |
    (
      region == "BC" &
        pop %in% BC_fresh
    ),
  .(
    region,
    snp,
    chr,
    pos,
    gene,
    pop,
    af,
    depth,
    focal_allele
  )
]


# ============================================================
# Join freshwater AF with its regional marine reference
# ============================================================

DT <- merge(
  FRESH,
  MAR,
  by = c(
    "region",
    "snp"
  ),
  all = FALSE
)

DT <- DT[
  is.finite(af) &
    is.finite(marine_af)
]


# ============================================================
# Calculate signed deltaAF
#
# AK:
#   freshwater AF - RS AF
#
# BC:
#   freshwater AF - SAY AF
# ============================================================

DT[, deltaAF := af - marine_af]

DT[, abs_deltaAF := abs(deltaAF)]


# ============================================================
# Reorder and save
# ============================================================

setcolorder(
  DT,
  c(
    "region",
    "snp",
    "chr",
    "pos",
    "gene",
    "pop",
    "deltaAF",
    "abs_deltaAF",
    "af",
    "marine_af",
    "marine_pop",
    "depth",
    "marine_depth",
    "focal_allele"
  )
)

setorder(
  DT,
  region,
  gene,
  snp,
  pop
)

fwrite(
  DT,
  OUT_DELTA,
  sep = "\t",
  compress = "gzip"
)


# ============================================================
# Population-level diagnostic summary
# ============================================================

SUM <- DT[
  ,
  .(
    n_snp_population_rows = .N,
    n_unique_snps = uniqueN(snp),
    mean_deltaAF = mean(
      deltaAF,
      na.rm = TRUE
    ),
    median_deltaAF = median(
      deltaAF,
      na.rm = TRUE
    ),
    mean_abs_deltaAF = mean(
      abs_deltaAF,
      na.rm = TRUE
    ),
    max_abs_deltaAF = max(
      abs_deltaAF,
      na.rm = TRUE
    )
  ),
  by = .(
    region,
    pop
  )
][
  order(
    region,
    pop
  )
]

fwrite(
  SUM,
  OUT_SUM,
  sep = "\t"
)


# ============================================================
# Checks
# ============================================================

cat(
  "[OK] wrote:\n",
  OUT_DELTA,
  "\n",
  OUT_SUM,
  "\n",
  sep = ""
)

cat(
  "\n[population counts]\n"
)

print(
  DT[
    ,
    .N,
    by = .(
      region,
      pop
    )
  ][
    order(
      region,
      pop
    )
  ]
)

cat(
  "\n[unique populations]\n"
)

print(
  sort(
    unique(DT$pop)
  )
)

if ("LB" %in% DT$pop) {
  stop("LB is still present in the deltaAF table.")
}

missing_AK <- setdiff(
  AK_fresh,
  unique(
    DT[
      region == "AK"
    ]$pop
  )
)

missing_BC <- setdiff(
  BC_fresh,
  unique(
    DT[
      region == "BC"
    ]$pop
  )
)

if (length(missing_AK) > 0) {

  stop(
    "Missing AK freshwater populations: ",
    paste(missing_AK, collapse = ", ")
  )
}

if (length(missing_BC) > 0) {

  stop(
    "Missing BC freshwater populations: ",
    paste(missing_BC, collapse = ", ")
  )
}

cat(
  "\n[OK] population check passed\n"
)

cat(
  "AK freshwater populations: ",
  length(AK_fresh),
  "\n",
  sep = ""
)

cat(
  "BC freshwater populations: ",
  length(BC_fresh),
  "\n",
  sep = ""
)
