#new af delta af
cat > /mnt/spareHD_2/nu_287/q2_parallelism/sync_to_deltaAF_noAMO_noLB.py <<'PY'
#!/usr/bin/env python3
import os, re, gzip, sys

SYNC_DIR = "/mnt/spareHD_2/nu_287/sync"
BAMLIST  = "/mnt/spareHD_2/oxphos_gene_tree/bamlist_nuclear.txt"
OUT_AF   = "/mnt/spareHD_2/nu_287/q2_parallelism/af_long.noAMO_noLB.tsv.gz"
OUT_META = "/mnt/spareHD_2/nu_287/q2_parallelism/snp_meta.noAMO_noLB.tsv.gz"

DROP = {"AMO", "LB"}
MIN_DEPTH = 10

def open_any(f):
    return gzip.open(f, "rt") if f.endswith(".gz") else open(f)

def norm_pop(x):
    x = os.path.basename(x)
    x = re.sub(r"_subset\.bam$", "", x)
    x = re.sub(r"^(?:\d+_)?([A-Za-z]+)(?:_S\d+)?$", r"\1", x)
    return x.upper()

def gene_from_fn(fn):
    fn = re.sub(r"\.sync(\.gz)?$", "", fn)
    return fn.split(".")[0]

def parse_sync_cell(cell):
    # PoPoolation2 sync order: A:T:C:G:N:del
    a = cell.split(":")
    A, T, C, G = map(int, a[:4])
    return A, T, C, G

pops_raw = [x.strip() for x in open(BAMLIST) if x.strip()]
pops = [norm_pop(x) for x in pops_raw]

keep_idx = [i for i,p in enumerate(pops) if p not in DROP]
keep_pops = [pops[i] for i in keep_idx]

print("[info] kept pops:", ",".join(keep_pops), file=sys.stderr)

outA = gzip.open(OUT_AF, "wt")
outM = gzip.open(OUT_META, "wt")

outA.write("chr\tpos\tgene\tpop\taf\tdepth\tfocal_allele\n")
outM.write("chr\tpos\tgene\tfocal_allele\n")

for fn in sorted(os.listdir(SYNC_DIR)):
    if not (fn.endswith(".sync") or fn.endswith(".sync.gz")):
        continue

    gene = gene_from_fn(fn)
    path = os.path.join(SYNC_DIR, fn)

    with open_any(path) as f:
        for line in f:
            if not line.strip():
                continue

            toks = line.rstrip("\n").split()
            if len(toks) < 3 + len(pops):
                continue

            chr_, pos, ref = toks[0], toks[1], toks[2]
            cells = toks[3:3+len(pops)]

            total = {"A":0, "T":0, "C":0, "G":0}
            per_pop = []

            for i in keep_idx:
                pop = pops[i]
                A,T,C,G = parse_sync_cell(cells[i])
                depth = A + T + C + G

                total["A"] += A
                total["T"] += T
                total["C"] += C
                total["G"] += G

                per_pop.append((pop, A, T, C, G, depth))

            focal = max(total.items(), key=lambda x: x[1])[0]

            if total[focal] == 0:
                continue

            outM.write(f"{chr_}\t{pos}\t{gene}\t{focal}\n")

            for pop,A,T,C,G,depth in per_pop:
                if depth < MIN_DEPTH:
                    continue

                count = {"A":A, "T":T, "C":C, "G":G}[focal]
                af = count / depth

                outA.write(
                    f"{chr_}\t{pos}\t{gene}\t{pop}\t{af:.6f}\t{depth}\t{focal}\n"
                )

outA.close()
outM.close()

print("[OK] wrote:", OUT_AF, file=sys.stderr)
print("[OK] wrote:", OUT_META, file=sys.stderr)
PY

python /mnt/spareHD_2/nu_287/q2_parallelism/sync_to_deltaAF_noAMO_noLB.py





#
cat > /mnt/spareHD_2/nu_287/q2_parallelism/calc_deltaAF_noAMO_noLB.R <<'RS'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

AF_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long.noAMO_noLB.tsv.gz"
OUTDIR  <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_noAMO_noLB"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_DELTA <- file.path(OUTDIR, "deltaAF_long.noAMO_noLB.tsv.gz")
OUT_SUM   <- file.path(OUTDIR, "deltaAF_pop_summary.noAMO_noLB.tsv")

MIN_DEPTH <- 10

AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","ROS","BOOT","ECHO","LAW","GOS","ROB")

AK_marine <- "RS"
BC_marine <- "SAY"

normalize_pop <- function(x){
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl=TRUE)
}

AF <- fread(cmd = paste("zcat", shQuote(AF_FILE)))
AF[, pop := normalize_pop(pop)]
AF[, gene := tolower(gene)]
AF <- AF[depth >= MIN_DEPTH]
AF <- AF[!pop %in% c("AMO","LB")]

AF[, snp := paste(chr, pos, gene, sep=":")]

keep_pops <- unique(c(AK_fresh, BC_fresh, AK_marine, BC_marine))
AF <- AF[pop %in% keep_pops]

AF[, region := fifelse(pop %in% c(AK_fresh, AK_marine), "AK",
                fifelse(pop %in% c(BC_fresh, BC_marine), "BC", NA_character_))]
AF <- AF[!is.na(region)]

marine_AK <- AF[pop == AK_marine, .(
  region = "AK",
  snp,
  marine_af = af,
  marine_depth = depth
)]

marine_BC <- AF[pop == BC_marine, .(
  region = "BC",
  snp,
  marine_af = af,
  marine_depth = depth
)]

MAR <- rbindlist(list(marine_AK, marine_BC), use.names = TRUE)

FRESH <- AF[
  (region == "AK" & pop %in% AK_fresh) |
  (region == "BC" & pop %in% BC_fresh),
  .(region, snp, chr, pos, gene, pop, af, depth, focal_allele)
]

DT <- merge(FRESH, MAR, by = c("region", "snp"), all.x = TRUE)
DT <- DT[is.finite(af) & is.finite(marine_af)]

DT[, deltaAF := af - marine_af]

setcolorder(DT, c(
  "region", "snp", "chr", "pos", "gene", "pop",
  "deltaAF", "af", "marine_af",
  "depth", "marine_depth", "focal_allele"
))

fwrite(DT, OUT_DELTA, sep = "\t", compress = "gzip")

SUM <- DT[, .(
  n_snps = .N,
  mean_deltaAF = mean(deltaAF, na.rm = TRUE),
  mean_abs_deltaAF = mean(abs(deltaAF), na.rm = TRUE),
  sd_deltaAF = sd(deltaAF, na.rm = TRUE),
  max_abs_deltaAF = max(abs(deltaAF), na.rm = TRUE)
), by = .(region, pop)][order(region, pop)]

fwrite(SUM, OUT_SUM, sep = "\t")

cat("[OK] wrote:\n")
cat("  ", OUT_DELTA, "\n")
cat("  ", OUT_SUM, "\n")
cat("\n[pop counts]\n")
print(DT[, .N, by = .(region, pop)][order(region, pop)])
RS

Rscript /mnt/spareHD_2/nu_287/q2_parallelism/calc_deltaAF_noAMO_noLB.R



#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Inputs
# ============================================================

AF_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"

OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_noAMO_noLB"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_DELTA <- file.path(OUTDIR, "deltaAF_long.noAMO_noLB.tsv.gz")
OUT_SUM   <- file.path(OUTDIR, "deltaAF_pop_summary.noAMO_noLB.tsv")

# ============================================================
# Parameters
# ============================================================

MIN_DEPTH <- 20

DROP_POPS <- c("AMO", "LB")

AK_fresh  <- c("FG","LG","SR","SL","TL","WB","WT","WK")
BC_fresh  <- c("SWA","THE","JOE","BEA","MUC","PYE","BOOT","ECHO","LAW","GOS","ROB")

AK_marine <- "RS"
BC_marine <- "SAY"

# ============================================================
# Helpers
# ============================================================

normalize_pop <- function(x) {
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl = TRUE)
}

# ============================================================
# Read AF
# ============================================================

AF <- fread(AF_FILE)

AF[, pop := normalize_pop(pop)]
AF[, gene := tolower(trimws(gene))]
AF <- AF[depth >= MIN_DEPTH]
AF <- AF[!pop %in% DROP_POPS]

AF[, snp := paste(chr, pos, gene, sep = ":")]

keep_pops <- unique(c(AK_fresh, BC_fresh, AK_marine, BC_marine))
AF <- AF[pop %in% keep_pops]

AF[, region := fifelse(
  pop %in% c(AK_fresh, AK_marine), "AK",
  fifelse(pop %in% c(BC_fresh, BC_marine), "BC", NA_character_)
)]

AF <- AF[!is.na(region)]

# ============================================================
# Marine references
# ============================================================

marine_AK <- AF[pop == AK_marine, .(
  region = "AK",
  snp,
  marine_af = af,
  marine_depth = depth
)]

marine_BC <- AF[pop == BC_marine, .(
  region = "BC",
  snp,
  marine_af = af,
  marine_depth = depth
)]

MAR <- rbindlist(list(marine_AK, marine_BC), use.names = TRUE)

# ============================================================
# Freshwater ΔAF
# ============================================================

FRESH <- AF[
  (region == "AK" & pop %in% AK_fresh) |
    (region == "BC" & pop %in% BC_fresh),
  .(region, snp, chr, pos, gene, pop, af, depth)
]

DT <- merge(FRESH, MAR, by = c("region", "snp"), all.x = TRUE)
DT <- DT[is.finite(af) & is.finite(marine_af)]

DT[, deltaAF := af - marine_af]
DT[, abs_deltaAF := abs(deltaAF)]

setcolorder(DT, c(
  "region", "snp", "chr", "pos", "gene", "pop",
  "deltaAF", "abs_deltaAF",
  "af", "marine_af",
  "depth", "marine_depth"
))

fwrite(DT, OUT_DELTA, sep = "\t", compress = "gzip")

SUM <- DT[, .(
  n_snps = .N,
  mean_deltaAF = mean(deltaAF, na.rm = TRUE),
  mean_abs_deltaAF = mean(abs(deltaAF), na.rm = TRUE),
  sd_deltaAF = sd(deltaAF, na.rm = TRUE),
  max_abs_deltaAF = max(abs(deltaAF), na.rm = TRUE)
), by = .(region, pop)][order(region, pop)]

fwrite(SUM, OUT_SUM, sep = "\t")

cat("[OK] wrote:\n")
cat("  ", OUT_DELTA, "\n", sep = "")
cat("  ", OUT_SUM, "\n", sep = "")

cat("\n[pop counts]\n")
print(DT[, .N, by = .(region, pop)][order(region, pop)])



#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

AF72 <- "/mnt/spareHD_2/nu_287/q2_parallelism/af_long_final_72genes_subunit_with_si.tsv.gz"
OUT  <- "/mnt/spareHD_2/nu_287/q2_parallelism/oxphos72_genes.list"

DT <- fread(AF72)

genes72 <- sort(unique(tolower(trimws(DT$gene))))

writeLines(genes72, OUT)

cat("[OK] wrote:", OUT, "\n")
cat("n genes:", length(genes72), "\n")




#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Inputs
# ============================================================

IN_DELTA <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_noAMO_noLB/deltaAF_long.noAMO_noLB.tsv.gz"
GENE_LIST_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/oxphos72_genes.list"

OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_noAMO_noLB"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_TOPSNP <- file.path(OUTDIR, "geneGlobalTopSNP_72OXPHOS.noAMO_noLB.tsv")
OUT_LONG   <- file.path(OUTDIR, "deltaAF_geneGlobalTopSNP_72OXPHOS.noAMO_noLB.tsv.gz")
OUT_SUM    <- file.path(OUTDIR, "deltaAF_geneGlobalTopSNP_72OXPHOS_regionSummary.noAMO_noLB.tsv")

# ============================================================
# Read
# ============================================================

DT <- fread(cmd = paste("zcat", shQuote(IN_DELTA)))
genes72 <- fread(GENE_LIST_FILE, header = FALSE)[[1]]

DT[, gene := tolower(trimws(gene))]
genes72 <- tolower(trimws(genes72))

DT <- DT[gene %in% genes72]
DT <- DT[is.finite(deltaAF)]

# ============================================================
# One fixed top SNP per gene
# Top SNP = largest max |ΔAF|
# Tie-breakers: mean |ΔAF|, then number of populations
# ============================================================

SNP_SCORE <- DT[, .(
  n_pop = uniqueN(pop),
  mean_abs_deltaAF = mean(abs(deltaAF), na.rm = TRUE),
  median_abs_deltaAF = median(abs(deltaAF), na.rm = TRUE),
  max_abs_deltaAF = max(abs(deltaAF), na.rm = TRUE),
  mean_depth = mean(depth, na.rm = TRUE),
  mean_marine_depth = mean(marine_depth, na.rm = TRUE),
  chr = chr[1],
  pos = pos[1]
), by = .(gene, snp)]

setorder(
  SNP_SCORE,
  gene,
  -max_abs_deltaAF,
  -mean_abs_deltaAF,
  -n_pop
)

TOP_SNP <- SNP_SCORE[, .SD[1], by = gene]

KEY <- TOP_SNP[, .(gene, snp)]

OUT <- merge(
  DT,
  KEY,
  by = c("gene", "snp"),
  all = FALSE
)

OUT[, abs_deltaAF := abs(deltaAF)]

setorder(OUT, region, pop, gene)

# ============================================================
# Region-level gene summary
# ============================================================

REGION_SUM <- OUT[, {
  n_pos <- sum(deltaAF > 0, na.rm = TRUE)
  n_neg <- sum(deltaAF < 0, na.rm = TRUE)
  maj_sign <- ifelse(n_pos >= n_neg, "+", "-")

  .(
    n_pop = .N,
    mean_deltaAF = mean(deltaAF, na.rm = TRUE),
    median_deltaAF = median(deltaAF, na.rm = TRUE),
    mean_abs_deltaAF = mean(abs(deltaAF), na.rm = TRUE),
    max_abs_deltaAF = max(abs(deltaAF), na.rm = TRUE),
    maj_sign = maj_sign,
    concordance = max(n_pos, n_neg) / .N,
    top_snp = snp[1],
    chr = chr[1],
    pos = pos[1]
  )
}, by = .(region, gene)]

setorder(REGION_SUM, region, gene)

# ============================================================
# Save
# ============================================================

fwrite(TOP_SNP, OUT_TOPSNP, sep = "\t")

fwrite(
  OUT,
  OUT_LONG,
  sep = "\t",
  compress = "gzip"
)

fwrite(
  REGION_SUM,
  OUT_SUM,
  sep = "\t"
)

# ============================================================
# Checks
# ============================================================

cat("[OK] wrote:\n")
cat("  ", OUT_TOPSNP, "\n", sep = "")
cat("  ", OUT_LONG, "\n", sep = "")
cat("  ", OUT_SUM, "\n", sep = "")

cat("\n[check]\n")
cat("genes kept:", uniqueN(OUT$gene), "\n")
cat("rows:", nrow(OUT), "\n")
cat("expected rows: 72 genes x 19 pops = 1368\n")

print(OUT[, .N, by = .(region, pop)][order(region, pop)])



#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# ============================================================
# Input / output
# ============================================================

IN_SUM <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_noAMO_noLB/deltaAF_geneGlobalTopSNP_72OXPHOS_regionSummary.noAMO_noLB.tsv"

OUTDIR <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_noAMO_noLB"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(OUTDIR, "Fig3_geneGlobalTopSNP_72OXPHOS_noAMO_noLB.png")
out_pdf <- file.path(OUTDIR, "Fig3_geneGlobalTopSNP_72OXPHOS_noAMO_noLB.pdf")
out_tab <- file.path(OUTDIR, "Fig3_geneGlobalTopSNP_72OXPHOS_AB_summary.tsv")

# ============================================================
# Parameters
# ============================================================

cut_dAF <- 0.5

col_shared <- "#5B7DB1"
col_region <- "#BFC9D9"
col_high   <- "#E49A8D"
col_bar    <- "#7F7F7F"
col_line   <- "#6A8DFF"

theme_pub <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10, colour = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      plot.tag = element_text(face = "bold", size = 15),
      plot.margin = margin(8, 10, 8, 10)
    )
}

# ============================================================
# Read and reshape
# ============================================================

G <- fread(IN_SUM)

AK <- G[region == "AK", .(
  gene,
  dAF_AK = mean_deltaAF,
  med_AK = median_deltaAF,
  abs_AK = mean_abs_deltaAF,
  max_AK = max_abs_deltaAF,
  sign_AK = maj_sign,
  concord_AK = concordance,
  top_snp_AK = top_snp
)]

BC <- G[region == "BC", .(
  gene,
  dAF_BC = mean_deltaAF,
  med_BC = median_deltaAF,
  abs_BC = mean_abs_deltaAF,
  max_BC = max_abs_deltaAF,
  sign_BC = maj_sign,
  concord_BC = concordance,
  top_snp_BC = top_snp
)]

AB <- merge(AK, BC, by = "gene")

# Same SNP check
if (!all(AB$top_snp_AK == AB$top_snp_BC)) {
  warning("Some genes have different AK/BC top SNP labels. Check input.")
}

# Parallel based on median direction
AB[, parallel := sign(med_AK) == sign(med_BC) & med_AK != 0 & med_BC != 0]
AB[, class := ifelse(parallel, "Shared parallel", "Discordant")]

AB[, high_effect :=
     parallel &
     abs(med_AK) >= cut_dAF &
     abs(med_BC) >= cut_dAF]

pearson <- cor.test(AB$dAF_AK, AB$dAF_BC, method = "pearson")

n_total <- nrow(AB)
n_shared <- sum(AB$parallel)
n_other <- n_total - n_shared
prop_shared <- n_shared / n_total
binom_p <- binom.test(n_shared, n_total, p = 0.5)$p.value
n_high <- sum(AB$high_effect)

fwrite(AB, out_tab, sep = "\t")

# ============================================================
# Panel A: mean ΔAF correlation
# ============================================================

panelA <- ggplot(AB, aes(dAF_AK, dAF_BC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey20") +
  geom_point(alpha = 0.75, size = 2.4, colour = "grey45") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, colour = col_line) +
  annotate(
    "text",
    x = min(AB$dAF_AK) + 0.05 * diff(range(AB$dAF_AK)),
    y = max(AB$dAF_BC) - 0.05 * diff(range(AB$dAF_BC)),
    hjust = 0,
    vjust = 1,
    label = sprintf(
      "n = %d genes\nPearson's r = %.2f\np = %.2g",
      n_total,
      unname(pearson$estimate),
      pearson$p.value
    ),
    size = 3.1
  ) +
  labs(
    x = expression("Mean AK " * Delta * AF * " per gene"),
    y = expression("Mean BC " * Delta * AF * " per gene")
  ) +
  theme_pub()

# ============================================================
# Panel B: median ΔAF shared direction
# ============================================================

panelB <- ggplot(AB, aes(med_AK, med_BC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey20") +
  geom_point(
    data = AB[class == "Discordant"],
    shape = 21,
    stroke = 0.5,
    size = 2.4,
    fill = "white",
    colour = col_region
  ) +
  geom_point(
    data = AB[class == "Shared parallel"],
    shape = 16,
    size = 2.5,
    colour = col_shared
  ) +
  annotate(
    "text",
    x = min(AB$med_AK) + 0.04 * diff(range(AB$med_AK)),
    y = max(AB$med_BC) - 0.04 * diff(range(AB$med_BC)),
    hjust = 0,
    vjust = 1,
    label = sprintf(
      "Genes = %d\nShared parallel = %d (%.1f%%)\nDiscordant = %d (%.1f%%)",
      n_total,
      n_shared,
      100 * prop_shared,
      n_other,
      100 * n_other / n_total
    ),
    size = 2.8
  ) +
  labs(
    x = expression("Median AK " * Delta * AF * " per gene"),
    y = expression("Median BC " * Delta * AF * " per gene")
  ) +
  theme_pub()

# ============================================================
# Panel C: counts
# ============================================================

bar_dt <- data.table(
  category = c("Discordant", "Shared parallel"),
  n = c(n_other, n_shared)
)

bar_dt[, prop := n / sum(n)]
bar_dt[, label := sprintf("%d\n(%.1f%%)", n, 100 * prop)]

panelC <- ggplot(bar_dt, aes(category, n)) +
  geom_col(width = 0.65, fill = col_bar) +
  geom_text(aes(label = label), vjust = -0.2, size = 3.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  annotate(
    "text",
    x = 0.78,
    y = max(bar_dt$n) * 1.3,
    hjust = 0,
    label = sprintf("Binomial p = %.2g", binom_p),
    size = 3.1
  ) +
  labs(
    x = NULL,
    y = "Number of genes"
  ) +
  theme_pub() +
  theme(legend.position = "none")

# ============================================================
# Panel D: high-effect shared genes
# ============================================================

panelD <- ggplot(AB, aes(med_AK, med_BC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey25") +
  geom_vline(xintercept = c(-cut_dAF, cut_dAF), linetype = "dotted", linewidth = 0.4, colour = "grey20") +
  geom_hline(yintercept = c(-cut_dAF, cut_dAF), linetype = "dotted", linewidth = 0.4, colour = "grey20") +
  geom_point(
    data = AB[class == "Discordant"],
    size = 2.0,
    colour = col_region,
    alpha = 0.7
  ) +
  geom_point(
    data = AB[class == "Shared parallel"],
    size = 2.2,
    colour = col_shared,
    alpha = 0.9
  ) +
  geom_point(
    data = AB[high_effect == TRUE],
    shape = 17,
    size = 3.3,
    colour = col_high
  ) +
  geom_text_repel(
    data = AB[high_effect == TRUE],
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25
  ) +
  annotate(
    "text",
    x = min(AB$med_AK) + 0.04 * diff(range(AB$med_AK)),
    y = max(AB$med_BC) - 0.04 * diff(range(AB$med_BC)),
    hjust = 0,
    vjust = 1,
    label = sprintf(
      "High-effect shared genes = %d\nabs(median ΔAF) >= %.1f in both regions",
      n_high,
      cut_dAF
    ),
    size = 2.8
  ) +
  labs(
    x = expression("Median AK " * Delta * AF * " per gene"),
    y = expression("Median BC " * Delta * AF * " per gene")
  ) +
  theme_pub()

# ============================================================
# Combine and save
# ============================================================

fig <- (panelA + panelB) / (panelC + panelD) +
  plot_annotation(tag_levels = "A")

print(fig)

ggsave(out_png, fig, width = 12, height = 9, dpi = 300)
ggsave(out_pdf, fig, width = 12, height = 9)

cat("[OK] wrote:\n")
cat("  ", out_png, "\n", sep = "")
cat("  ", out_pdf, "\n", sep = "")
cat("  ", out_tab, "\n", sep = "")