#norway
/work/cyu/poolseq/raw_data
#trim

#!/bin/bash

INPUT_DIR="/work/cyu/poolseq/raw_data"
CUTADAPT_DIR="/work/cyu/poolseq/PPalign_output/cutadapt"
THREADS=48

mkdir -p "$CUTADAPT_DIR"

for R1_FILE in "$INPUT_DIR"/*_1.fastq.gz; do
    R2_FILE="${R1_FILE/_1.fastq.gz/_2.fastq.gz}"
    PREFIX=$(basename "$R1_FILE" | sed 's/_1.fastq.gz//')

    OUT_R1="$CUTADAPT_DIR/cutadapt_R1_${PREFIX}.fastq.gz"
    OUT_R2="$CUTADAPT_DIR/cutadapt_R2_${PREFIX}.fastq.gz"
    LOG="$CUTADAPT_DIR/${PREFIX}_cutadapt.log"

    echo "Running cutadapt for $PREFIX..."

    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "$OUT_R1" \
        -p "$OUT_R2" \
        "$R1_FILE" "$R2_FILE" \
        -j "$THREADS" \
        > "$LOG" 2>&1
done

#!/bin/bash

INPUT_DIR="/work/cyu/poolseq/PPalign_output/cutadapt"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
THREADS=48

PREFIX="SRR27891779"

R1_FILE="$INPUT_DIR/cutadapt_R1_${PREFIX}.fastq.gz"
R2_FILE="$INPUT_DIR/cutadapt_R2_${PREFIX}.fastq.gz"

OUT_R1="$OUTPUT_DIR/trimmed_R1_${PREFIX}.fastq.gz"
OUT_R2="$OUTPUT_DIR/trimmed_R2_${PREFIX}.fastq.gz"
LOG_FILE="$OUTPUT_DIR/${PREFIX}_bbduk_log.txt"

mkdir -p "$OUTPUT_DIR"

bbduk.sh \
    in1="$R1_FILE" \
    in2="$R2_FILE" \
    out1="$OUT_R1" \
    out2="$OUT_R2" \
    trimq=20 \
    minlength=25 \
    ftl=10 \
    tossbrokenreads=t \
    threads="$THREADS" \
    > "$LOG_FILE" 2>&1

#mt bam
#!/bin/bash

# ==============================================================================
# Map Norway worm SRR27891779 trimmed reads to mt reference
# Then convert SAM -> BAM -> sorted BAM -> index
# ==============================================================================

TRIMMED_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
MAPPED_DIR="/work/cyu/poolseq/PPalign_output/mapped"
REFERENCE="/work/cyu/chrM_index"
THREADS=48

PREFIX="SRR27891779"

R1_FILE="$TRIMMED_DIR/trimmed_R1_${PREFIX}.fastq.gz"
R2_FILE="$TRIMMED_DIR/trimmed_R2_${PREFIX}.fastq.gz"

SAM_FILE="$MAPPED_DIR/${PREFIX}.sam"
BAM_FILE="$MAPPED_DIR/${PREFIX}.bam"
SORTED_BAM_FILE="$MAPPED_DIR/${PREFIX}_sorted.bam"

BOWTIE_LOG="$MAPPED_DIR/${PREFIX}_bowtie2.log"

mkdir -p "$MAPPED_DIR"

echo "Mapping $PREFIX..."

bowtie2 -x "$REFERENCE" \
    -1 "$R1_FILE" \
    -2 "$R2_FILE" \
    -p "$THREADS" \
    --very-sensitive-local \
    --no-mixed \
    --no-discordant \
    -X 2000 \
    -S "$SAM_FILE" \
    > "$BOWTIE_LOG" 2>&1

echo "Converting SAM to BAM with MAPQ >= 20..."

samtools view -b -q 20 "$SAM_FILE" > "$BAM_FILE"

echo "Sorting BAM..."

samtools sort -o "$SORTED_BAM_FILE" "$BAM_FILE"

echo "Indexing sorted BAM..."

samtools index "$SORTED_BAM_FILE"

echo "Finished processing $PREFIX"

echo "Mapping rate:"
grep "overall alignment rate" "$BOWTIE_LOG"

echo "Output files:"
ls -lh "$MAPPED_DIR"/${PREFIX}*




#
mkdir -p /work/cyu/poolseq/PPalign_output/norway_mt_sync

samtools mpileup -B \
  -f /work/cyu/sequence.fasta \
  /work/cyu/poolseq/PPalign_output/mapped/SRR27891779_sorted.bam \
  -q 30 \
  -Q 30 \
  -d 5000 \
  > /work/cyu/poolseq/PPalign_output/norway_mt_sync/SRR27891779.mt.mpileup

perl mpileup2sync.pl \
  --input /work/cyu/poolseq/PPalign_output/norway_mt_sync/SRR27891779.mt.mpileup \
  --output /work/cyu/poolseq/PPalign_output/norway_mt_sync/SRR27891779.mt.sync \
  --fastq-type sanger \
  --min-qual 20

cd /work/cyu/poolseq/PPalign_output/norway_mt_sync


awk '{
  print $1 "\t" $2-1 "\t" $2 "\t" $0
}' OFS='\t' SRR27891779.mt.sync > SRR27891779.mt.sync.bed

bedtools intersect \
  -a SRR27891779.mt.sync.bed \
  -b /work/cyu/poolseq/PPalign_output/ann/dloop.bed \
  -v \
  > SRR27891779.mt.noDloop.bed.sync

cut -f4- SRR27891779.mt.noDloop.bed.sync > SRR27891779.mt.noDloop.sync

awk 'NR<=5{print NF}' /work/cyu/poolseq/PPalign_output/norway_mt_sync/SRR27891779.mt.noDloop.sync
head /work/cyu/poolseq/PPalign_output/norway_mt_sync/SRR27891779.mt.noDloop.sync


#
# 1. rebuild poolsize_withNorway.txt
cp /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
   /work/cyu/poolseq/PPalign_output/fst/poolsize_withNorway.txt

echo "Norway,50" >> /work/cyu/poolseq/PPalign_output/fst/poolsize_withNorway.txt

# 2. check
tail /work/cyu/poolseq/PPalign_output/fst/poolsize_withNorway.txt
wc -l /work/cyu/poolseq/PPalign_output/fst/poolsize_withNorway.txt

# 3. rerun mt FST
conda activate /home/cyu/y/envs/grenedalf

SYNC="/work/cyu/poolseq/PPalign_output/ann/fish_noDloop_withNorway.sync"
POOL="/work/cyu/poolseq/PPalign_output/fst/poolsize_withNorway.txt"
RENAME="/work/cyu/poolseq/PPalign_output/fst/rename_mt_noDloop_withNorway.txt"
OUTDIR="/work/cyu/poolseq/PPalign_output/fst_mt_noDloop_withNorway"

mkdir -p "$OUTDIR"

awk -F',' '{print "fish_noDloop_withNorway." NR "," $1}' "$POOL" > "$RENAME"

cd "$OUTDIR"

grenedalf fst \
  --method unbiased-nei \
  --sync-path "$SYNC" \
  --rename-samples-list "$RENAME" \
  --pool-sizes "$POOL" \
  --window-type chromosomes \
  --window-average-policy valid-loci \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --no-extra-columns \
  --allow-file-overwriting

mv fst.csv mtgenome_noDloop_withNorway_fst.csv

head mtgenome_noDloop_withNorway_fst.csv





#nu gene merged 

cd /mnt/spareHD_2/nuclear_marked_duplicates



samtools mpileup -B \
    -f /work/cyu/stickleback_nuclear_only.fa \
    /mnt/spareHD_2/nuclear_marked_duplicates/SRR27891779_dedup.bam \
    -q 30 \
    -Q 30 \
    -d 5000 \
    | gzip > /work/cyu/norway_nuclear.mpileup.gz

gunzip norway_nuclear.mpileup.gz


perl mpileup2sync.pl \
  --input /mnt/spareHD_2/nuclear_marked_duplicates/norway_nuclear.mpileup \
  --output /work/cyu/norway_nuclear.sync \
  --fastq-type sanger \
  --min-qual 20



#!/bin/bash
set -euo pipefail

OLD="/mnt/spareHD_2/nuclear_marked_duplicates/nuclear.sync"
NOR="/work/cyu/norway_nuclear.sync"
OUT="/work/cyu/nuclear_withNorway.merged.sync"

echo "Merging Norway into old nuclear.sync..."

awk '
BEGIN{OFS="\t"}
FNR==NR{
    key=$1 FS $2 FS $3
    norway[key]=$4
    next
}
{
    key=$1 FS $2 FS $3
    if (key in norway) {
        print $0, norway[key]
    } else {
        print $0, "0:0:0:0:0:0"
    }
}
' "$NOR" "$OLD" > "$OUT"

echo "Done:"
ls -lh "$OUT"

echo "Check line numbers:"
wc -l "$OLD" "$NOR" "$OUT"

echo "Check columns:"
awk 'NR<=5{print NF}' "$OLD"
awk 'NR<=5{print NF}' "$OUT"

echo "Check Norway last column:"
head "$OUT" | awk '{print NF, $1, $2, $3, $NF}'

echo "Check Norway missing ratio:"
awk '
{
    total++
    if ($NF=="0:0:0:0:0:0") miss++
}
END{
    print "total =", total
    print "missing =", miss
    print "missing_ratio =", miss/total
}
' "$OUT"




#!/bin/bash
set -euo pipefail

# ============================================================
# Nuclear gene-level FST with Norway
# ============================================================

SYNC="/work/cyu/nuclear_withNorway.merged.sync"
GENE_BED="/work/cyu/gene_regions.bed"
WORKDIR="/work/cyu/gene_fst_work_withNorway"
POOL="/work/cyu/poolseq/PPalign_output/fst/poolsize_withNorway.txt"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# 1. sync -> bed-like
awk '{
  print $1 "\t" $2-1 "\t" $2 "\t" $0
}' OFS='\t' "$SYNC" > nuclear_withNorway.sync.bed

# 2. intersect with gene regions
bedtools intersect \
  -a nuclear_withNorway.sync.bed \
  -b "$GENE_BED" \
  -wa -wb > nuclear_gene_withNorway.bed.sync

# 3. split by gene
mkdir -p gene_sync_raw

awk '{
  gene=$NF
  print $0 >> "gene_sync_raw/" gene ".bed.sync"
}' nuclear_gene_withNorway.bed.sync

# 4. recover sync columns dynamically: chr pos ref + all samples
mkdir -p gene_sync_fixed

for f in gene_sync_raw/*.bed.sync; do
  base=$(basename "$f" .bed.sync)

  awk '{
    for (i=4; i<=NF-4; i++) {
      printf "%s%s", $i, (i<NF-4 ? OFS : ORS)
    }
  }' OFS='\t' "$f" > "gene_sync_fixed/${base}.sync"
done

# 5. deduplicate sites
mkdir -p gene_sync_dedup

for f in gene_sync_fixed/*.sync; do
  base=$(basename "$f" .sync)
  awk '!seen[$1 FS $2]++' "$f" > "gene_sync_dedup/${base}.sync"
done

# 6. run grenedalf FST
conda activate /home/cyu/y/envs/grenedalf

OUTDIR="$WORKDIR/fst_out_all_withNorway"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

for f in "$WORKDIR"/gene_sync_dedup/*.sync; do
  gene=$(basename "$f" .sync)

  awk -F',' -v g="$gene" '{print g "." NR "\t" $1}' \
    "$POOL" > "${gene}.rename.txt"

  grenedalf fst \
    --method unbiased-nei \
    --sync-path "$f" \
    --rename-samples-list "${gene}.rename.txt" \
    --pool-sizes "$POOL" \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting

  mv fst.csv "${gene}_fst.csv"
done

echo "Done. Number of gene FST files:"
ls "$OUTDIR"/*_fst.csv | wc -l


conda activate /home/cyu/y/envs/grenedalf

WORKDIR="/work/cyu/gene_fst_work_withNorway"
POOL="/work/cyu/poolseq/PPalign_output/fst/poolsize_withNorway.txt"
OUTDIR="$WORKDIR/fst_out_all_withNorway"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

for f in "$WORKDIR"/gene_sync_dedup/*.sync; do
  gene=$(basename "$f" .sync)

  echo "Running $gene"

  awk -F',' -v g="$gene" '{print g "." NR "\t" $1}' \
    "$POOL" > "${gene}.rename.txt"

  grenedalf fst \
    --method unbiased-nei \
    --sync-path "$f" \
    --rename-samples-list "${gene}.rename.txt" \
    --pool-sizes "$POOL" \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --threads 32 \
    --allow-file-overwriting

  mv fst.csv "${gene}_fst.csv"
done


#pbs
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# paths: with Norway
# =========================
MT_FST_FILE <- "/work/cyu/poolseq/PPalign_output/fst_mt_noDloop_withNorway/mtgenome_noDloop_withNorway_fst.csv"
NU_FST_DIR  <- "/work/cyu/gene_fst_work_withNorway/fst_out_all_withNorway"

OUTDIR <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# =========================
# population design
# =========================
AK_FW <- c("FG","LG","SR","SL","TL","WB","WT","WK","LB")
BC_FW <- c("SWA","THE","JOE","BEA","MUC","PYE","AMO","BOOT","ECHO","LAW","GOS","ROB")

design <- rbindlist(list(
  data.table(focal = AK_FW, region = "AK", marine = "RS",  outgroup = "Norway"),
  data.table(focal = BC_FW, region = "BC", marine = "SAY", outgroup = "Norway")
))

design <- unique(design)

# =========================
# helper functions
# =========================
fst_to_matrix <- function(file) {
  df <- fread(file)
  if (nrow(df) < 1 || ncol(df) < 4) stop("Bad FST file: ", file)

  row <- df[1]
  pair_cols <- colnames(df)[4:ncol(df)]
  vals <- as.numeric(unlist(row[, 4:ncol(df), with = FALSE]))

  pair_names <- sub("\\.fst$", "", pair_cols)
  sp <- tstrsplit(pair_names, ":", fixed = TRUE)

  pop1 <- sp[[1]]
  pop2 <- sp[[2]]
  pops <- sort(unique(c(pop1, pop2)))

  mat <- matrix(NA_real_, nrow = length(pops), ncol = length(pops),
                dimnames = list(pops, pops))
  diag(mat) <- 0

  for (i in seq_along(vals)) {
    mat[pop1[i], pop2[i]] <- vals[i]
    mat[pop2[i], pop1[i]] <- vals[i]
  }

  mat
}

get_fst <- function(mat, a, b) {
  if (!(a %in% rownames(mat)) || !(b %in% colnames(mat))) return(NA_real_)
  mat[a, b]
}

calc_pbs_one <- function(mat, focal, marine, outgroup) {
  fst_AB <- get_fst(mat, focal, marine)
  fst_AC <- get_fst(mat, focal, outgroup)
  fst_BC <- get_fst(mat, marine, outgroup)

  fst_AB <- pmin(pmax(fst_AB, 0), 0.999999)
  fst_AC <- pmin(pmax(fst_AC, 0), 0.999999)
  fst_BC <- pmin(pmax(fst_BC, 0), 0.999999)

  if (!all(is.finite(c(fst_AB, fst_AC, fst_BC)))) {
    return(data.table(
      fst_focal_marine = fst_AB,
      fst_focal_outgroup = fst_AC,
      fst_marine_outgroup = fst_BC,
      T_focal_marine = NA_real_,
      T_focal_outgroup = NA_real_,
      T_marine_outgroup = NA_real_,
      PBS = NA_real_
    ))
  }

  T_AB <- -log(1 - fst_AB)
  T_AC <- -log(1 - fst_AC)
  T_BC <- -log(1 - fst_BC)

  PBS <- (T_AB + T_AC - T_BC) / 2

  data.table(
    fst_focal_marine = fst_AB,
    fst_focal_outgroup = fst_AC,
    fst_marine_outgroup = fst_BC,
    T_focal_marine = T_AB,
    T_focal_outgroup = T_AC,
    T_marine_outgroup = T_BC,
    PBS = PBS
  )
}

calc_pbs_table <- function(mat, design_dt) {
  rbindlist(lapply(seq_len(nrow(design_dt)), function(i) {
    d <- design_dt[i]
    out <- calc_pbs_one(mat, d$focal, d$marine, d$outgroup)
    cbind(d, out)
  }))
}

# =========================
# 1. mt PBS
# =========================
cat("[read] mt FST with Norway\n")
mt_mat <- fst_to_matrix(MT_FST_FILE)

design_use <- design[
  focal %in% rownames(mt_mat) &
    marine %in% rownames(mt_mat) &
    outgroup %in% rownames(mt_mat)
]

cat("[info] design rows used:", nrow(design_use), "\n")
print(design_use)

mt_pbs <- calc_pbs_table(mt_mat, design_use)
setnames(mt_pbs, "PBS", "mt_PBS")

fwrite(mt_pbs, file.path(OUTDIR, "mt_noDloop_PBS_NorwayOutgroup_by_population.tsv"), sep = "\t")

# =========================
# 2. nuclear gene PBS
# =========================
files <- list.files(NU_FST_DIR, pattern = "_fst\\.csv$", full.names = TRUE)
cat("[info] nuclear gene FST files:", length(files), "\n")

nu_list <- vector("list", length(files))

for (i in seq_along(files)) {
  f <- files[i]
  gene <- sub("_fst\\.csv$", "", basename(f))

  cat("[run] ", i, "/", length(files), " ", gene, "\n", sep = "")

  nu_list[[i]] <- tryCatch({
    mat <- fst_to_matrix(f)
    tmp <- calc_pbs_table(mat, design_use)
    tmp[, gene := gene]
    tmp
  }, error = function(e) {
    data.table(
      focal = NA_character_,
      region = NA_character_,
      marine = NA_character_,
      outgroup = NA_character_,
      gene = gene,
      fst_focal_marine = NA_real_,
      fst_focal_outgroup = NA_real_,
      fst_marine_outgroup = NA_real_,
      T_focal_marine = NA_real_,
      T_focal_outgroup = NA_real_,
      T_marine_outgroup = NA_real_,
      PBS = NA_real_,
      error = conditionMessage(e)
    )
  })
}

nu_pbs <- rbindlist(nu_list, use.names = TRUE, fill = TRUE)
setnames(nu_pbs, "PBS", "nu_PBS")

fwrite(nu_pbs, file.path(OUTDIR, "nuclear_gene_PBS_NorwayOutgroup_by_population.tsv"), sep = "\t")

# =========================
# 3. merge mt and nuclear PBS
# =========================
merged <- merge(
  nu_pbs,
  mt_pbs[, .(focal, region, marine, outgroup, mt_PBS)],
  by = c("focal", "region", "marine", "outgroup"),
  all.x = TRUE
)

fwrite(merged, file.path(OUTDIR, "merged_mt_nuclear_PBS_NorwayOutgroup_by_gene_population.tsv"), sep = "\t")

# =========================
# 4. gene-wise correlation
# =========================
gene_res <- merged[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))
      data.table(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    } else {
      data.table(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    }
  },
  by = gene
]

gene_res[, p_bh := p.adjust(p, method = "BH")]
setorder(gene_res, p, -cor_r)

fwrite(gene_res, file.path(OUTDIR, "gene_correlation_nuPBS_vs_mtPBS_NorwayOutgroup.tsv"), sep = "\t")

# =========================
# 5. region-specific correlation
# =========================
gene_region_res <- merged[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 4 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))
      data.table(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    } else {
      data.table(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS, na.rm = TRUE),
        mean_mt_PBS = mean(mt_PBS, na.rm = TRUE)
      )
    }
  },
  by = .(gene, region)
]

gene_region_res[, p_bh := p.adjust(p, method = "BH"), by = region]
setorder(gene_region_res, region, p, -cor_r)

fwrite(gene_region_res, file.path(OUTDIR, "gene_region_correlation_nuPBS_vs_mtPBS_NorwayOutgroup.tsv"), sep = "\t")

# =========================
# 6. plots
# =========================
plot_dt <- gene_res[is.finite(cor_r) & is.finite(p)]

if (nrow(plot_dt) > 0) {
  p1 <- ggplot(plot_dt, aes(x = cor_r, y = -log10(p))) +
    geom_point() +
    theme_classic(base_size = 14) +
    labs(
      x = "Correlation between nuclear gene PBS and mt PBS",
      y = expression(-log[10](p)),
      title = "Genes whose PBS trajectories match mitochondrial PBS\nOutgroup = Norway"
    )

  ggsave(
    file.path(OUTDIR, "Fig_gene_nuPBS_vs_mtPBS_NorwayOutgroup_volcano.png"),
    p1, width = 7, height = 5, dpi = 300
  )
}

cat("\n[OK] done\n")
cat("[outdir] ", OUTDIR, "\n", sep = "")

#exclude am in oxphos 72 pbs
cd /work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup

awk -F'\t' 'NR==1 || $1!="AMO"' \
OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_by_gene_population.tsv \
> OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_noAMO.tsv




library(data.table)

dt <- fread("/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_noAMO.tsv")

res <- dt[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      ct <- suppressWarnings(cor.test(nu_PBS, mt_PBS, method = "pearson"))
      .(
        n_pop = .N,
        cor_r = unname(ct$estimate),
        p = ct$p.value,
        mean_nu_PBS = mean(nu_PBS),
        mean_mt_PBS = mean(mt_PBS)
      )
    } else {
      .(
        n_pop = .N,
        cor_r = NA_real_,
        p = NA_real_,
        mean_nu_PBS = mean(nu_PBS),
        mean_mt_PBS = mean(mt_PBS)
      )
    }
  },
  by = gene
]

res[, p_bh := p.adjust(p, method = "BH")]
setorder(res, p, -cor_r)

fwrite(
  res,
  "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/OXPHOS72_gene_correlation_nuPBS_vs_mtPBS_NorwayOutgroup_noAMO.tsv",
  sep = "\t"
)

print(head(res, 20))


#perm
library(data.table)

# =========================
# INPUT
# =========================
FILE <- "OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_by_gene_population.tsv"
OUT  <- "OXPHOS72_gene_correlation_perm_noAMO.tsv"

dt <- fread(FILE)

# =========================
# remove AMO
# =========================
dt <- dt[focal != "AMO"]

set.seed(123)

# =========================
# permutation function
# =========================
perm_test <- function(x, y, nperm = 10000) {
  r_obs <- cor(x, y)

  r_perm <- replicate(nperm, cor(x, sample(y)))
  p <- mean(abs(r_perm) >= abs(r_obs))

  list(r = r_obs, p = p)
}

# =========================
# run per gene
# =========================
res <- dt[
  is.finite(nu_PBS) & is.finite(mt_PBS),
  {
    if (.N >= 5 && var(nu_PBS) > 0 && var(mt_PBS) > 0) {
      pt <- perm_test(nu_PBS, mt_PBS)

      .(
        n_pop = .N,
        cor_r = pt$r,
        p_perm = pt$p,
        mean_nu_PBS = mean(nu_PBS),
        mean_mt_PBS = mean(mt_PBS)
      )
    } else {
      .(
        n_pop = .N,
        cor_r = NA_real_,
        p_perm = NA_real_,
        mean_nu_PBS = mean(nu_PBS),
        mean_mt_PBS = mean(mt_PBS)
      )
    }
  },
  by = gene
]

# =========================
# FDR correction
# =========================
res[, p_bh := p.adjust(p_perm, method = "BH")]

setorder(res, p_perm, -cor_r)

# =========================
# save
# =========================
fwrite(res, OUT, sep = "\t")

cat("\n=== TOP genes (perm, no AMO) ===\n")
print(head(res, 20))





#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ==============================================================================
# 1. Paths: Norway outgroup, OXPHOS72, no AMO
# ==============================================================================
INPUT_FILE <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/OXPHOS72_merged_mt_nuclear_PBS_NorwayOutgroup_noAMO.tsv"
ANNOT_FILE <- "/work/cyu/codeml_sites_summary.merged.tsv"

OUTDIR <- "/work/cyu/gene_fst_work_withNorway/pbs_mt_vs_nuclear_NorwayOutgroup/validation_results_OXPHOS72_noAMO"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. Read data
# ==============================================================================
cat("[Loading] Reading merged PBS data...\n")
dt <- fread(INPUT_FILE)

cat("[Loading] Reading gene annotation...\n")
annot <- fread(ANNOT_FILE)

annot_gene <- unique(annot[, .(gene, role, complex)])

# ==============================================================================
# 3. Calculate region-wise mean PBS
# ==============================================================================
gene_summary <- dt[is.finite(nu_PBS), .(
  mean_PBS = mean(nu_PBS, na.rm = TRUE),
  n_pop = .N
), by = .(gene, region)]

ak_dt <- gene_summary[region == "AK", .(
  gene,
  mean_PBS_AK = mean_PBS,
  n_pop_AK = n_pop
)]

bc_dt <- gene_summary[region == "BC", .(
  gene,
  mean_PBS_BC = mean_PBS,
  n_pop_BC = n_pop
)]

compare_dt <- merge(ak_dt, bc_dt, by = "gene")

# ==============================================================================
# 4. Add annotation and keep subunits only
# ==============================================================================
compare_dt <- merge(compare_dt, annot_gene, by = "gene", all.x = TRUE)

sub_dt <- compare_dt[role == "subunit"]

cat("[Info] Number of subunit genes:", nrow(sub_dt), "\n")

# ==============================================================================
# 5. Define shared vs region-specific status among subunits only
# ==============================================================================
cutoff_ak <- quantile(sub_dt$mean_PBS_AK, 0.95, na.rm = TRUE)
cutoff_bc <- quantile(sub_dt$mean_PBS_BC, 0.95, na.rm = TRUE)

cat("[Info] AK 95% cutoff:", cutoff_ak, "\n")
cat("[Info] BC 95% cutoff:", cutoff_bc, "\n")

sub_dt[, status := "Background"]

sub_dt[
  mean_PBS_AK > cutoff_ak & mean_PBS_BC > cutoff_bc,
  status := "Shared_Selection"
]

sub_dt[
  mean_PBS_AK > cutoff_ak & mean_PBS_BC <= cutoff_bc,
  status := "AK_Specific"
]

sub_dt[
  mean_PBS_AK <= cutoff_ak & mean_PBS_BC > cutoff_bc,
  status := "BC_Specific"
]

cat("[Info] Status counts:\n")
print(sub_dt[, .N, by = status])

# ==============================================================================
# 6. Save table
# ==============================================================================
OUT_TABLE <- file.path(OUTDIR, "OXPHOS72_noAMO_shared_vs_specific_PBS_NorwayOutgroup.tsv")
fwrite(sub_dt, OUT_TABLE, sep = "\t")
cat("[Write] Table:", OUT_TABLE, "\n")

# ==============================================================================
# 7. Plot
# ==============================================================================
p1 <- ggplot(sub_dt, aes(x = mean_PBS_BC, y = mean_PBS_AK, color = status)) +
  geom_point(size = 2.4, alpha = 0.75) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_text_repel(
    data = sub_dt[status != "Background"],
    aes(label = gene),
    size = 3.2,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "AK_Specific" = "blue",
    "BC_Specific" = "green",
    "Shared_Selection" = "purple",
    "Background" = "grey80"
  )) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Mean PBS (British Columbia, SAY vs Norway)",
    y = "Mean PBS (Alaska, RS vs Norway)",
    color = "status",
  )

OUT_PNG <- file.path(OUTDIR, "Fig_OXPHOS72_noAMO_AK_vs_BC_PBS_NorwayOutgroup.png")
OUT_PDF <- file.path(OUTDIR, "Fig_OXPHOS72_noAMO_AK_vs_BC_PBS_NorwayOutgroup.pdf")

ggsave(OUT_PNG, p1, width = 8, height = 7, dpi = 300)
ggsave(OUT_PDF, p1, width = 8, height = 7)

cat("[Write] Figure PNG:", OUT_PNG, "\n")
cat("[Write] Figure PDF:", OUT_PDF, "\n")

cat("\n[Done]\n")
print(p1)




#pbsn1
#max pbs
