#1. create bam.list
ls -1 /work/cyu/poolseq/PPalign_output/mtDNA_bam/*.bam > /work/cyu/poolseq/PPalign_output/mtDNA_bam/bamlist.txt

#2. merge BAM files (in the order of the file paths in BAMlist.txt) in a MPILEUP file only retaining nucleotides with BQ >20 and reads with MQ > 20
samtools mpileup -B \
    -f /work/cyu/sequence.fasta \
    -b bamlist.txt \
    -q 30 \
    -Q 30 \
    -d 5000 \
    | gzip > fish.mpileup.gz

perl mpileup2sync.pl \
    --input /work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.mpileup \
    --output /work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.sync \
    --fastq-type sanger \
    --min-qual 20

#3. ann and create bed files
/work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.sync

cd /work/cyu/poolseq/PPalign_output/ann



# grep gene/window/region/you want from gff file
grep -P "\tCDS\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | \
awk '{print $1, $4, $5}' OFS='\t' > pcg.bed

grep -P "\tgene\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | grep "rRNA" | \
awk '{print $1, $4, $5}' OFS='\t' > rrna.bed

grep -P "\tgene\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | grep "tRNA" | \
awk '{print $1, $4, $5}' OFS='\t' > trna.bed

grep -P "\tD[-_]loop\t|\tcontrol_region\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | awk '{print $1, $4, $5}' OFS='\t' > dloop.bed

grep -P "\tCDS\t" /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff | \
grep -w "CYTB" | \
awk '{print $1, $4, $5}' OFS='\t' > complex3.bed

sed -i 's/NC_041244.1/MH205729.1/g' rrna.bed
awk '{print $1, $2, $2+1, $0}' OFS='\t' /work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.sync > fish.sync.bed


#change header make bed header same as sync files
sed -i 's/NC_041244.1/MH205729.1/g' rrna.bed
#intersect bed and sync
bedtools intersect -a fish.sync.bed -b trna.bed > fish_trna.bed.sync
bedtools intersect -a fish.sync.bed -b ND6.bed > ND6.bed.sync
bedtools intersect -a fish.sync.bed -b complex1.bed > complex1.bed.sync

#make final sync files
#There are some errors might occur at this step, please deduplicate the sync files, some position might ouccr multiple times
awk '{print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33}' OFS='\t' complex1.bed.sync > complex1_fixed.sync
awk '{print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33}' OFS='\t' fish_rrna.bed.sync > fish_rrna_fixed.sync

#must read instruction of grenedalf and change conda env before running it
conda activate /home/cyu/y/envs/grenedalf

#grenedalf should change rename.txt header as same as the sync files name. I list some example scripts.

grenedalf fst \
    --method unbiased-nei \
    --sync-path /work/cyu/poolseq/PPalign_output/ann/fish_pcg_fixed.sync \
    --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting


grenedalf fst \
    --method unbiased-nei \
    --sync-path /work/cyu/poolseq/PPalign_output/mtDNA_bam/fish.sync\
    --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting

 grenedalf fst \
    --method unbiased-nei \
    --sync-path /work/cyu/poolseq/PPalign_output/ann/fish_rrna_fixed.sync\
    --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting   



 grenedalf fst \
    --method unbiased-nei \
    --sync-path /work/cyu/poolseq/PPalign_output/ann/complex5_fixed.sync\
    --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting 



/work/cyu/poolseq/PPalign_output/fst/
mv /work/cyu/poolseq/PPalign_output/fst/fst.csv /work/cyu/poolseq/PPalign_output/fst/grenedalf_ND6_fst.csv
grenedalf fst \
    --method unbiased-nei \
    --sync-path /work/cyu/poolseq/PPalign_output/ann/fish_trna_fixed.sync \
    --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting

grenedalf fst \
    --method unbiased-nei \
    --sync-path /work/cyu/poolseq/PPalign_output/ann/fish_dloop_fixed.sync \
    --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting

grenedalf fst \
    --method unbiased-hudson \
    --sync-path /work/cyu/poolseq/PPalign_output/ann/fish_pcg_fixed.sync \
    --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting



grenedalf fst \
    --method unbiased-nei \
    --sync-path /work/cyu/poolseq/PPalign_output/ann/ND6_fixed.sync \
    --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
    --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
    --window-type chromosomes \
    --window-average-policy valid-loci \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 4 \
    --no-extra-columns \
    --allow-file-overwriting

# scripts for diversity calculation involved pi, theta, tajima D
grenedalf diversity \
  --sync-path /work/cyu/poolseq/PPalign_output/ann/fish_pcg_fixed.sync \
  --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --window-type interval \
  --window-interval-width 500 \
  --window-interval-stride 100 \
  --window-average-policy valid-loci \
  --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
  --out-dir /work/cyu/poolseq/PPalign_output/diversity \
  --allow-file-overwriting 


grenedalf diversity \
  --sync-path /work/cyu/poolseq/PPalign_output/diversity/mt_genes_fixed.sync \
  --rename-samples-list /work/cyu/poolseq/PPalign_output/fst/rename.txt \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --window-type regions \
  --window-region-bed /work/cyu/poolseq/PPalign_output/diversity/mt_genes_new.bed \
  --window-average-policy valid-loci \
  --pool-sizes /work/cyu/poolseq/PPalign_output/fst/poolsize.txt \
  --out-dir /work/cyu/poolseq/PPalign_output/diversity \
  --allow-file-overwriting 


#pi plot and table
library(tidyverse)

# 1. read file
df <- read.csv(
  "/work/cyu/Poolinfo_with_13gene_details.csv",
  header = TRUE,
  check.names = FALSE
)

# 2. make long-format gene-level pi table
gene_pi_long <- df %>%
  dplyr::select(
    Population, Region, Habitat,
    pi_ND1, pi_ND2, pi_ND3, pi_ND4, pi_ND4L, pi_ND5, pi_ND6,
    pi_CYTB,
    pi_COX1, pi_COX2, pi_COX3,
    pi_ATP6, pi_ATP8
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("pi_"),
    names_to = "Gene",
    values_to = "pi"
  ) %>%
  dplyr::mutate(
    Gene = sub("^pi_", "", Gene),
    Complex = dplyr::case_when(
      Gene %in% c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") ~ "Complex I",
      Gene %in% c("CYTB") ~ "Complex III",
      Gene %in% c("COX1", "COX2", "COX3") ~ "Complex IV",
      Gene %in% c("ATP6", "ATP8") ~ "Complex V",
      TRUE ~ NA_character_
    )
  )

# 3. set biological order:
# left = Alaska (marine -> recent -> freshwater)
# right = BC (marine -> recent -> freshwater)
pop_order <- c(
  # Alaska
  "RS",                  # marine
  "SC", "CH",            # recent colonized
  "FG", "LG", "SR", "SL", "TL", "WB", "WT", "WK", "LB",  # freshwater
  
  # BC
  "SAY",                 # marine
  "FRED", "PACH",        # recent colonized
  "SWA", "THE", "JOE", "BEA", "MUC", "PYE", "AMO", "GOS", "ROB", "BOOT", "ECHO", "LAW"  # freshwater
)

gene_pi_long <- gene_pi_long %>%
  dplyr::mutate(
    Region = factor(Region, levels = c("Alaska", "BC")),
    Habitat = factor(Habitat, levels = c("Marine", "Recent Colonized", "Freshwater")),
    Population = factor(Population, levels = pop_order),
    Complex = factor(Complex, levels = c("Complex I", "Complex III", "Complex IV", "Complex V"))
  )

# 4. plot
p <- ggplot(gene_pi_long, aes(x = Population, y = pi, fill = Complex)) +
  geom_boxplot(
    position = position_dodge(width = 0.60),
    width = 0.70,
    linewidth = 0.30,
    outlier.size = 0.50
  ) +
  geom_jitter(
    aes(color = Complex),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.60),
    size = 0.70,
    alpha = 0.60,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(
    "Complex I"   = "#D55E00",
    "Complex III" = "#009E73",
    "Complex IV"  = "#0072B2",
    "Complex V"   = "#CC79A7"
  )) +
  scale_color_manual(values = c(
    "Complex I"   = "#D55E00",
    "Complex III" = "#009E73",
    "Complex IV"  = "#0072B2",
    "Complex V"   = "#CC79A7"
  )) +
  facet_grid(~Region, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(0, 0.018)) +
  labs(
    x = "Population",
    y = expression(pi),
    fill = "Complex"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    plot.title = element_text(size = 16, face = "bold"),
    
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "right"
  )
# 5. show plot
print(p)

# 6. save flat 2:1 figure
ggsave("complex_pi_flat_2to1.pdf", p, width = 12, height = 6)
ggsave("complex_pi_flat_2to1.png", p, width = 12, height = 6, dpi = 300)


