3 genes shift

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================================================
# 1. paths
# =========================================================
DELTA_FILE   <- "/mnt/spareHD_2/nu_287/q2_parallelism/q2_deltaAF_SNPlevel_vs_mtPC_noAMO/deltaAF_long.noAMO.tsv.gz"
CLUSTER_FILE <- "/mnt/spareHD_2/nu_287/q2_parallelism/mtCluster_manual.tsv"
OUTDIR       <- "/mnt/spareHD_2/nu_287/q2_parallelism/_figs_table6_gene_shift_v2"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

target_genes <- c("ndufa4", "hccsb", "ndufb4")

# =========================================================
# 2. helper
# =========================================================
normalize_pop <- function(x){
  x <- toupper(x)
  gsub("^(\\d+_)?([A-Z]+)(?:_S\\d+)?$", "\\2", x, perl = TRUE)
}

# =========================================================
# 3. read data
# =========================================================
DEL <- fread(cmd = paste("zcat", shQuote(DELTA_FILE)))
DEL[, pop := normalize_pop(pop)]

CL <- fread(CLUSTER_FILE, sep = "\t", header = FALSE, fill = TRUE, strip.white = TRUE)
CL <- CL[nzchar(V1) & nzchar(V2), .(pop = V1, cluster = V2)]
CL <- CL[!(pop %in% c("pop", "POP", "sample"))]
CL[, pop := normalize_pop(pop)]
CL <- unique(CL, by = "pop")

DEL <- merge(DEL, CL, by = "pop", all.x = TRUE)

# =========================================================
# 4. filter genes
# =========================================================
D <- DEL[gene %in% target_genes]
D <- D[!is.na(deltaAF) & !is.na(pop) & !is.na(gene)]

# 只保留 table 6 对应 region 也行；先全保留
D[, gene := factor(gene, levels = target_genes)]
D[, pop_label := ifelse(is.na(cluster), pop, paste0(pop, " (", cluster, ")"))]

cat("Rows kept:", nrow(D), "\n")
print(D[, .N, by = .(region, gene)][order(region, gene)])

# =========================================================
# 5. summary
# =========================================================
SUM <- D[, .(
  mean_deltaAF   = mean(deltaAF, na.rm = TRUE),
  median_deltaAF = median(deltaAF, na.rm = TRUE),
  sd_deltaAF     = sd(deltaAF, na.rm = TRUE),
  n_snps         = .N,
  se_deltaAF     = sd(deltaAF, na.rm = TRUE) / sqrt(.N)
), by = .(region, gene, pop_label)]

fwrite(
  SUM[order(region, gene, mean_deltaAF)],
  file = file.path(OUTDIR, "Table6_top_genes_deltaAF_by_population_summary.tsv"),
  sep = "\t"
)

# =========================================================
# 6. function: one gene per plot
# =========================================================
make_gene_plot <- function(D_gene, SUM_gene, gene_name, outfile_prefix){
  
  # 排序 population
  ord <- SUM_gene[order(region, mean_deltaAF), pop_label]
  ord <- unique(ord)
  
  D_gene[, pop_label := factor(pop_label, levels = ord)]
  SUM_gene[, pop_label := factor(pop_label, levels = ord)]
  
  p <- ggplot(D_gene, aes(x = deltaAF, y = pop_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    
    geom_point(
      alpha = 0.15,
      size = 0.8,
      position = position_jitter(height = 0.12, width = 0)
    ) +
    
    # mean +/- SE as horizontal segment
    geom_segment(
      data = SUM_gene,
      aes(
        x = mean_deltaAF - se_deltaAF,
        xend = mean_deltaAF + se_deltaAF,
        y = pop_label,
        yend = pop_label
      ),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    
    geom_point(
      data = SUM_gene,
      aes(x = mean_deltaAF, y = pop_label),
      inherit.aes = FALSE,
      size = 2
    ) +
    
    facet_grid(region ~ ., scales = "free_y", space = "free_y") +
    
    labs(
      x = "deltaAF",
      y = "Population",
      title = paste0(gene_name, ": population-specific allele-frequency shifts"),
      subtitle = "Grey points = SNP-level deltaAF; large points and lines = mean +/- SE"
    ) +
    
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold")
    )
  
  print(p)
  
  ggsave(
    filename = file.path(OUTDIR, paste0(outfile_prefix, ".pdf")),
    plot = p, width = 7.5, height = 8.5, units = "in"
  )
  ggsave(
    filename = file.path(OUTDIR, paste0(outfile_prefix, ".png")),
    plot = p, width = 7.5, height = 8.5, units = "in", dpi = 300
  )
}

# =========================================================
# 7. single-gene plots
# =========================================================
for(g in target_genes){
  Dg <- D[gene == g]
  SUMg <- SUM[gene == g]
  
  cat("\nGene:", g, "\n")
  cat("  Dg rows =", nrow(Dg), "\n")
  cat("  SUMg rows =", nrow(SUMg), "\n")
  print(SUMg[order(region, mean_deltaAF)][1:min(6, .N)])
  
  if(nrow(Dg) == 0 || nrow(SUMg) == 0){
    cat("  Skipped: no data\n")
  } else {
    make_gene_plot(
      D_gene = Dg,
      SUM_gene = SUMg,
      gene_name = g,
      outfile_prefix = paste0("deltaAF_by_population_", g)
    )
  }
}

# =========================================================
# 8. summary plot
# =========================================================
ord_all <- SUM[order(region, gene, mean_deltaAF), pop_label]
ord_all <- unique(ord_all)

D[, pop_label := factor(pop_label, levels = ord_all)]
SUM[, pop_label := factor(pop_label, levels = ord_all)]

p_summary <- ggplot(D, aes(x = deltaAF, y = pop_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  
  geom_point(
    alpha = 0.12,
    size = 0.7,
    position = position_jitter(height = 0.12, width = 0)
  ) +
  
  geom_segment(
    data = SUM,
    aes(
      x = mean_deltaAF - se_deltaAF,
      xend = mean_deltaAF + se_deltaAF,
      y = pop_label,
      yend = pop_label
    ),
    inherit.aes = FALSE,
    linewidth = 0.55
  ) +
  
  geom_point(
    data = SUM,
    aes(x = mean_deltaAF, y = pop_label),
    inherit.aes = FALSE,
    size = 1.8
  ) +
  
  facet_grid(region ~ gene, scales = "free_y", space = "free_y") +
  
  labs(
    x = "deltaAF",
    y = "Population",
    title = "Population-specific allele-frequency shifts for top mitonuclear candidate genes",
    subtitle = "Grey points = SNP-level deltaAF; large points and lines = mean +/- SE"
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8.5),
    plot.title = element_text(face = "bold")
  )

print(p_summary)

ggsave(
  filename = file.path(OUTDIR, "deltaAF_by_population_summary_3genes.pdf"),
  plot = p_summary, width = 11, height = 8.5, units = "in"
)
ggsave(
  filename = file.path(OUTDIR, "deltaAF_by_population_summary_3genes.png"),
  plot = p_summary, width = 11, height = 8.5, units = "in", dpi = 300
)

cat("\nDone.\n")
cat("Output dir: ", OUTDIR, "\n")
print(list.files(OUTDIR))
#