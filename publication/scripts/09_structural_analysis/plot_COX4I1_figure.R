#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# COX4I1 replacement of the previous NDUFS2-focused Figure 7.
# The project reference protein carries the freshwater-associated states
# V58/Y124/L133. The population table reports frequencies of the alternative
# marine-associated amino acids A58/F124/F133, so freshwater-state frequency
# is calculated as 1 - altAA frequency.

freq_file <- Sys.getenv(
  "COX4I1_FREQ_FILE",
  "/path/to/workspace/OXPHOS_structural_validation/287_recalled_nonsyn/nonsynonymous_variant_population_frequencies.tsv"
)
protein_file <- Sys.getenv(
  "COX4I1_PROTEIN_FILE",
  "/path/to/workspace/cox4i1_coordinate_validation/cox4i1.reference_protein.fa"
)
outdir <- Sys.getenv(
  "COX4I1_FIGURE_DIR",
  "/path/to/workspace/OXPHOS_structural_validation/cox4i1_Figure7"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(freq_file), file.exists(protein_file))

pop_order <- c(
  "RS", "FG", "LG", "SL", "SR", "TL", "WB", "WK", "WT",
  "SAY", "AMO", "BEA", "BOOT", "ECHO", "GOS", "JOE", "LAW",
  "MUC", "PYE", "ROB", "SWA", "THE"
)

variant_key <- data.table(
  genomic_position = c(12508430L, 12508082L, 12507940L),
  genomic_ref = c("A", "T", "G"),
  genomic_alt = c("G", "A", "A"),
  codon_position = c(58L, 124L, 133L),
  biological_variant = c("A58V", "F124Y", "F133L"),
  freshwater_AA = c("V", "Y", "L"),
  panel_title = c(
    "A58V: frequency of V58",
    "F124Y: frequency of Y124",
    "F133L: frequency of L133"
  ),
  variant_color = c("#D95F8D", "#6A3D9A", "#C23B3B")
)

raw <- fread(freq_file)
setnames(raw, names(raw), trimws(names(raw)))

x <- merge(
  raw[
    tolower(gene) == "cox4i1" &
      position %in% variant_key$genomic_position
  ],
  variant_key,
  by.x = c("position", "genomic_ref", "genomic_alt", "codon_position"),
  by.y = c("genomic_position", "genomic_ref", "genomic_alt", "codon_position"),
  all = FALSE
)

if (nrow(x) != 60L) {
  stop("Expected 60 COX4I1 freshwater comparison rows (3 variants x 20 populations); observed ", nrow(x))
}

# Add each marine population once per variant.
marine <- unique(x[, .(
  biological_variant,
  panel_title,
  variant_color,
  population = marine_population,
  region,
  habitat = "Marine",
  freshwater_state_frequency = 1 - marine_altAA_frequency
)])

freshwater <- x[, .(
  biological_variant,
  panel_title,
  variant_color,
  population = freshwater_population,
  region,
  habitat = "Freshwater",
  freshwater_state_frequency = 1 - freshwater_altAA_frequency
)]

plot_dt <- rbindlist(list(marine, freshwater), use.names = TRUE)
plot_dt[, population := factor(population, levels = pop_order)]
plot_dt[, panel_title := factor(panel_title, levels = variant_key$panel_title)]
plot_dt[, group := fifelse(
  habitat == "Marine", "Marine",
  fifelse(region == "AK", "Alaska freshwater", "British Columbia freshwater")
)]

qa <- plot_dt[, .(
  total_populations = uniqueN(population),
  marine_populations = uniqueN(population[habitat == "Marine"]),
  freshwater_populations = uniqueN(population[habitat == "Freshwater"])
), by = biological_variant]
print(qa)
if (any(qa$total_populations != 22L) ||
    any(qa$marine_populations != 2L) ||
    any(qa$freshwater_populations != 20L)) {
  stop("Population audit failed for COX4I1 Panel A")
}

fwrite(
  plot_dt[order(panel_title, population)],
  file.path(outdir, "COX4I1_population_freshwater_state_frequencies.tsv"),
  sep = "\t"
)

group_cols <- c(
  "Marine" = "#303030",
  "Alaska freshwater" = "#1594C5",
  "British Columbia freshwater" = "#00A878"
)
group_shapes <- c(
  "Marine" = 17,
  "Alaska freshwater" = 16,
  "British Columbia freshwater" = 16
)

theme_fig <- theme_bw(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.55),
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    strip.background = element_rect(fill = "grey94", color = "black", linewidth = 0.55),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.margin = margin(5, 10, 5, 5)
  )

pA <- ggplot(
  plot_dt,
  aes(x = population, y = freshwater_state_frequency, color = group, shape = group)
) +
  geom_hline(yintercept = c(0, 1), color = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 9.5, linetype = "dashed", color = "grey45", linewidth = 0.55) +
  geom_point(size = 3.6, stroke = 0.25) +
  facet_wrap(~panel_title, ncol = 1, scales = "fixed", strip.position = "top") +
  scale_color_manual(values = group_cols) +
  scale_shape_manual(values = group_shapes) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Frequency of freshwater-associated amino acid"
  ) +
  theme_fig +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    legend.text = element_text(size = 12)
  )

# ---------- Panel B: sequence windows ----------
read_fasta_sequence <- function(path) {
  z <- readLines(path, warn = FALSE)
  paste0(z[!grepl("^>", z)], collapse = "")
}

fresh_seq <- strsplit(read_fasta_sequence(protein_file), "", fixed = TRUE)[[1]]
if (length(fresh_seq) != 169L) stop("Expected a 169-aa COX4I1 sequence; observed ", length(fresh_seq))
if (!identical(fresh_seq[c(58, 124, 133)], c("V", "Y", "L"))) {
  stop("Reference COX4I1 does not carry the validated freshwater states V58/Y124/L133")
}
marine_seq <- fresh_seq
marine_seq[c(58, 124, 133)] <- c("A", "F", "F")

aa_class <- function(aa) {
  ifelse(aa %in% c("D", "E"), "acidic",
    ifelse(aa %in% c("K", "R", "H"), "basic",
      ifelse(aa %in% c("S", "T", "N", "Q", "C"), "polar", "hydrophobic")
    )
  )
}
aa_cols <- c(
  acidic = "#2F71C1",
  basic = "#C43AA9",
  polar = "#159447",
  hydrophobic = "#E33D3D"
)

make_sequence_window <- function(start, end, title) {
  pos <- start:end
  d <- rbindlist(list(
    data.table(state = "Marine A58–F124–F133", position = pos, aa = marine_seq[pos], y = 2),
    data.table(state = "Freshwater V58–Y124–L133", position = pos, aa = fresh_seq[pos], y = 1)
  ))
  d[, aa_class := aa_class(aa)]
  highlights <- intersect(pos, c(58L, 124L, 133L))
  # Draw one small yellow cell around each sequence character, matching the
  # NDUFS2 figure. Do not join the two sequence rows into one tall rectangle.
  highlight_cells <- CJ(position = highlights, y = c(1, 2))

  ggplot(d, aes(x = position, y = y)) +
    geom_rect(
      data = highlight_cells,
      aes(
        xmin = position - 0.45,
        xmax = position + 0.45,
        ymin = y - 0.38,
        ymax = y + 0.38
      ),
      inherit.aes = FALSE,
      fill = "#FFF2A8",
      color = "#D9A500",
      linewidth = 0.55
    ) +
    geom_text(aes(label = aa, color = aa_class), family = "mono", fontface = "bold", size = 4.5) +
    scale_color_manual(values = aa_cols, guide = "none") +
    scale_y_continuous(
      breaks = c(2, 1),
      labels = c(
        "Marine A58–F124–F133",
        "Freshwater V58–Y124–L133"
      ),
      limits = c(0.55, 2.45)
    ) +
    scale_x_continuous(breaks = sort(unique(c(start, end, highlights)))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = "black", face = "bold", size = 10),
      axis.text.x = element_text(color = "grey35", size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.margin = margin(4, 6, 4, 6)
    )
}

pB1 <- make_sequence_window(50, 66, "COX4I1 residues 50–66")
pB2 <- make_sequence_window(116, 141, "COX4I1 residues 116–141")
pB_inner <- pB1 / pB2 + plot_layout(heights = c(1, 1))
# Treat the two sequence windows as one figure panel so that patchwork assigns
# one shared panel tag (B), as in the NDUFS2 version.
pB <- wrap_elements(full = pB_inner)

# ---------- Panel C: protein map and interface classification ----------
site_dt <- copy(variant_key)
site_dt[, interface_class := c("Distant from mt interface", "Near mt interface", "Direct mt interface")]
site_dt[, distance_label := c("21–23 Å", "~7.6 Å", "3.2–3.9 Å")]
site_dt[, label_y := c(1.00, 1.00, 1.23)]
site_dt[, label_hjust := c(0.5, 1.0, 0.0)]

pC <- ggplot() +
  annotate("rect", xmin = 1, xmax = 169, ymin = 0.42, ymax = 0.68,
           fill = "#9BC8E5", color = "black", linewidth = 0.55) +
  annotate("text", x = 85, y = 0.55, label = "COX4I1 (169 aa)",
           size = 4.8, fontface = "bold") +
  geom_segment(
    data = site_dt,
    aes(x = codon_position, xend = codon_position, y = 0.68, yend = 0.95, color = biological_variant),
    linewidth = 0.8
  ) +
  geom_point(
    data = site_dt,
    aes(x = codon_position, y = 0.95, color = biological_variant),
    size = 4
  ) +
  geom_text(
    data = site_dt,
    aes(
      x = codon_position,
      y = label_y,
      label = biological_variant,
      color = biological_variant
    ),
    hjust = site_dt$label_hjust,
    size = 4.0,
    fontface = "bold",
    lineheight = 0.9
  ) +
  scale_color_manual(values = setNames(variant_key$variant_color, variant_key$biological_variant)) +
  scale_x_continuous(limits = c(1, 169), breaks = c(1, 58, 100, 124, 133, 169)) +
  scale_y_continuous(limits = c(0.30, 1.40), expand = c(0, 0)) +
  labs(x = "Amino-acid position", y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold"),
    legend.position = "none",
    plot.margin = margin(12, 10, 4, 10)
  )

combined <- pA / pB / pC +
  plot_layout(heights = c(2.9, 1.55, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 22))

pdf_file <- file.path(outdir, "Figure7_COX4I1_ABC.pdf")
png_file <- file.path(outdir, "Figure7_COX4I1_ABC.png")

ggsave(pdf_file, combined, width = 12.5, height = 14.2, device = cairo_pdf)
ggsave(png_file, combined, width = 12.5, height = 14.2, dpi = 400, bg = "white")

cat("[OK] COX4I1 Figure 7 generated\n")
cat("[PDF] ", pdf_file, "\n", sep = "")
cat("[PNG] ", png_file, "\n", sep = "")
