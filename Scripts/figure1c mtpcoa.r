#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(grid)
  library(scales)
})

# ======================
# Input and output
# ======================
mt_tree <- paste0(
  "/work/cyu/poolseq/PPalign_output/overlap.vcf/",
  "consensus/mt_noDloop_noAMO.iqtree.treefile"
)

out_dir <- paste0(
  "/mnt/spareHD_2/nu_287/_assoc72_subunit/",
  "mitoPC_cluster_circle_noAMO"
)

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ======================
# Read mitochondrial tree
# ======================
tr <- read.tree(mt_tree)

if ("AMO" %in% tr$tip.label) {
  tr <- drop.tip(tr, "AMO")
}

# ======================
# PCoA
# ======================
D <- cophenetic(tr)

pco <- cmdscale(
  D,
  k = 5,
  eig = TRUE
)

expl <- 100 * pco$eig / sum(pco$eig)

PC <- data.frame(
  pop = rownames(pco$points),
  mitoPC1 = pco$points[, 1],
  mitoPC2 = pco$points[, 2],
  mitoPC3 = pco$points[, 3],
  mitoPC4 = pco$points[, 4],
  mitoPC5 = pco$points[, 5],
  stringsAsFactors = FALSE
)

# ======================
# Region and mt cluster
# ======================
PC <- PC %>%
  mutate(
    region = case_when(
      pop %in% c(
        "FG", "LG", "SR", "SL", "TL",
        "WB", "WT", "WK", "LB"
      ) ~ "AK",
      
      # PACH and FRED are BC populations
      pop %in% c(
        "SWA", "THE", "JOE", "BEA", "MUC",
        "PYE", "BOOT", "ECHO", "LAW", "GOS",
        "ROB", "PACH", "FRED"
      ) ~ "BC",
      
      pop %in% c(
        "RS", "SAY", "SC", "CH"
      ) ~ "Marine",
      
      TRUE ~ NA_character_
    ),
    
    mtCluster = case_when(
      pop %in% c(
        "SR", "TL", "WK", "LB", "MUC", "SWA"
      ) ~ "C1",
      
      pop %in% c(
        "PACH", "FRED", "BEA", "THE"
      ) ~ "C2",
      
      pop %in% c(
        "RS", "SAY", "SC", "CH",
        "ROB", "WB", "LG", "SL",
        "LAW", "BOOT", "JOE", "FG"
      ) ~ "C3",
      
      pop %in% c(
        "GOS", "PYE", "ECHO", "WT"
      ) ~ "C4",
      
      TRUE ~ NA_character_
    )
  )

# ======================
# Check assignments
# ======================
# Stop instead of creating an "Other" category
if (anyNA(PC$region)) {
  missing_region <- PC$pop[is.na(PC$region)]
  
  stop(
    "Region not assigned for: ",
    paste(missing_region, collapse = ", ")
  )
}

if (anyNA(PC$mtCluster)) {
  missing_cluster <- PC$pop[is.na(PC$mtCluster)]
  
  stop(
    "mtCluster not assigned for: ",
    paste(missing_cluster, collapse = ", ")
  )
}

PC <- PC %>%
  mutate(
    region = factor(
      region,
      levels = c("AK", "BC", "Marine")
    ),
    
    mtCluster = factor(
      mtCluster,
      levels = c("C1", "C2", "C3", "C4")
    )
  )

# ======================
# Colors
# Use the same colors as Panels A and B
# ======================
region_cols <- c(
  "AK" = "#2B2BFF",
  "BC" = "#7EE69A",
  "Marine" = "#FF2D2D"
)

cluster_cols <- c(
  "C1" = "#2B2BFF",
  "C2" = "#FFA04D",
  "C3" = "#FF2D2D",
  "C4" = "#7EE69A"
)

# ======================
# Panel C
# ======================
p <- ggplot(
  PC,
  aes(
    x = mitoPC1,
    y = mitoPC2
  )
) +
  
  # Ellipses are calculated automatically from the points
  # in each mitochondrial cluster.
  geom_mark_ellipse(
    aes(
      group = mtCluster,
      colour = mtCluster
    ),
    fill = NA,
    linetype = 2,
    linewidth = 0.8,
    
    # Visual padding around the points
    expand = grid::unit(10, "mm"),
    
    alpha = 1,
    show.legend = TRUE
  ) +
  
  # Points are colored by geographic region
  geom_point(
    aes(fill = region),
    shape = 21,
    size = 4,
    colour = "black",
    stroke = 0.4,
    show.legend = TRUE
  ) +
  
  # Labels may be outside the ellipses.
  # Connecting lines show which point each label represents.
  geom_text_repel(
    aes(label = pop),
    size = 3.4,
    fontface = "bold",
    seed = 1,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.colour = "black",
    segment.linewidth = 0.3,
    show.legend = FALSE
  ) +
  
  # Region legend: points only, without dashed boxes
  scale_fill_manual(
    name = "Region",
    values = region_cols,
    breaks = c("AK", "BC", "Marine"),
    drop = TRUE,
    guide = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        size = 4,
        colour = "black",
        linetype = 0,
        linewidth = 0
      )
    )
  ) +
  
  # mt cluster legend: dashed lines only
  scale_colour_manual(
    name = "mt cluster",
    values = cluster_cols,
    breaks = c("C1", "C2", "C3", "C4"),
    drop = TRUE,
    guide = guide_legend(
      order = 2,
      override.aes = list(
        shape = NA,
        fill = NA,
        linetype = 2,
        linewidth = 0.8
      )
    )
  ) +
  
  labs(
    x = sprintf(
      "mitoPC1 (%.1f%%)",
      expl[1]
    ),
    y = sprintf(
      "mitoPC2 (%.1f%%)",
      expl[2]
    )
  ) +
  
  # Labels are allowed to extend outside the plotting panel
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.9
    ),
    
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(
      colour = "black",
      size = 14
    ),
    
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    
    # Put legends outside on the right
    legend.position = "right",
    
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    legend.key = element_rect(
      fill = NA,
      colour = NA
    ),
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.spacing.y = grid::unit(
      0.05,
      "cm"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )
# ======================
# Output files
# ======================
png_file <- file.path(
  out_dir,
  "mitoPC_noAMO_cluster_circle.png"
)

pdf_file <- file.path(
  out_dir,
  "mitoPC_noAMO_cluster_circle.pdf"
)

coordinate_file <- file.path(
  out_dir,
  "mitoPC_noAMO_cluster_coordinates.tsv"
)

explained_file <- file.path(
  out_dir,
  "mitoPC_noAMO_explained.tsv"
)

# ======================
# Save plot
# ======================
ggsave(
  filename = png_file,
  plot = p,
  width = 7.2,
  height = 5.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 7.2,
  height = 5.8,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ======================
# Save coordinates
# ======================
fwrite(
  PC,
  coordinate_file,
  sep = "\t"
)

# ======================
# Save explained variation
# ======================
fwrite(
  data.table(
    PC = paste0("PC", 1:5),
    explained_percent = expl[1:5]
  ),
  explained_file,
  sep = "\t"
)

# ======================
# Print results
# ======================
cat(
  "[OK] Output directory:",
  out_dir,
  "\n"
)

cat(
  sprintf(
    "mitoPC1 = %.2f%%; mitoPC2 = %.2f%%\n",
    expl[1],
    expl[2]
  )
)

cat(
  "[OK] PNG:",
  png_file,
  "\n"
)

cat(
  "[OK] PDF:",
  pdf_file,
  "\n"
)

print(p)