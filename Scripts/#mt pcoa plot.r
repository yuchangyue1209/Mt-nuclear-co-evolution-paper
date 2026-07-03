#mt pcoa plot
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# ======================
# input / output
# ======================
mt_tree <- "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_noDloop_noAMO.iqtree.treefile"
out_dir <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_cluster_circle_noAMO"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ======================
# read tree
# ======================
tr <- read.tree(mt_tree)

if ("AMO" %in% tr$tip.label) {
  tr <- drop.tip(tr, "AMO")
}

# ======================
# PCoA
# ======================
D <- cophenetic(tr)
pco <- cmdscale(D, k = 5, eig = TRUE)

expl <- 100 * pco$eig / sum(pco$eig)

PC <- data.frame(
  pop = rownames(pco$points),
  mitoPC1 = pco$points[,1],
  mitoPC2 = pco$points[,2],
  mitoPC3 = pco$points[,3],
  mitoPC4 = pco$points[,4],
  mitoPC5 = pco$points[,5]
)

# ======================
# region + mt cluster label
# ======================
PC <- PC %>%
  mutate(
    region = case_when(
      pop %in% c("FG","LG","SR","SL","TL","WB","WT","WK","LB") ~ "AK",
      pop %in% c("SWA","THE","JOE","BEA","MUC","PYE","ROS","BOOT","ECHO","LAW","GOS","ROB") ~ "BC",
      pop %in% c("RS","SAY") ~ "Marine",
      pop %in% c("PACH","FRED","SC","CH") ~ "Recent",
      TRUE ~ "Other"
    ),
    mtCluster = case_when(
      pop %in% c("SR","TL","WK","LB","MUC","SWA") ~ "C1",
      pop %in% c("PACH","FRED","BEA","THE") ~ "C2",
      pop %in% c("RS","SAY","SC","CH","ROB","WB","LG","SL","LAW","BOOT","JOE","FG") ~ "C3",
      pop %in% c("GOS","PYE","ECHO","WT") ~ "C4",
      TRUE ~ "Other"
    )
  )

# ======================
# colors
# ======================
region_cols <- c(
  "BC" = "#7EE69A",
  "AK" = "#2B2BFF",
  "Marine" = "#FF1F1F",
  "Recent" = "#FFA04D",
  "Other" = "grey60"
)

cluster_cols <- c(
  "C1" = "#2B2BFF",
  "C2" = "#FFA04D",
  "C3" = "#FF1F1F",
  "C4" = "#7EE69A",
  "Other" = "grey60"
)

# ======================
# plot
# ======================
p <- ggplot(PC, aes(mitoPC1, mitoPC2)) +
  stat_ellipse(
    aes(group = mtCluster, color = mtCluster),
    type = "norm",
    linetype = 2,
    linewidth = 0.8,
    alpha = 0.8,
    show.legend = TRUE
  ) +
  geom_point(aes(fill = region), shape = 21, size = 3.8, color = "black", stroke = 0.4) +
  geom_text_repel(
    aes(label = pop),
    size = 3.4,
    fontface = "bold",
    max.overlaps = Inf,
    seed = 1,
    box.padding = 0.35,
    point.padding = 0.25
  ) +
  scale_fill_manual(values = region_cols, name = "Region") +
  scale_color_manual(values = cluster_cols, name = "mt cluster") +
  labs(
    x = sprintf("mitoPC1 (%.1f%%)", expl[1]),
    y = sprintf("mitoPC2 (%.1f%%)", expl[2]),
    title = "Mitochondrial PCoA after removing AMO"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# ======================
# save
# ======================
ggsave(file.path(out_dir, "mitoPC_noAMO_cluster_circle.png"),
       p, width = 7.8, height = 5.6, dpi = 300)

ggsave(file.path(out_dir, "mitoPC_noAMO_cluster_circle.pdf"),
       p, width = 7.8, height = 5.6)

fwrite(PC, file.path(out_dir, "mitoPC_noAMO_cluster_coordinates.tsv"), sep = "\t")

fwrite(
  data.table(
    PC = paste0("PC", 1:5),
    explained_percent = expl[1:5]
  ),
  file.path(out_dir, "mitoPC_noAMO_explained.tsv"),
  sep = "\t"
)

cat("[OK] written to:", out_dir, "\n")
cat(sprintf("PC1 = %.2f%%; PC2 = %.2f%%\n", expl[1], expl[2]))
print(p)





#2026-0629 version
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(grid)
})

#======================
# input / output
#======================
mt_tree <- "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/mt_noDloop_noAMO.iqtree.treefile"
out_dir <- "/mnt/spareHD_2/nu_287/_assoc72_subunit/mitoPC_cluster_circle_noAMO"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#======================
# read tree
#======================
tr <- read.tree(mt_tree)

if ("AMO" %in% tr$tip.label)
  tr <- drop.tip(tr, "AMO")

#======================
# PCoA
#======================
D <- cophenetic(tr)
pco <- cmdscale(D, k = 5, eig = TRUE)

expl <- 100 * pco$eig / sum(pco$eig)

PC <- data.frame(
  pop = rownames(pco$points),
  mitoPC1 = pco$points[,1],
  mitoPC2 = pco$points[,2],
  mitoPC3 = pco$points[,3],
  mitoPC4 = pco$points[,4],
  mitoPC5 = pco$points[,5]
)

#======================
# Region
#======================
PC <- PC %>%
  mutate(
    region = case_when(
      pop %in% c("FG","LG","SR","SL","TL","WB","WT","WK","LB") ~ "AK",
      pop %in% c("SWA","THE","JOE","BEA","MUC","PYE","BOOT","ECHO","LAW","GOS","ROB") ~ "BC",
      pop %in% c("RS","SAY","PACH","FRED","SC","CH") ~ "Marine",
      TRUE ~ "Other"
    ),
    mtCluster = case_when(
      pop %in% c("SR","TL","WK","LB","MUC","SWA") ~ "C1",
      pop %in% c("PACH","FRED","BEA","THE") ~ "C2",
      pop %in% c("RS","SAY","SC","CH","ROB","WB","LG","SL","LAW","BOOT","JOE","FG") ~ "C3",
      pop %in% c("GOS","PYE","ECHO","WT") ~ "C4",
      TRUE ~ "Other"
    )
  )

#======================
# colors
#======================
region_cols <- c(
  AK = "#2B2BFF",
  BC = "#7EE69A",
  Marine = "#FF2D2D",
  Other = "grey60"
)

cluster_cols <- c(
  C1 = "#2B2BFF",
  C2 = "#FFA04D",
  C3 = "#FF2D2D",
  C4 = "#7EE69A"
)

#======================
# plot
#======================
p <- ggplot(PC, aes(mitoPC1, mitoPC2)) +

  geom_mark_ellipse(
    aes(group = mtCluster, colour = mtCluster),
    fill = NA,
    linewidth = 0.9,
    linetype = 2,
    expand = unit(3, "mm"),
    alpha = 1,
    show.legend = TRUE
  ) +

  geom_point(
    aes(fill = region),
    shape = 21,
    size = 4,
    colour = "black",
    stroke = 0.4
  ) +

  geom_text_repel(
    aes(label = pop),
    size = 3.4,
    fontface = "bold",
    seed = 1,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25
  ) +

  scale_fill_manual(values = region_cols, name = "Region") +
  scale_colour_manual(values = cluster_cols, name = "mt cluster") +

  labs(
    x = sprintf("mitoPC1 (%.1f%%)", expl[1]),
    y = sprintf("mitoPC2 (%.1f%%)", expl[2])
  ) +

  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
    panel.grid.major = element_line(colour = "grey88"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

#======================
# save
#======================
ggsave(
  file.path(out_dir, "mitoPC_noAMO_cluster_circle.png"),
  p,
  width = 8,
  height = 5.8,
  dpi = 300
)

ggsave(
  file.path(out_dir, "mitoPC_noAMO_cluster_circle.pdf"),
  p,
  width = 8,
  height = 5.8
)

fwrite(
  PC,
  file.path(out_dir, "mitoPC_noAMO_cluster_coordinates.tsv"),
  sep = "\t"
)

fwrite(
  data.table(
    PC = paste0("PC", 1:5),
    explained_percent = expl[1:5]
  ),
  file.path(out_dir, "mitoPC_noAMO_explained.tsv"),
  sep = "\t"
)

cat("[OK] written to:", out_dir, "\n")
cat(sprintf("PC1 = %.2f%%; PC2 = %.2f%%\n", expl[1], expl[2]))

print(p)
