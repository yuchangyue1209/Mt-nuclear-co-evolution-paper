#!/usr/bin/env Rscript

# Figure S1: sampling locations for 27 stickleback populations.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggspatial)
  library(patchwork)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(sf)
})

# Set OUTPUT_DIR when running outside the repository root.
out_dir <- Sys.getenv(
  "OUTPUT_DIR",
  unset = file.path("results", "FigureS1_sampling_map")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pop <- read.table(
  text = "
Population Region n_individuals Latitude Longitude Habitat Watershed
FG Alaska 137 61.6056 -149.2792 Freshwater 'Mat-Su valley'
LG Alaska 138 61.577305 -149.767506 Freshwater 'Mat-Su valley'
SR Alaska 94 61.667545 -150.135689 Freshwater 'Mat-Su valley'
SL Alaska 74 60.596104 -150.997664 Freshwater 'Kenai peninsula'
TL Alaska 56 60.533128 -149.547014 Freshwater 'Kenai peninsula'
WB Alaska 96 61.62 -149.213 Freshwater 'Mat-Su valley'
WT Alaska 92 60.539 -150.467 Freshwater 'Kenai peninsula'
WK Alaska 99 60.72 -151.25 Freshwater 'Kenai peninsula'
SWA BC 72 50.661395 -128.193026 Freshwater 'Northern Vancouver Island'
THE BC 78 50.517955 -126.983808 Freshwater 'Nimkish region'
JOE BC 77 50.623841 -127.483829 Freshwater 'Northern Vancouver Island'
BEA BC 55 50.60117 -127.310019 Freshwater 'Northern Vancouver Island'
MUC BC 58 49.874151 -126.164685 Freshwater 'Gold River'
PYE BC 29 50.297904 -125.584358 Freshwater 'Pye River'
AMO BC 61 50.158382 -125.580664 Freshwater 'Amor River'
SAY BC 100 50.377329 -125.950915 Marine 'Salmon River'
GOS BC 200 50.068186 -125.505989 Freshwater 'Campbell River'
ROB BC 200 50.21606 -125.54475 Freshwater 'Amor River'
BOOT BC 48 50.045048 -125.525494 Freshwater 'Campbell River'
ECHO BC 48 49.984894 -125.414318 Freshwater 'Campbell River'
FRED BC 46 48.857105 -125.02186 Recent_Colonized 'Bamfield peninsula'
LAW BC 48 50.038143 -125.577298 Freshwater 'Campbell River'
PACH BC 53 48.838426 -125.023866 Recent_Colonized 'Bamfield peninsula'
RS Alaska 200 61.5337007 -149.266752 Marine 'Mat-Su valley'
SC Alaska 100 60.535203 -150.831276 Recent_Colonized 'Mat-Su valley'
LB Alaska 100 61.559183 -149.258019 Recent_Colonized 'Mat-Su valley'
CH Alaska 100 61.2023128 -149.761914 Recent_Colonized Anchorage
",
  header = TRUE
)

pop$Habitat <- factor(
  pop$Habitat,
  levels = c("Marine", "Recent_Colonized", "Freshwater")
)

# Draw freshwater populations first and marine references last so that marine
# symbols remain visible where sampling locations are close.
pop_plot <- pop %>%
  mutate(
    draw_order = match(
      as.character(Habitat),
      c("Freshwater", "Recent_Colonized", "Marine")
    )
  ) %>%
  arrange(draw_order)

world <- ne_countries(scale = "medium", returnclass = "sf")

theme_map_pub <- theme_void(base_size = 15) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 18, face = "bold", hjust = 0),
    plot.margin = margin(5, 5, 5, 5),
    aspect.ratio = 1,
    panel.background = element_rect(fill = "#EAF4F7", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

habitat_cols <- c(
  Marine = "#F28E2B",
  Recent_Colonized = "#2CA02C",
  Freshwater = "#2B8CBE"
)

habitat_shapes <- c(
  Marine = 24,
  Recent_Colonized = 22,
  Freshwater = 21
)

habitat_labels <- c(
  Marine = "Marine",
  Recent_Colonized = "Recent",
  Freshwater = "Freshwater"
)

add_population_points <- function(region_name) {
  geom_point(
    data = filter(pop_plot, Region == region_name),
    aes(x = Longitude, y = Latitude, fill = Habitat, shape = Habitat),
    size = 3.8,
    color = "black",
    stroke = 0.4
  )
}

add_population_labels <- function(region_name) {
  geom_text_repel(
    data = filter(pop, Region == region_name),
    aes(x = Longitude, y = Latitude, label = Population),
    seed = 2026,
    size = 4.0,
    max.overlaps = Inf,
    max.time = 5,
    max.iter = 50000,
    force = 2,
    force_pull = 0.15,
    box.padding = 0.65,
    point.padding = 0.45,
    min.segment.length = 0,
    segment.color = "grey45",
    segment.linewidth = 0.25
  )
}

add_habitat_scales <- function(plot) {
  plot +
    scale_fill_manual(
      values = habitat_cols,
      breaks = c("Marine", "Recent_Colonized", "Freshwater"),
      labels = habitat_labels
    ) +
    scale_shape_manual(
      values = habitat_shapes,
      breaks = c("Marine", "Recent_Colonized", "Freshwater"),
      labels = habitat_labels
    )
}

map_ak <- ggplot() +
  geom_sf(
    data = world,
    fill = "#E8E8E8",
    color = "#B8B8B8",
    linewidth = 0.25
  ) +
  add_population_points("Alaska") +
  add_population_labels("Alaska") +
  coord_sf(
    xlim = c(-152.2, -148.5),
    ylim = c(60.2, 62.0),
    expand = FALSE
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    unit_category = "metric",
    style = "bar",
    text_cex = 0.90,
    line_width = 0.6,
    pad_x = grid::unit(0.15, "in"),
    pad_y = grid::unit(0.15, "in")
  ) +
  labs(title = "Alaska", fill = NULL, shape = NULL) +
  theme_map_pub

map_ak <- add_habitat_scales(map_ak)

map_bc <- ggplot() +
  geom_sf(
    data = world,
    fill = "#E8E8E8",
    color = "#B8B8B8",
    linewidth = 0.25
  ) +
  add_population_points("BC") +
  add_population_labels("BC") +
  coord_sf(
    xlim = c(-128.8, -124.6),
    ylim = c(48.5, 51.0),
    expand = FALSE
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    unit_category = "metric",
    style = "bar",
    text_cex = 0.90,
    line_width = 0.6,
    pad_x = grid::unit(0.15, "in"),
    pad_y = grid::unit(0.15, "in")
  ) +
  labs(title = "British Columbia", fill = NULL, shape = NULL) +
  theme_map_pub

map_bc <- add_habitat_scales(map_bc)

final_map <- map_ak + map_bc +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.justification = "center",
    legend.key = element_blank(),
    legend.spacing.y = grid::unit(0.08, "in")
  )

print(final_map)

ggsave(
  filename = file.path(out_dir, "stickleback_sampling_map_NEE_style.pdf"),
  plot = final_map,
  width = 12,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "stickleback_sampling_map_NEE_style.png"),
  plot = final_map,
  width = 12,
  height = 5.5,
  dpi = 600
)

message("[OK] PDF: ", file.path(out_dir, "stickleback_sampling_map_NEE_style.pdf"))
message("[OK] PNG: ", file.path(out_dir, "stickleback_sampling_map_NEE_style.png"))
