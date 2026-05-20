# Publication-style map
# 27 stickleback populations
# =========================

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(patchwork)
library(dplyr)

# -------------------------
# 1. Input data
# -------------------------
pop <- read.table(text = "
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
LB Alaska 100 61.559183 -149.258019 Freshwater 'Mat-Su valley'
CH Alaska 100 61.2023128 -149.761914 Recent_Colonized Anchorage
", header = TRUE)

pop$Habitat <- factor(
  pop$Habitat,
  levels = c("Marine", "Freshwater", "Recent_Colonized")
)

# -------------------------
# 2. Optional mtCluster
# 如果你不想按 mtCluster 上色，可以跳过这一段
# -------------------------
pop <- pop %>%
  mutate(
    mtCluster = case_when(
      Population %in% c("SR", "TL", "WK", "LB", "MUC", "SWA") ~ "C1_AK",
      Population %in% c("PACH", "FRED", "BEA", "THE") ~ "C2_Recent",
      Population %in% c("RS", "SAY", "SC", "CH", "ROB", "WB", "LG", "SL",
                        "LAW", "BOOT", "JOE", "FG") ~ "C3_MarineLike",
      Population %in% c("GOS", "PYE", "ECHO", "WT") ~ "C4_GOS",
      Population == "AMO" ~ "AMO",
      TRUE ~ "Other"
    )
  )

# -------------------------
# 3. Map background
# -------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

theme_map_pub <- theme_void(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = c(0.52, 0.52),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

cluster_cols <- c(
  "C1_AK" = "#4E79A7",
  "C2_Recent" = "#59A14F",
  "C3_MarineLike" = "#F28E2B",
  "C4_GOS" = "#B07AA1",
  "AMO" = "#E15759",
  "Other" = "grey40"
)

hab_shapes <- c(
  "Marine" = 24,
  "Freshwater" = 21,
  "Recent_Colonized" = 22
)

# -------------------------
# 4. Alaska inset
# -------------------------
map_ak <- ggplot() +
  geom_sf(data = world, fill = "grey92", color = "grey70", linewidth = 0.25) +
  geom_point(
    data = filter(pop, Region == "Alaska"),
    aes(x = Longitude, y = Latitude, fill = mtCluster, shape = Habitat),
    size = 3.8,
    color = "black",
    stroke = 0.35
  ) +
  geom_text_repel(
    data = filter(pop, Region == "Alaska"),
    aes(x = Longitude, y = Latitude, label = Population),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey50",
    segment.linewidth = 0.25
  ) +
  coord_sf(
    xlim = c(-152.2, -148.5),
    ylim = c(60.2, 62.0),
    expand = FALSE
  ) +
  scale_fill_manual(values = cluster_cols) +
  scale_shape_manual(values = hab_shapes) +
  labs(title = "Alaska", fill = "mtDNA cluster", shape = "Habitat") +
  theme_map_pub

# -------------------------
# 5. British Columbia inset
# -------------------------
map_bc <- ggplot() +
  geom_sf(data = world, fill = "grey92", color = "grey70", linewidth = 0.25) +
  geom_point(
    data = filter(pop, Region == "BC"),
    aes(x = Longitude, y = Latitude, fill = mtCluster, shape = Habitat),
    size = 3.8,
    color = "black",
    stroke = 0.35
  ) +
  geom_text_repel(
    data = filter(pop, Region == "BC"),
    aes(x = Longitude, y = Latitude, label = Population),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey50",
    segment.linewidth = 0.25
  ) +
  coord_sf(
    xlim = c(-128.8, -124.6),
    ylim = c(48.5, 51.0),
    expand = FALSE
  ) +
  scale_fill_manual(values = cluster_cols) +
  scale_shape_manual(values = hab_shapes) +
  labs(title = "British Columbia", fill = "mtDNA cluster", shape = "Habitat") +
  theme_map_pub +
  theme(legend.position = "none")
# -------------------------
# 6. Combine figure
# -------------------------
final_map <- map_ak + map_bc +
  plot_layout(ncol = 2) &
  theme(
    plot.title = element_text(size = 15, face = "bold")
  )

print(final_map)

# -------------------------
# 7. Save output
# -------------------------
ggsave(
  filename = "stickleback_sampling_map_NEE_style.pdf",
  plot = final_map,
  width = 11,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  filename = "stickleback_sampling_map_NEE_style.png",
  plot = final_map,
  width = 11,
  height = 5.5,
  dpi = 600
)