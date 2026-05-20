#fig1 panel e within diversity pi 
library(tidyverse)
library(rstatix)

# ============================================================
# 0. Input / output
# ============================================================

input_file <- "/work/cyu/Poolinfo_with_13gene_details.csv"

out_dir <- "/work/cyu/within_pi_complex_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Read data
# ============================================================

df <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE
)

# 修正 CH habitat
df$Habitat[df$Population == "CH"] <- "Recent Colonized"

# ============================================================
# 2. Convert gene-level pi to long format
# ============================================================

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

# ============================================================
# 3. Set population order
#    Alaska: marine -> recent -> freshwater
#    BC: marine -> recent -> freshwater
# ============================================================

pop_order <- c(
  # Alaska
  "RS",
  "SC", "CH",
  "FG", "LG", "SR", "SL", "TL", "WB", "WT", "WK", "LB",
  
  # BC
  "SAY",
  "FRED", "PACH",
  "SWA", "THE", "JOE", "BEA", "MUC", "PYE", "AMO",
  "GOS", "ROB", "BOOT", "ECHO", "LAW"
)

gene_pi_long <- gene_pi_long %>%
  dplyr::mutate(
    Region = factor(Region, levels = c("Alaska", "BC")),
    Habitat = factor(Habitat, levels = c("Marine", "Recent Colonized", "Freshwater")),
    Population = factor(Population, levels = pop_order),
    Complex = factor(
      Complex,
      levels = c("Complex I", "Complex III", "Complex IV", "Complex V")
    )
  )

# 保存整理后的长表
write.csv(
  gene_pi_long,
  file.path(out_dir, "gene_level_pi_long.csv"),
  row.names = FALSE
)

# ============================================================
# 4. Plot: within-population pi per complex
# ============================================================

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
    title = "Nucleotide Diversity per Complex across Populations",
    x = "Population",
    y = expression(pi),
    fill = "Complex"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 17, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

print(p)

ggsave(
  file.path(out_dir, "within_population_pi_complex_flat_2to1.pdf"),
  p,
  width = 12,
  height = 6
)

ggsave(
  file.path(out_dir, "within_population_pi_complex_flat_2to1.png"),
  p,
  width = 12,
  height = 6,
  dpi = 300
)

# ============================================================
# 5. Table S1: mean pi per population x complex
# ============================================================

table1 <- gene_pi_long %>%
  group_by(Population, Region, Habitat, Complex) %>%
  summarise(
    mean_pi = mean(pi, na.rm = TRUE),
    sd_pi = sd(pi, na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  arrange(
    Region,
    Habitat,
    Population,
    Complex
  ) %>%
  mutate(
    mean_pi = round(mean_pi, 5),
    sd_pi = round(sd_pi, 5)
  )

print(table1)

write.csv(
  table1,
  file.path(out_dir, "Table_S1_mean_pi_population_complex_long.csv"),
  row.names = FALSE
)

# Wide version for supplement readability
table1_wide <- table1 %>%
  dplyr::select(Population, Region, Habitat, Complex, mean_pi) %>%
  mutate(
    Complex = recode(
      Complex,
      "Complex I" = "pi_ComplexI",
      "Complex III" = "pi_ComplexIII",
      "Complex IV" = "pi_ComplexIV",
      "Complex V" = "pi_ComplexV"
    )
  ) %>%
  pivot_wider(
    names_from = Complex,
    values_from = mean_pi
  ) %>%
  arrange(Region, Habitat, Population)

print(table1_wide)

write.csv(
  table1_wide,
  file.path(out_dir, "Table_S1_mean_pi_population_complex_wide.csv"),
  row.names = FALSE
)

# ============================================================
# 6. Table S2: AK vs BC statistics by complex
#    Population-level test to avoid pseudo-replication
# ============================================================

mean_pi_pop_complex <- gene_pi_long %>%
  group_by(Population, Region, Complex) %>%
  summarise(
    mean_pi = mean(pi, na.rm = TRUE),
    .groups = "drop"
  )

table2 <- mean_pi_pop_complex %>%
  group_by(Complex) %>%
  summarise(
    AK_mean = mean(mean_pi[Region == "Alaska"], na.rm = TRUE),
    BC_mean = mean(mean_pi[Region == "BC"], na.rm = TRUE),
    p_value = wilcox.test(mean_pi ~ Region)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    fold_change = AK_mean / BC_mean,
    direction = case_when(
      AK_mean > BC_mean ~ "AK > BC",
      AK_mean < BC_mean ~ "BC > AK",
      TRUE ~ "equal"
    )
  ) %>%
  mutate(
    AK_mean = round(AK_mean, 5),
    BC_mean = round(BC_mean, 5),
    fold_change = round(fold_change, 2),
    p_value = signif(p_value, 3),
    p_adj = signif(p_adj, 3)
  )

print(table2)

write.csv(
  table2,
  file.path(out_dir, "Table_S2_population_level_AK_vs_BC_stats.csv"),
  row.names = FALSE
)

# ============================================================
# 7. Table S3: population-level pi ranking
# ============================================================

table3 <- gene_pi_long %>%
  group_by(Population, Region, Habitat) %>%
  summarise(
    mean_pi = mean(pi, na.rm = TRUE),
    sd_pi = sd(pi, na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_pi)) %>%
  mutate(
    rank = row_number(),
    top10 = ifelse(rank <= ceiling(n() * 0.1), "Top 10%", "Other")
  ) %>%
  mutate(
    mean_pi = round(mean_pi, 5),
    sd_pi = round(sd_pi, 5)
  )

print(table3)

write.csv(
  table3,
  file.path(out_dir, "Table_S3_population_pi_ranking.csv"),
  row.names = FALSE
)

# ============================================================
# 8. Statistic A:
#    Is Complex I the highest-diversity complex more often than expected?
# ============================================================

complex_winner <- mean_pi_pop_complex %>%
  group_by(Population) %>%
  slice_max(mean_pi, n = 1, with_ties = FALSE) %>%
  ungroup()

complexI_wins <- sum(complex_winner$Complex == "Complex I")
n_populations <- n_distinct(complex_winner$Population)

binom_complexI <- binom.test(
  complexI_wins,
  n_populations,
  p = 0.25,
  alternative = "greater"
)

table_stat_complexI <- tibble(
  test = "Binomial test: Complex I highest-diversity complex",
  n_populations = n_populations,
  complexI_wins = complexI_wins,
  proportion = round(complexI_wins / n_populations, 3),
  expected_probability = 0.25,
  p_value = signif(binom_complexI$p.value, 3)
)

print(table_stat_complexI)

write.csv(
  complex_winner,
  file.path(out_dir, "Table_S4_highest_complex_per_population.csv"),
  row.names = FALSE
)

write.csv(
  table_stat_complexI,
  file.path(out_dir, "Table_S5_binomial_test_complexI_dominance.csv"),
  row.names = FALSE
)

# ============================================================
# 9. Statistic B:
#    Are LG and SL higher than other populations?
# ============================================================

pop_mean <- gene_pi_long %>%
  group_by(Population, Region, Habitat) %>%
  summarise(
    mean_pi = mean(pi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    top_group = ifelse(Population %in% c("LG", "SL"), "LG_SL", "Other")
  )

wilcox_LG_SL <- wilcox.test(
  mean_pi ~ top_group,
  data = pop_mean,
  exact = TRUE
)

table_stat_population_outlier <- tibble(
  test = "Wilcoxon test: LG and SL vs other populations",
  top_group = "LG_SL",
  n_top = sum(pop_mean$top_group == "LG_SL"),
  n_other = sum(pop_mean$top_group == "Other"),
  top_mean = mean(pop_mean$mean_pi[pop_mean$top_group == "LG_SL"], na.rm = TRUE),
  other_mean = mean(pop_mean$mean_pi[pop_mean$top_group == "Other"], na.rm = TRUE),
  fold_change = top_mean / other_mean,
  p_value = wilcox_LG_SL$p.value
) %>%
  mutate(
    top_mean = round(top_mean, 5),
    other_mean = round(other_mean, 5),
    fold_change = round(fold_change, 2),
    p_value = signif(p_value, 3)
  )

print(table_stat_population_outlier)

write.csv(
  table_stat_population_outlier,
  file.path(out_dir, "Table_S6_LG_SL_vs_other_populations_wilcox.csv"),
  row.names = FALSE
)

# Optional: overall Kruskal-Wallis across populations
# Note: population-level table has one value per population, so this is not very informative.
kw_population <- kruskal.test(mean_pi ~ Population, data = pop_mean)

table_stat_kw <- tibble(
  test = "Kruskal-Wallis across populations",
  statistic = unname(kw_population$statistic),
  df = unname(kw_population$parameter),
  p_value = kw_population$p.value
) %>%
  mutate(
    statistic = round(statistic, 3),
    p_value = signif(p_value, 3)
  )

print(table_stat_kw)

write.csv(
  table_stat_kw,
  file.path(out_dir, "Table_S7_population_kruskal_wallis.csv"),
  row.names = FALSE
)

# ============================================================
# 10. Print summary text
# ============================================================

cat("\nSummary:\n")
cat("Complex I was highest in", complexI_wins, "out of", n_populations,
    "populations; binomial p =", signif(binom_complexI$p.value, 3), "\n")
cat("LG and SL vs other populations Wilcoxon p =",
    signif(wilcox_LG_SL$p.value, 3), "\n")
cat("All outputs saved to:", out_dir, "\n")