## Current Figure 6 update

Use `scripts/08_geo_mtlineage/03_run_Figure6_associations.R` with LD-pruned data, metadata, gene list, output directory, and permutation count. No nuclear-PC input. A uses geography versus intercept; B uses mitochondrial cluster versus intercept. Both use 19 long-term freshwater populations excluding AMO. See scripts/08_geo_mtlineage/README.md for current commands. Any treePC-controlled commands below are historical alternatives.

# Current manuscript figure code

These plotting updates change presentation and panel organization. Saved statistical results remain the source of significance annotations. Replotting does not establish a new independent test of selection or causality.

| Figure | Contents | Current entry point |
|---|---|---|
| 2 | Divergence: OXPHOS roles, RP, ARS (A–C) | `scripts/04_divergence_diversity/R/12_export_manuscript_Figure2_S2_S5.R` |
| S2 | Diversity: OXPHOS roles, RP, ARS (A–C) | Same exporter |
| S3 | OXPHOS complexes; divergence A, diversity B | Same exporter |
| S4 | Core/noncore; divergence A, diversity B | Same exporter |
| S5 | Direct/indirect/non-n-mt; divergence A, diversity B | Same exporter |
| 5 A–B | Shared representative SNP, all-population regional shifts and population allele frequencies | `scripts/07_delta_af/Figure5_all_populations.R` (analysis), `Figure5_replot_colors_largefont.R` (saved data) |
| 5 C | Recent versus established regional median shifts across the ten candidates | `scripts/07_delta_af/Figure5C_recent_vs_established.R` |
| 6 A–C | Geography and treePC-controlled mitochondrial-cluster GLMs; within-gene signed profiles | `scripts/08_geo_mtlineage/03_run_treePC_controlled_Figure6.R` |
| S9 (legacy) | Archived regional distance contrasts; not yet regenerated under the current Figure 6 workflow | `scripts/08_geo_mtlineage/Figure6_and_S9_largefont.R` |
| 7 | COX4I1 amino-acid frequencies, sequence windows, protein map | `scripts/09_structural_analysis/plot_COX4I1_figure.R` |

S6 is the PBS schematic, S7 the non-subunit PBS comparisons, and S8 the LD heatmap. Their source plotting workflows are unchanged by this update. Figures 1, 3, and 4 also retain their existing scripts; a blanket font change to these unreviewed figures has not been applied.

## Replotting commands

Run from the publication directory. Substitute actual input and output locations. R packages: data.table, ggplot2, ggrepel, patchwork, ggsignif; watershed analysis additionally uses readxl.

```bash
Rscript scripts/04_divergence_diversity/R/12_export_manuscript_Figure2_S2_S5.R \
  "$NUCLEAR_CODEML_TSV" "$MT_CODEML_TSV" "$PINPIS_ROOT" "$FIGURE2_OUT"
Rscript scripts/07_delta_af/Figure5_replot_colors_largefont.R "$FIGURE5_RESULTS" "$FIGURE5_OUT"
Rscript scripts/07_delta_af/Figure5C_recent_vs_established.R "$FIGURE5_RESULTS" "$FIGURE5_OUT"
Rscript scripts/08_geo_mtlineage/03_run_treePC_controlled_Figure6.R \
  "$LD_PRUNED_TSV" scripts/08_geo_mtlineage/inputs/freshwater19_metadata.tsv \
  "$NUCLEAR_PC_TSV" scripts/08_geo_mtlineage/inputs/OXPHOS72_genes.list \
  "$FIGURE6_OUT" 9999
COX4I1_FREQ_FILE="$COX4I1_FREQ_FILE" COX4I1_PROTEIN_FILE="$COX4I1_PROTEIN_FILE" \
  COX4I1_FIGURE_DIR="$FIGURE7_OUT" Rscript scripts/09_structural_analysis/plot_COX4I1_figure.R
```

The Figure 2 exporter evaluates the archived row-building/statistical code before its legacy export stage. Explicit arguments override its input paths. It combines the raw rows with new manuscript labels and a fixed 0.8-mm frame. It does not execute the legacy A–F image exports. The pi root must contain `results/genomewide_nuclear_mt13_pinpis_combined.tsv` and `results/pinpis_group_tests.tsv`. The archived pi panel F retains its display cap and dot subsampling; significance tests use all eligible genes.

Figure 5 A candidates are pink; B established freshwater populations are blue and recent populations orange, using circles throughout. B black segments summarize **all** freshwater populations in each region, including recent populations; purple segments show marine references. Signed delta-AF numbers have been removed. C is square, with green filled AK and open BC circles; its medians use established and recent populations separately. It also writes `Figure5C_correlation_statistics.tsv` with descriptive Pearson/Spearman correlations and directional/magnitude counts. Cross-sectional medians are not longitudinal trajectories. SNP selection prioritizes concordant signs, so the same-direction fraction is not an independent binomial test of parallelism.

The current Figure 6 is 18 × 18 inches, with separate A/B panels 9 × 9 and C 13 × 9. It uses 19 long-term freshwater populations excluding AMO. B controls for nuclear treePC1/treePC2; C retains multiple SNPs within each gene. Separate one-SNP-per-gene cross-gene profiles are exported as additional statistics. Earlier Figure 6/S9 plotting settings are legacy and must not be used to imply the new model was fitted. Figure 7 uses a 9 × 18-inch layout. The duplicate Figure 7 entry point in `scripts/10_structure/01_plot_cox4i1_figure.R` is kept identical.

## Watershed analysis

```bash
Rscript scripts/03_variant_calling/test_watershed_mtcluster.R \
  "$TABLE_S1_XLSX" "$MT_CLUSTER_TSV" "$WATERSHED_OUT" 9999
```

This accepts blank lines in the cluster TSV. AMO is excluded from the four-cluster comparison; marine populations are excluded and recent freshwater populations retained. Watersheds are nested within region. Permutations shuffle cluster assignments within regions. Leave-one-population-out accuracy uses fractional credit for ties and compares watershed and region predictions on the same evaluable populations; singleton watersheds are not scored. The manuscript analysis used 24 freshwater populations and 19 cross-validation-evaluable populations. Keep the output audits with the analysis records; an additional manuscript table or figure is not necessary.

## Verification scope

The new and edited scripts have been syntax-checked locally. Figure 5C numerical summaries were checked against saved data. Full rendering of the new Figure 2 exporter and tall Figure 7 requires the server dependencies and inputs. Figure 6/S9 final enlarged layout should be reviewed at manuscript display size after server rendering. No raw-data recomputation is implied by successful syntax checks.

The current treePC-controlled Figure 6 entry point and its input/filter definitions are documented in `../scripts/08_geo_mtlineage/README.md`. Historical S8/S9 numbers and figure outputs are not automatically recomputed.
