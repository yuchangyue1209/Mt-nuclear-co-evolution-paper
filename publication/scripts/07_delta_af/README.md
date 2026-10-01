# Figure 5: regional allele-frequency change

`Figure5_all_populations.R` is the current primary all-population implementation. It uses 25 freshwater populations (11 AK, 14 BC), chooses one shared representative SNP per nuclear OXPHOS gene with the same focal allele in both regions, and defines candidates by absolute regional median delta-AF ≥ 0.5 in at least one region. SNP selection prioritizes same-direction shifts and then maximizes the smaller absolute regional shift; this is not a largest mean per-SNP shift statistic. Recent populations contribute to primary regional medians.

```bash
Rscript Figure5_all_populations.R --af <ALLELE_FREQUENCY_TABLE> \
  --genes <OXPHOS72_GENE_LIST> --outdir <RESULTS_DIR>
Rscript Figure5_replot_colors_largefont.R <RESULTS_DIR> <PLOT_DIR>
Rscript Figure5C_recent_vs_established.R <RESULTS_DIR> <PLOT_DIR>
```

The saved-data replot produces A and B with large fonts, pink A candidates, blue established populations, orange recent populations, and no numeric delta-AF annotations above B. All B points are circles. Black segments are all-freshwater medians; purple segments are regional marine references. Panel C uses green filled AK and open BC circles and compares established/recent medians across the same ten candidate genes, using the unchanged shared SNPs. Its square plot and numerical summaries are written separately.

The recent populations are SC, CH, LB, PACH, and FRED. Do not add LY to this set. See `../../docs/FIGURE_REPLOTTING.md` for figure mapping and interpretation limits. The numbered established-first scripts remain for provenance and do not define current primary results.
