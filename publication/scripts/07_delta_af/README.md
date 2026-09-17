# Figure 5: regional allele-frequency change

## Current manuscript workflow

`Figure5_all_populations.R` is the standalone implementation used for the current Figure 5 analysis. It:

- uses all 25 freshwater populations (11 Alaska and 14 British Columbia);
- selects one shared representative SNP for each of the 72 nuclear OXPHOS genes;
- uses the same focal allele in both regions;
- defines high-effect genes using the all-population regional median delta-AF values;
- produces panels A and B, the combined figure, numerical plot data, and audit tables; and
- labels the five recently colonized populations and summarizes their trajectories separately after the primary all-population analysis.

Run it as:

```bash
Rscript --vanilla Figure5_all_populations.R \
  --af <ALLELE_FREQUENCY_TABLE> \
  --genes <OXPHOS72_GENE_LIST> \
  --outdir <FIGURE5_OUTPUT>
```

Optional arguments are `--min-depth`, `--threshold`, and `--dpi`.

## Legacy workflow

The numbered scripts `01_select_shared_snp_and_plot_established.R`, `02_evaluate_recent_populations.R`, and `03_plot_established_and_recent.R` preserve the earlier established-first analysis. They are retained for provenance but do not define the current primary Figure 5 results.
