# ERC analysis

Reproducible genome-wide mitonuclear evolutionary-rate-correlation workflow for the stickleback manuscript.

## Run order on `stickleback`

1. `run_mt13_codeml.py` — estimate branch lengths for all 13 mitochondrial PCGs using the fixed 27-population topology.
2. `extract_mt13_dnds.py` — collect mt codeml results.
3. `run_genomewide_erc.py` — calculate gene-specific nuclear and mt relative evolutionary rates and genome-wide ERC.
4. `R/00_gene_level_screening.R` — empirical non-n-mt background percentiles and category-level screening.
5. `R/01_weaver_gene_set_erc.R` — Weaver-style gene-set composites and 10,000 gene bootstrap replicates.
6. `R/02_complex_permutation_and_CIV_robustness.R` — 100,000 branch permutations for CI–CV and CIV leave-one-out tests.
7. `R/03_plot_ERC_main_figure.R` — final two-panel ERC figure.

The fixed phylogenetic topology is shared by all genes, whereas branch lengths are estimated separately for each alignment. The main figure uses the Weaver-style bootstrap criterion: red denotes a positive ERC whose 95% gene-bootstrap confidence interval excludes zero. Complex II is completely nuclear encoded and is compared with all 13 mtPCGs as a negative control (`†`).

## Server execution

```bash
bash /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/run_erc_server.sh

Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/R/00_gene_level_screening.R
Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/R/01_weaver_gene_set_erc.R
Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/R/02_complex_permutation_and_CIV_robustness.R
Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/R/03_plot_ERC_main_figure.R
```

Primary results are written to:

`/path/to/data/genomewide_codeml_kuster/10_erc_results/mtPCG_composite_t_spearman`

The final figure is written to:

`/path/to/data/genomewide_codeml_kuster/11_erc_figures/Figure_ERC_Weaver_style`
