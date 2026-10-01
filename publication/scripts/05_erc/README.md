# ERC analysis

Reproducible genome-wide mitonuclear evolutionary-rate-correlation workflow for the stickleback manuscript.

## Run order on `stickleback`

1. `run_mt13_codeml.py` — estimate branch lengths for all 13 mitochondrial PCGs using the fixed 27-population topology.
2. `extract_mt13_dnds.py` — collect mt codeml results.
3. `run_genomewide_erc.py` — calculate gene-specific nuclear and mt relative evolutionary rates and genome-wide ERC.
4. `final_R/01_erc_category_statistics.R` — gene-level ERC screening and category summaries.
5. `final_R/02_weaver_gene_set_erc.R` — Weaver-style gene-set composites and 10,000 gene-bootstrap replicates. This main analysis includes the non-n-mt and genome-wide nuclear sets as descriptive empirical backgrounds in Figure 3A.
6. `final_R/03_validate_CIV_ERC.R` — leave-one-out robustness checks for Complex IV.
7. `final_R/04_all_complex_permutation_and_membership_audit.R` — branch-permutation tests and complex-membership audit.
8. `final_R/05_plot_Weaver_style_ERC.R` — final two-panel ERC figure.

The fixed phylogenetic topology is shared by all genes, whereas branch lengths are estimated separately for each alignment. The main figure uses the Weaver-style bootstrap criterion: red denotes a positive ERC whose 95% gene-bootstrap confidence interval excludes zero. Panel A includes all non-n-mt genes and all quality-eligible nuclear genes as complementary empirical backgrounds; these sets overlap substantially and should not be interpreted as independent controls. Complex II is completely nuclear encoded and is compared with all 13 mtPCGs as a negative control (`†`).

`final_R/06_complex_vs_genome_background.R` is retained as an exploratory sensitivity analysis. It tests whether each complex ERC exceeds size-matched random nuclear gene sets from the non-n-mt and genome-wide backgrounds. It is not part of the primary Figure 3 workflow and should not be used to claim that Complex IV exceeds the genomic background unless that comparison is supported by the empirical permutation test.

## Server execution

```bash
bash /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/run_erc_server.sh

Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/final_R/01_erc_category_statistics.R
Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/final_R/02_weaver_gene_set_erc.R
Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/final_R/03_validate_CIV_ERC.R
Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/final_R/04_all_complex_permutation_and_membership_audit.R
Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/final_R/05_plot_Weaver_style_ERC.R
```

Optional exploratory comparison:

```bash
N_PERM=100000 Rscript /path/to/workspace/Kuster2026_genomewide_reanalysis/genomewide_erc/final_R/06_complex_vs_genome_background.R
```

Primary results are written to:

`/path/to/data/genomewide_codeml_kuster/10_erc_results/mtPCG_composite_t_spearman`

The final figure is written to:

`/path/to/data/genomewide_codeml_kuster/11_erc_figures/Figure_ERC_Weaver_style`
