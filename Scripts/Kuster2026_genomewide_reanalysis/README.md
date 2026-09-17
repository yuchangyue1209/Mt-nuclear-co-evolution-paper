# Genome-wide mitonuclear dN/dS and piN/piS reanalysis

Clean code archive for the Kuster-style stickleback genome-wide coding-divergence,
standing-variation, and updated manuscript Figure 2 analyses.

## Scope

- 20,426 classified nuclear genes: 151 direct n-mt, 862 indirect n-mt, 19,413 non-n-mt.
- 20,347 codon-ready genes after excluding 79 incomplete/duplicated reference CDS models.
- 17,965 primary genes after alignment completeness filtering.
- 15,241 variable genes successfully analyzed by codeml; 2,724 invariant genes are restored as `dN=0`, `dS=0`, `omega=NA`.
- Updated Figure 2 rows A–F, including core/noncore and Kuster categories.
- 18,696 piN/piS-eligible nuclear genes; 17,918 overlap the dN/dS primary set.
- Corrected piN/piS estimator includes strand correction and separate N/S opportunities.
- Corrected mitochondrial piN/piS estimates for all 13 protein-coding genes.

## Directory layout

```text
config/paths.sh                 server paths and constants
docs/RUN_ORDER.md               ordered workflow and verified counts
docs/ANALYSIS_RULES.md          zero-rate and omega rules
docs/Figure2_caption.md         current caption
python/08_extract_codeml_sites.py
python/pinpis/01_genomewide_pinpis.py
python/pinpis/04_mt13_pinpis_corrected.py
R/09_build_codeml_master_table.R
R/10_primary_codeml_tests.R
R/11_plot_Figure2_complete.R
R/pinpis/03_build_gene_level_pinpis.R       synced from server
R/pinpis/05_build_combined_pinpis_table.R   synced from server
R/pinpis/06_test_pinpis_groups.R            synced from server
R/pinpis/07_plot_pinpis_Figure2_large_fonts.R
shell/run_09_to_11.sh           convenience driver for finished inputs
shell/pinpis/sync_verified_pinpis_server_scripts.sh
docs/DNDS_PINPIS_WORKFLOW.md
```

The computationally expensive server products remain under
`/mnt/spareHD_2/genomewide_codeml_kuster`. This repository contains code and small metadata only, not VCFs, masked genomes, alignments, or per-gene `mlc` files.

## Start here

1. Open this folder in VS Code.
2. Review `config/paths.sh`.
3. Read `docs/RUN_ORDER.md` for dN/dS and `docs/DNDS_PINPIS_WORKFLOW.md`
   for the shared dN/dS-piN/piS analysis rules.
4. On the server, run `shell/run_09_to_11.sh` only after all listed inputs pass the documented counts.

## Important

The `poolseq_env` R 4.3.3 installation reports a base Graphics API mismatch.
Statistical scripts are unaffected. Render figures with `/opt/conda/envs/r_env/bin/Rscript`;
do not repair or upgrade `poolseq_env` in place.
