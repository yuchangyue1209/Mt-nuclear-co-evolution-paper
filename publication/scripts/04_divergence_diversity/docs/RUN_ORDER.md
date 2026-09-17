# Verified run order

## Upstream products already completed on the server

1. Human/mouse Kuster accessions mapped to stickleback proteins.
2. Final classification: 20,426 genes (151 direct, 862 indirect, 19,413 non).
3. Canonical CDS curation: 20,347 genes; 79 reference models excluded.
4. Haploid Pool-seq variant calls: 27/27 populations.
5. Per-population callability masks: `10 <= DP <= 2 x median DP`.
6. Filtered variants: biallelic SNP, QUAL >= 30, DP >= 10, ALT fraction >= 0.80.
7. Masked consensus and canonical CDS: 27 files x 20,347 sequences.
8. Gene alignments: 20,347 genes x 27 populations.
9. Primary completeness set: 17,965 genes.
10. codeml M0: 15,241 variable genes, zero failures; 2,724 invariant genes.

## Final lightweight stages

```bash
source config/paths.sh

python3 python/08_extract_codeml_sites.py
Rscript R/09_build_codeml_master_table.R
Rscript R/10_primary_codeml_tests.R
Rscript R/11_plot_Figure2_complete.R
```

Expected master table size: 17,966 lines including the header.

## Figure 2 structure

- A: mtOXPHOS / nuOXPHOS / assembly factors
- B: the same classes within complexes I–V
- C: cyto-RP / Nmt-RP
- D: cyto-ARS / Nmt-ARS
- E: mtOXPHOS / nu-core / nu-noncore
- F: direct n-mt / indirect n-mt / non-n-mt

Each row is saved separately as PDF and PNG, followed by a combined 6 x 3 figure. Every panel contains a dotted genome-wide median.

## Reproducibility caution

Rows A, B, and E currently combine the legacy mtOXPHOS M0 table with newly estimated nuclear rates. Before publication, verify that the mt and nuclear analyses use the same 27 populations, topology, genetic code, and M0 settings. Otherwise rerun the mt genes with the finalized topology.

