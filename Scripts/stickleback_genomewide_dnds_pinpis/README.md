# Stickleback genome-wide dN/dS and piN/piS

Reproducible coding-sequence analyses for the stickleback mitonuclear
coevolution manuscript. The archive is organized by analysis rather than by
the source publication used to define gene classes.

## Modules

```text
dnds/scripts/      lineage-level coding divergence, stages 00-11
pinpis/scripts/    within-population standing coding variation, stages 01-07
shared/config/     server paths and expected dimensions
docs/              run order and analysis rules
```

## Final analysis sets

- Classified nuclear genes: 20,426.
- Codon-ready genes: 20,347.
- dN/dS primary genes: 17,965.
- piN/piS primary nuclear genes: 18,696.
- Matched dN/dS-piN/piS genes: 17,918.
- Corrected mitochondrial piN/piS: 13 protein-coding genes.

The expensive data products remain on the server. This directory archives
code and documentation only.

## Plotting environment

Use `/opt/conda/envs/r_env/bin/Rscript` for figures. Do not use the
`poolseq_env` R installation for graphics because its graphics API is broken.
