# Figure 2 and Figure S2 sensitivity analysis without recent populations

The primary analysis retains all 27 populations. This sensitivity analysis
excludes SC, CH, LB, PACH, and FRED and writes to a separate `no_recent`
directory. It must not overwrite the primary tables or figures.

For piN/piS, run `01_filter_pinpis_populations.R`, then rebuild gene-level
summaries, group tests, and Figure S2 from the filtered population tables.
Require at least 18 of 22 eligible populations, preserving approximately the
same eligibility fraction as the primary criterion of 22 of 27.

For dN/dS, run `02_prune_codeml_alignments.py` on both nuclear and mitochondrial
codon alignments and run `03_prune_tree.R` on the fixed topology. Re-run codeml
with these 22-tip inputs before rebuilding the master table, group tests, and
Figure 2. Filtering the existing 27-tip dN/dS table is not valid.

Compare group medians, effect directions, BH-adjusted significance, and gene
counts between the primary and no-recent runs.
