# Run order

## dN/dS

Run the numbered scripts under `dnds/scripts` in ascending order. Stage 02 is
intentionally absent because filtering and mask generation are consolidated in
stage 03. The completed workflow is:

1. Prepare canonical targets.
2. Call haploid population variants.
3. Build callability masks and filtered VCFs.
4. Build masked consensus CDSs.
5. Build per-gene 27-population alignments.
6. Prepare and pilot codeml M0.
7. Run genome-wide codeml.
8. Extract synonymous/nonsynonymous site information.
9. Restore invariant genes and build the master table.
10. Run primary group tests.
11. Plot the dN/dS Figure 2 panels.

## piN/piS

Run the numbered scripts under `pinpis/scripts`:

1. Stream the nuclear sync once and calculate corrected population estimates.
3. Build the nuclear gene-level table.
4. Calculate corrected mt13 estimates.
5. Combine nuclear and mitochondrial tables.
6. Run group tests and BH correction.
7. Plot the piN/piS Figure 2 panels.

The missing stage 02 reflects the final streamlined implementation, not a
missing script.
