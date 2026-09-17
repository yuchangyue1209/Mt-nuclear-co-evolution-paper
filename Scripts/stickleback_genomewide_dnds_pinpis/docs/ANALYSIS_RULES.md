# Analysis rules

## dN/dS

- Invariant genes: dN=0, dS=0, ratio undefined.
- dN=0 and dS>0: ratio=0.
- dN>0 and dS=0: ratio undefined.
- No pseudocounts.
- Primary ratio estimates require at least one expected synonymous change.
- Main background contains all chromosomes; autosome-only sensitivity excludes
  chrUn and the sex chromosome chrXIX.

## piN/piS

- Sync columns 1-27 are analyzed; column 28 (Norway) is excluded.
- Negative-strand allele counts are converted to coding-strand bases.
- piN and piS are normalized by their respective site opportunities.
- Gene-population callable fraction must be at least 0.70.
- Primary genes require at least 22 of 27 eligible populations.
- Gene ratio is mean(piN)/mean(piS), not the mean of population ratios.
- Zero-denominator ratios are undefined; no pseudocounts.

## Figures

- Dotted lines show genome-wide nuclear medians.
- Wilcoxon tests are BH-adjusted within panel and metric.
- Main panels omit sample-size labels; exact n is retained in result tables.
- Panel F is visually capped at the 99th percentile, while statistical tests
  retain all observations.
