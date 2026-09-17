# Genome-wide dN/dS and piN/piS workflow

## Scope

This directory contains the reproducible coding-sequence analyses used for the
stickleback mitonuclear manuscript. The two analyses use the same canonical
stickleback gene classification but measure different evolutionary timescales.

- dN/dS: lineage-level coding divergence among 27 populations.
- piN/piS: within-population standing coding variation from Pool-seq.

## Shared gene classifications

The nuclear genes are classified as direct n-mt, indirect n-mt, or non-n-mt.
The original functional labels are retained independently: OXPHOS subunit,
assembly factor, mitochondrial/cytosolic ribosomal protein,
mitochondrial/cytosolic ARS, OXPHOS complex, and core/noncore status.

## dN/dS primary rules

- Canonical codon-ready genes: 20,347.
- Alignment-complete primary set: 17,965 genes.
- Invariant genes contribute dN=0 and dS=0; omega is undefined.
- dN=0 and dS>0 gives omega=0.
- dN>0 and dS=0 gives undefined omega; no pseudocount is used.
- The primary omega analysis requires at least one expected synonymous change.
- Main genome-wide background includes chrUn and chrXIX.
- Autosome-only sensitivity excludes chrUn and the sex chromosome chrXIX:
  16,932 genes. Genome-wide medians change by no more than 2%.

## piN/piS primary rules

- The 28-column sync map is explicit. Columns 1-27 are the manuscript
  populations; column 28 is Norway and is excluded.
- Minimum coverage is 10 reads; nuclear maximum depth is population-specific.
- Negative-strand alleles are complemented before codon classification.
- piN and piS are normalized separately by nonsynonymous and synonymous site
  opportunities. This corrects the earlier targeted-gene estimator.
- A gene-population estimate requires callable fraction >=0.70.
- A gene enters the primary analysis with at least 22 of 27 eligible populations.
- Gene piN and piS are arithmetic means across eligible populations.
- Gene piN/piS is mean(piN)/mean(piS), not the mean of population ratios.
- piN=piS=0 and piN>0,piS=0 both have undefined ratios; no pseudocount is used.
- Nuclear primary set: 18,696 genes; dN/dS-matched set: 17,918 genes.
- The 13 mt protein-coding genes are recalculated with the vertebrate
  mitochondrial code; ND6 is treated as a negative-strand gene.

## Figure 2 conventions

Each row contains piN, piS, and piN/piS (or dN, dS, and dN/dS):

- A: mtOXPHOS, nuOXPHOS, and nuclear assembly factors.
- B: OXPHOS classes partitioned by Complex I-V.
- C: cytosolic and mitochondrial ribosomal proteins.
- D: cytosolic and mitochondrial ARSs.
- E: mtOXPHOS, nuclear core, and nuclear noncore subunits.
- F: direct n-mt, indirect n-mt, and non-n-mt.

Dotted lines are genome-wide nuclear medians. Main figures omit sample-size
labels; exact n values are retained in result tables. Stars are Wilcoxon tests
with BH correction within biological panel and metric. Panel F point display is
subsampled and its visible range is capped at the 99th percentile, while all
genes remain in statistical tests.
