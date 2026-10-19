# Genotype-checked, masked mitochondrial consensus

This optional replacement for historical steps 04 and 05 reuses existing caller
VCFs and mt BAMs. It does not rerun mapping, HaplotypeCaller, or joint genotyping.
Historical scripts and outputs remain unchanged. Nuclear analyses are unchanged.

From the publication directory, after editing your local config/paths.sh:

```bash
bash scripts/03_variant_calling/04b_mtdna_intersect_genotype_checked.sh config/paths.sh
bash scripts/03_variant_calling/05b_mtdna_masked_consensus.sh config/paths.sh
```

Requirements: bash, bcftools, samtools, python3. The implementation was smoke-tested
with bcftools 1.20 and samtools 1.19.2; run on a pilot population and inspect logs
with your production versions before processing all populations.

## Inputs

- `$DATA_ROOT/mt_variants/bcftools/SAMPLE.vcf.gz`
- `$DATA_ROOT/mt_variants/gatk_vcf/SAMPLE.vcf.gz` (joint-genotyped, per-population VCF;
  not a GVCF)
- `$DATA_ROOT/mt_bam/SAMPLE.mt.sorted.bam`
- `MITO_REFERENCE_FASTA`, `MITO_CONTIG`, and `DATA_ROOT` in the config.

All VCF sample names must match the file sample name. Reference sequences and
coordinates must match the BAM/VCF. Use the BAM source appropriate to the analysis
being compared (for example, joint-reference validation BAMs in a separate run).

## What changes

04b selects the named population, splits/normalizes records, and retains only
records whose genotype contains that ALT. Thus GATK reference/missing genotypes
cannot count as population-specific ALT support. Exact REF/ALT intersections are
filtered using BCFtools QUAL/DP, retaining biallelic SNPs only. GATK support here
means at least one called ALT copy, not agreement on frequency or predominant
allele. Raw caller QUAL/DP are not interchangeable between programs.

05b re-pileups the full mt contig, including reference-supported positions, using
explicit MAPQ/base-quality thresholds. Every unobserved position starts masked.
A base is masked when:

- usable DP or summed AD is below MIN_DP, or DP exceeds MAX_DP;
- no canonical base has at least MIN_MAJOR_AF of total AD;
- the predominant nonreference base is absent from the passing intersection;
- a passing SNP disagrees with the predominant pileup base;
- an ALT candidate from either caller fails the intersection/SNP/quality filters;
- it overlaps a candidate indel reference span or its INDEL_PAD flanks.

The same majority-fraction rule is checked for reference and alternative bases.
All masks are applied as N in the reference coordinate system. No indel is applied
to the consensus, so reference length and downstream annotation coordinates are
preserved. This is a conservative representative sequence, not an individual
haplotype or an estimator of pool haplotype frequencies.

Defaults (configurable through exported environment variables):

| Variable | Default |
|---|---:|
| MIN_DP | 10 |
| MAX_DP | 5000 |
| MIN_QUAL (04b) | 30 |
| MIN_MAPQ (05b) | 20 |
| MIN_BASEQ (05b) | 13 |
| MIN_MAJOR_AF (05b) | 0.80 |
| INDEL_PAD (05b) | 5 bp |
| PILEUP_CAP (05b) | 100000 |

The 0.80 majority threshold and 5 bp indel flanks are explicit new conservative
choices, not claims about the historical analysis. Assess their sensitivity before
updating manuscript results. The high pileup cap allows the MAX_DP mask to detect
excessive coverage rather than truncating all sites to MAX_DP. Depths at the cap
are necessarily masked because the cap must be greater than MAX_DP.

Example overrides (use the same exported settings for both steps):

```bash
export MIN_DP=10 MAX_DP=5000 MIN_QUAL=30
export MIN_MAPQ=20 MIN_BASEQ=13 MIN_MAJOR_AF=0.80 INDEL_PAD=5
```

## Outputs and comparison

- `$DATA_ROOT/mt_variants/genotype_checked/SAMPLE/`: normalized ALT-supported
  caller files, shared records, and passing SNP VCF.
- `$DATA_ROOT/mt_consensus_masked/SAMPLE.fasta`
- `$DATA_ROOT/mt_consensus_masked/all_populations.fasta`
- `$DATA_ROOT/mt_consensus_masked/SAMPLE/mask.bed`, `mask_qa.json`, intermediate
  pileup/candidate tables, and `parameters.txt`.

Reason counts can overlap; their sum is not the total masked length. Inspect
callable_bp/masked_bp and the BED distribution per population. Empty pass sets
and entirely masked sequences are reported rather than assumed successful
biological reconstructions.

Existing step 06 still reads the historical mt_consensus directory. Do not run it
expecting masked input. Use mt_consensus_masked/all_populations.fasta explicitly
in a separate alignment/tree run, apply the same control-region exclusion as in
the original analysis, and write to a separate results directory. Compare base
changes on jointly callable sites, new missingness, diagnostic SNPs, and C1-C4
assignments. Update mt CDS/dN/dS inputs if they use these sequences. Add the
masking procedure to Methods only after running and checking the actual results.

## Local validation

Synthetic checks covered reference/ALT majority support, mixed sites, low/high
and zero coverage, rejected variants, indel spans, empty pass sets, and BED
coordinates. An end-to-end synthetic BAM/VCF run confirmed that a true shared
ALT is applied, a GATK 0/0 record is excluded and masked, absent coverage becomes
N, and sequence length is preserved. Real population data were not run locally.
