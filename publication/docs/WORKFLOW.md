# Workflow and verification status

| Analysis | Manuscript output | Script status |
|---|---|---|
| Raw-read trimming and QC | Analysis inputs | Curated BBDuk and FastQC workflow present |
| Nuclear raw reads to BAM | Nuclear variant inputs | Curated Bowtie2, samtools, Picard, and indexing workflow present |
| Mitochondrial raw reads to BAM/sync | Mitochondrial variant and diversity inputs | Curated Bowtie2, samtools, and Popoolation2 workflow present |
| Mitochondrial dual-caller validation | Figure 1 consensus inputs | GATK and bcftools branches, joint genotyping, overlap filtering, consensus, and phylogeny scripts present |
| Divergence and diversity | Figure 2; Figure S2; Tables S3–S4 | Final server scripts retrieved and archived; curated copies present; path cleanup pending |
| Divergence/diversity no-recent sensitivity | Figure 2 and Figure S2 sensitivity versions | 22-population filter and 22-tip alignment/tree preparation scripts present; codeml and summary reruns required on the analysis server |
| Evolutionary-rate covariation | Figure 3 | Final five-step server workflow retrieved and archived; path cleanup pending |
| Regional maximum delta-AF | Figure 5A–B | Standalone all-population analysis script present; selects representative SNPs and candidates using all 25 freshwater populations |
| Recent-population historical follow-up | Figure 5B | Implemented within the standalone Figure 5 script after the all-population analysis; recent populations are outlined and summarized separately |
| COX4I1 structural summary | COX4I1 figure | Final plotting script copied; model-generation commands and numerical distance table pending |
| Mitochondrial population analyses | Figure 1 | Server script identification pending |
| PBS and controlled mitonuclear associations | Figure 4; Figure S4 | Jesse-revised plotting and controlled-PBS scripts retrieved; path cleanup pending |
| Geographic and mitochondrial-lineage GLMs | Figure 6 | Joint latitude/longitude and mitochondrial-lineage GLM script retrieved; path cleanup pending |
| Nonsynonymous candidate discovery | Candidate table and protein figure | Reconstruction, prioritization, SIFT4G, and BLOSUM62 scripts retrieved; path cleanup pending |

The server inventory is retained in `docs/mitonuclear_server_inventory.tsv`. Selected paths are recorded in `docs/SERVER_FILES_TO_RETRIEVE.txt`; this is a reviewable retrieval manifest, not an assertion that every legacy dependency belongs in the final release.

## Figure 5 dependency rule

The primary Figure 5 analysis uses all 25 freshwater populations (11 Alaska and 14 British Columbia) to select one shared representative SNP per gene and define high-effect candidates. The five recently colonized populations (SC, CH, LB, PACH, and FRED) therefore contribute to the primary regional estimates. Their allele-frequency trajectories are subsequently labeled and summarized as a separate historical follow-up without redefining the representative SNPs or candidate set. The current standalone implementation is `scripts/07_delta_af/Figure5_all_populations.R`; the numbered scripts in the same directory document the superseded established-first workflow.

## Release criteria

A workflow is marked verified only when:

1. its exact manuscript input and output files are identified;
2. hard-coded author and server paths are replaced with configuration variables;
3. software versions and random seeds are recorded;
4. the script runs in a clean output directory;
5. its numerical output matches the manuscript figure or table.
