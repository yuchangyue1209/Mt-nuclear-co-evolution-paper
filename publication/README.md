# Repeated freshwater adaptation and mitonuclear evolution in stickleback

This repository contains the analysis code used to examine mitochondrial variation, molecular evolution, evolutionary-rate covariation, population differentiation, parallel nuclear OXPHOS allele-frequency change, recurrent nonsynonymous substitutions, and structural predictions in threespine stickleback.

## Repository status

The publication repository is being reconstructed from the verified scripts used for the manuscript. Scripts are copied here only after their role and inputs have been confirmed. Legacy and exploratory scripts remain in the private project archive and are not part of this release.

## Workflow

1. Read processing and variant calling
2. Mitochondrial consensus sequences, phylogeny, and diversity
3. Functional gene-set and orthology curation
4. Gene-level divergence and polymorphism
5. Evolutionary-rate covariation
6. Population branch statistic analyses
7. Regional allele-frequency change and recent-population evaluation
8. Recurrent nonsynonymous substitutions
9. COX4I1 structural analysis
10. Manuscript figures and tables

See `docs/WORKFLOW.md` for the current script map, `docs/PROVENANCE.md` for the relationship between archived and curated code, and `config/paths.example.sh` for external data configuration.

## Data availability

Large sequencing files, BAM files, VCF files, genome assemblies, and intermediate alignments are not stored in GitHub. Public accession numbers and download instructions will be listed in `docs/DATA_AVAILABILITY.md`. Small metadata and final numerical tables required to reproduce figures will be included when redistribution is permitted.

## Reproducibility

Personal workstation and server paths have been replaced with generic placeholders such as `/path/to/repository`, `/path/to/workspace`, and `/path/to/data`. Users should adapt these aliases through the configuration template before running the workflows. Original server copies and the full server inventory are retained separately in the private project backup and are not part of this public release.

The release includes separate workflows from paired raw reads through indexed nuclear and mitochondrial BAM files. Nuclear BAMs follow the historical Pool-seq workflow with MAPQ filtering, duplicate removal, and read-group assignment. Mitochondrial reads are mapped independently and are not deduplicated; indexed mitochondrial BAMs can subsequently be converted to mpileup and sync allele-count formats.
