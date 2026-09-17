# Server script manifest

The following scripts were created and used interactively on `user@server`. Copy their final server versions into this archive before treating stages 00–07 as frozen:

```text
/path/to/workspace/Kuster2026_genomewide_reanalysis/00_prepare_genomewide_targets.sh
/path/to/workspace/Kuster2026_genomewide_reanalysis/01_call_variants_genomewide_hap1.sh
/path/to/workspace/Kuster2026_genomewide_reanalysis/03_make_masks_and_filtered_vcf.sh
/path/to/workspace/Kuster2026_genomewide_reanalysis/04_build_masked_consensus_cds.sh
/path/to/workspace/Kuster2026_genomewide_reanalysis/05_build_gene_alignments.py
/path/to/workspace/Kuster2026_genomewide_reanalysis/05_build_gene_alignments.sh
/path/to/workspace/Kuster2026_genomewide_reanalysis/06_prepare_codeml.py
/path/to/workspace/Kuster2026_genomewide_reanalysis/06_run_codeml_pilot.sh
/path/to/workspace/Kuster2026_genomewide_reanalysis/07_codeml_worker.sh
/path/to/workspace/Kuster2026_genomewide_reanalysis/07_run_codeml_genomewide.sh
```

Do not substitute early chat versions for these server copies: several were corrected for incomplete CDSs, duplicate transcript IDs, the PAML tree header (`27 1`), and restart-safe output checks.

