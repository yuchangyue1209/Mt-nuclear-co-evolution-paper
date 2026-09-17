# Server script manifest

The following scripts were created and used interactively on `cyu@144.92.58.154`. Copy their final server versions into this archive before treating stages 00–07 as frozen:

```text
/work/cyu/Kuster2026_genomewide_reanalysis/00_prepare_genomewide_targets.sh
/work/cyu/Kuster2026_genomewide_reanalysis/01_call_variants_genomewide_hap1.sh
/work/cyu/Kuster2026_genomewide_reanalysis/03_make_masks_and_filtered_vcf.sh
/work/cyu/Kuster2026_genomewide_reanalysis/04_build_masked_consensus_cds.sh
/work/cyu/Kuster2026_genomewide_reanalysis/05_build_gene_alignments.py
/work/cyu/Kuster2026_genomewide_reanalysis/05_build_gene_alignments.sh
/work/cyu/Kuster2026_genomewide_reanalysis/06_prepare_codeml.py
/work/cyu/Kuster2026_genomewide_reanalysis/06_run_codeml_pilot.sh
/work/cyu/Kuster2026_genomewide_reanalysis/07_codeml_worker.sh
/work/cyu/Kuster2026_genomewide_reanalysis/07_run_codeml_genomewide.sh
```

Do not substitute early chat versions for these server copies: several were corrected for incomplete CDSs, duplicate transcript IDs, the PAML tree header (`27 1`), and restart-safe output checks.

