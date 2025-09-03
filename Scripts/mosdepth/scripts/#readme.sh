#readme
## OXPHOS CNV QC with mosdepth

This folder contains scripts to compute per-gene depth and copy-number (CN) for nuclear OXPHOS genes across pooled BAMs, to flag loci with potential paralogy/CNV that could bias dN/dS.

### Requirements
- bash, coreutils, awk
- samtools ≥ 1.10
- mosdepth ≥ 0.3.5
- Python ≥ 3.8

### Inputs
- `OXPHOS_merged_status.backfilled.tsv`: master status table with columns `present_in_fish`, `stickleback_name`
- `nuOXPHOS_v5.v3.cds_by_gene.sorted.bed` (or `.bed`): 6-column exon/CDS parts BED
- `pool_info.txt`: header + `SAMPLE SIZE ...`
- BAMs: `${SAMPLE}_rg.bam` in a single directory
- Nuclear-only reference FASTA (for mosdepth `-f`)

### Run
```bash
bash scripts/run_mosdepth_parts.sh \
  --outdir /work/cyu/oxphos_from_ref_no_biomart/06_igv \
  --backfilled /work/cyu/oxphos_from_ref_no_biomart/06_igv/OXPHOS_merged_status.backfilled.tsv \
  --parts /work/cyu/oxphos_from_ref_no_biomart/06_igv/nuOXPHOS_v5.v3.cds_by_gene.sorted.bed \
  --pool-info /work/cyu/poolseq/pool_info.txt \
  --bam-dir /mnt/spareHD_2/nuclear_with_readgroups \
  --reference /work/cyu/stickleback_nuclear_only.fa \
  --threads 8
