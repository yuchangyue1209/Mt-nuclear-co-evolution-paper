# RUN_ALL.sh  —— 一键顺序跑完
#!/usr/bin/env bash
set -euo pipefail
bash 00_build_present_parts.sh
bash 01_add_putative_and_diamond_spans.sh
bash 02_run_mosdepth.sh
python 03_collect_regions_to_tables.py --regions-glob "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/*.COMB.regions.bed.gz" \
  --out-depth "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.comb_gene_depth.tsv"
python 04_depth_to_cn.py \
  --depth "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.comb_gene_depth.tsv" \
  --pool-info "/work/cyu/poolseq/pool_info.txt" \
  --out-cn "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.comb_gene_cn.tsv"
python 05_join_cnv_into_status.py \
  "/work/cyu/oxphos_from_ref_no_biomart/06_igv/OXPHOS_merged_status.backfilled.tsv" \
  "/work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.comb_gene_cn.tsv" \
  > "/work/cyu/oxphos_from_ref_no_biomart/06_igv/OXPHOS_merged_status.backfilled.withCN.plus.tsv"
bash 06_make_safe_subunits.sh
