# 先跑覆盖
bash run_mosdepth_combined.sh

# 汇总 depth / CN
python build_gene_depth_cn_combined.py \
  /work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear \
  /work/cyu/poolseq/pool_info.txt \
  /work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/parts.present+putative.4col.bed \
  /work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.COMB.gene_depth.tsv \
  /work/cyu/oxphos_from_ref_no_biomart/06_igv/cnv_qc_backfilled_nuclear/oxphos_cnv.COMB.gene_cn.tsv
