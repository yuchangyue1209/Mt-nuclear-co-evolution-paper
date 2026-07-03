#!/usr/bin/env bash


merge_one() {
    local fig="$1"
    local gtf="$2"
    local pi="$3"
    local out="$4"

    echo "[MERGE] $out"

    # 先检查原始 fig 和 gtf 的数据行数是否一致
    local n_fig
    local n_gtf
    n_fig=$(zcat "$fig" | tail -n +2 | wc -l)
    n_gtf=$(wc -l < "$gtf")

    echo "  fig rows: $n_fig"
    echo "  gtf rows: $n_gtf"

    if [ "$n_fig" -ne "$n_gtf" ]; then
        echo "[ERROR] row count mismatch: fig=$n_fig gtf=$n_gtf" >&2
        return 1
    fi

    {
        # header 加 pi
        zcat "$fig" | head -n 1 | awk 'BEGIN{OFS="\t"}{print $0,"pi"}'

        # data: 原始 fig 数据行 + gtf 对应的 gene_id
        paste \
          <(zcat "$fig" | tail -n +2) \
          <(awk -F'\t' 'match($9,/gene_id "([^"]+)"/,m){print m[1]}' "$gtf") \
        | awk -F'\t' 'BEGIN{OFS="\t"}
            NR==FNR{
                pi[$1]=$NF
                next
            }
            {
                key=$NF
                out=$1
                for(i=2;i<NF;i++) out=out OFS $i
                if(key in pi){
                    print out, pi[key]
                } else {
                    print out, "NA"
                }
            }' "$pi" -
    } | gzip > "$out"

    echo "[OK] wrote $out"
}

# -------------------------
# GOS
merge_one \
  /work/cyu/meth/Gosling/Gos_AllSites_Fig.tsv.gz \
  /work/cyu/meth/Gosling/Gos_AllSites_Fig.gtf \
  /work/cyu/meth/Gosling/Gos_AllSites_Fig.pi.txt \
  /work/cyu/meth/Gosling/Gos_AllSites_Fig_with_pi.tsv.gz

# ROB
merge_one \
  /work/cyu/meth/Roberts/Rob_AllSites_Fig.tsv.gz \
  /work/cyu/meth/Roberts/Rob_AllSites_Fig.gtf \
  /work/cyu/meth/Roberts/Rob_AllSites_Fig.pi.txt \
  /work/cyu/meth/Roberts/Rob_AllSites_Fig_with_pi.tsv.gz

# SAY
merge_one \
  /work/cyu/meth/Sayward/Say_AllSites_Fig.tsv.gz \
  /work/cyu/meth/Sayward/Say_AllSites_Fig.gtf \
  /work/cyu/meth/Sayward/Say_AllSites_Fig.pi.txt \
  /work/cyu/meth/Sayward/Say_AllSites_Fig_with_pi.tsv.gz

# WT
merge_one \
  /work/cyu/meth/Watson/WGBS_WT_AllSites_Fig.tsv.gz \
  /work/cyu/meth/Watson/WGBS_WT_AllSites_Fig.gtf \
  /work/cyu/meth/Watson/WGBS_WT_AllSites_Fig.pi.txt \
  /work/cyu/meth/Watson/WGBS_WT_AllSites_Fig_with_pi.tsv.gz

# WK
merge_one \
  /work/cyu/meth/Wik/WGBS_WK_AllSites_Fig.tsv.gz \
  /work/cyu/meth/Wik/WGBS_WK_AllSites_Fig.gtf \
  /work/cyu/meth/Wik/WGBS_WK_AllSites_Fig.pi.txt \
  /work/cyu/meth/Wik/WGBS_WK_AllSites_Fig_with_pi.tsv.gz