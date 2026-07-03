nature aging
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ERC (Evolutionary Rate Covariation) — Nature Aging (2024) 三步法复刻 + 仅限 PAML 基因
新增功能：
- 读取 PAML 汇总表 (--paml-tsv)，只对其中出现的基因跑 ERC
- 可选按 role / model / 是否有有效 ω 过滤 (--roles, --models, --require-omega)
- 其他与上一版相同：成对修剪共享 TaxonNamespace、mt 基因名归一、Pearson/Spearman、BH-FDR、LOO/perm

示例：
python erc_mt_vs_nu_py.py \
  --ref-tree /mnt/spareHD_2/oxphos_gene_tree/species_astral.tre \
  --mt-dir   /mnt/spareHD_2/mt_gene_tree/pruned_trees_te \
  --nu-dir   /mnt/spareHD_2/oxphos_gene_tree/pruned_trees_te \
  --paml-tsv /mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/codeml_sites_summary.merged.tsv \
  --roles subunit --models M0 --require-omega \
  --out-dir  /mnt/spareHD_2/oxphos_gene_tree/_erc_outputs_py_pamlOnly \
  --jackknife --perm --perm-n 1000 --debug
"""

import argparse
import os
import sys
import glob
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import dendropy

# ============ 常量（13 个 mtPCG 的标准名） ============
CANON_MT = {
    "ATP6","ATP8","COX1","COX2","COX3",
    "CYTB","ND1","ND2","ND3","ND4","ND4L","ND5","ND6"
}

def canonicalize_mt_name(name: str) -> str:
    x = name.upper()
    x = x.replace("COXI", "COX1").replace("COI", "COX1").replace("CO1", "COX1")
    x = x.replace("COII", "COX2").replace("CO2", "COX2")
    x = x.replace("COIII", "COX3").replace("CO3", "COX3")
    return x

def read_tree(path: str, tax_ns: Optional[dendropy.TaxonNamespace]=None) -> dendropy.Tree:
    return dendropy.Tree.get(path=path, schema="newick",
                             rooting="force-unrooted",
                             preserve_underscores=True,
                             taxon_namespace=tax_ns)

def load_fixed_trees_relaxed(dir_path: str, is_mt: bool=False, debug: bool=False) -> Dict[str, dendropy.Tree]:
    trees = {}
    files = glob.glob(os.path.join(dir_path, "*_fixed.treefile"))
    if debug:
        print(f"[DEBUG] scanning {dir_path}, found {len(files)} files", file=sys.stderr)
    for p in files:
        base = os.path.basename(p).replace("_fixed.treefile", "")
        gene = canonicalize_mt_name(base) if is_mt else base
        try:
            tr = read_tree(p)
            trees[gene] = tr
        except Exception as e:
            print(f"[WARN] skip {p}: {e}", file=sys.stderr)
    if not trees:
        raise RuntimeError(f"No *_fixed.treefile found in {dir_path}")
    if is_mt:
        got = set(trees.keys())
        keep = {k: trees[k] for k in got & CANON_MT}
        if debug:
            print("[DEBUG] mt detected:", sorted(got), file=sys.stderr)
            print("[DEBUG] mt kept    :", sorted(keep.keys()), file=sys.stderr)
            miss = sorted(list(CANON_MT - set(keep.keys())))
            if miss:
                print("[DEBUG] mt missing :", miss, file=sys.stderr)
        trees = keep
    return trees

# ============ 成对修剪 & 分裂位映射 ============
def co_prune_pair(ref_full: dendropy.Tree, gene_full: dendropy.Tree, taxa_to_keep: List[str]) -> Tuple[Optional[dendropy.Tree], Optional[dendropy.Tree]]:
    keep = sorted(set(taxa_to_keep))
    if len(keep) < 4:
        return None, None
    ns = dendropy.TaxonNamespace(keep)

    r = ref_full.clone(depth=2); r.retain_taxa_with_labels(keep)
    g = gene_full.clone(depth=2); g.retain_taxa_with_labels(keep)

    r_new = dendropy.Tree.get(
        data=r.as_string(schema="newick", suppress_rooting=True),
        schema="newick", rooting="force-unrooted",
        preserve_underscores=True, taxon_namespace=ns
    )
    g_new = dendropy.Tree.get(
        data=g.as_string(schema="newick", suppress_rooting=True),
        schema="newick", rooting="force-unrooted",
        preserve_underscores=True, taxon_namespace=ns
    )
    return r_new, g_new

def bipartition_length_map(tree: dendropy.Tree, min_len: float = 0.0) -> Dict[int, float]:
    m = {}
    for edge in tree.postorder_edge_iter():
        if edge.bipartition is None:
            continue
        # 排除根边（seed/root edge）
        if edge.tail_node is None:
            continue
        bl = edge.length if edge.length is not None else 0.0
        if bl >= min_len:
            m[edge.bipartition.split_bitmask] = bl
    return m

def relative_rates_vec(gene_tr: dendropy.Tree, ref_tr: dendropy.Tree, min_ref_len: float = 1e-6) -> np.ndarray:
    ref_tr.encode_bipartitions()
    gene_tr.encode_bipartitions()
    ref_map  = bipartition_length_map(ref_tr,  min_len=0.0)
    gene_map = bipartition_length_map(gene_tr, min_len=0.0)
    rel, eps = [], 1e-10
    for k, ref_len in ref_map.items():
        if ref_len <= min_ref_len:
            continue
        if k in gene_map:
            g = gene_map[k]
            if np.isfinite(g) and np.isfinite(ref_len):
                rel.append(g / max(ref_len, eps))
    return np.array(rel, dtype=float)

def mt_super_for_taxa(ref_full: dendropy.Tree, mt_trees: Dict[str, dendropy.Tree], taxa: List[str], min_ref_len: float = 1e-6) -> np.ndarray:
    vecs = []
    for mt_name in CANON_MT:
        if mt_name not in mt_trees:
            continue
        ref_sub, mt_sub = co_prune_pair(ref_full, mt_trees[mt_name], taxa)
        if ref_sub is None:
            continue
        rr = relative_rates_vec(mt_sub, ref_sub, min_ref_len=min_ref_len)
        if rr.size > 0:
            vecs.append(rr)
    if not vecs:
        return np.array([])
    L = min(v.size for v in vecs)
    if L < 3:
        return np.array([])
    M = np.vstack([v[:L] for v in vecs])
    return M.mean(axis=0)

def cor_two(x: np.ndarray, y: np.ndarray):
    pe = stats.pearsonr(x, y)
    sp = stats.spearmanr(x, y)
    return pe.statistic, pe.pvalue, sp.statistic, sp.pvalue

# ============ PAML 基因白名单相关 ============
def load_paml_genes(tsv_path: str,
                    roles: Optional[List[str]] = None,
                    models: Optional[List[str]] = None,
                    require_omega: bool = False,
                    uppercase: bool = False,
                    debug: bool = False) -> Set[str]:
    """
    从 PAML 汇总表读取满足条件的基因名集合。
    - roles: 仅保留这些 role（如 ["subunit"]）
    - models: 仅保留这些 model（如 ["M0"]）
    - require_omega: 要求 omega 非 NA 且 >0（可按需调整为 >=0）
    - uppercase: 是否将 gene 名转为大写（用于与树文件名对齐）
    """
    df = pd.read_csv(tsv_path, sep="\t")
    if roles:
        df = df[df["role"].isin(roles)]
    if models:
        df = df[df["model"].isin(models)]
    if require_omega:
        df = df[~df["omega"].isna()]
        # 如果你只想保留有有效 dN/dS 的行，也可以再加一个条件：
        # df = df[df["omega"] > 0]
    genes = df["gene"].astype(str)
    if uppercase:
        genes = genes.str.upper()
    wl = set(genes.unique().tolist())
    if debug:
        print(f"[DEBUG] PAML whitelist genes: {len(wl)}", file=sys.stderr)
    return wl

# ============ 核心：单个基因的 ERC ============
def erc_for_one_nu(nu_name: str,
                   nu_tr: dendropy.Tree,
                   ref_full: dendropy.Tree,
                   mt_trees: Dict[str, dendropy.Tree],
                   min_ref_len: float = 1e-6,
                   do_jackknife: bool = False,
                   do_perm: bool = False,
                   perm_n: int = 1000,
                   rng_seed: int = 1,
                   debug: bool = False) -> Dict:
    taxa = list(set([t.label for t in ref_full.taxon_namespace]) &
                set([t.label for t in nu_tr.taxon_namespace]))
    if len(taxa) < 4:
        if debug:
            print(f"[DBG] {nu_name}: taxa<4, skip", file=sys.stderr)
        return {}

    r_mt = mt_super_for_taxa(ref_full, mt_trees, taxa, min_ref_len=min_ref_len)
    ref_sub, nu_sub = co_prune_pair(ref_full, nu_tr, taxa)
    if r_mt.size == 0 or ref_sub is None or nu_sub is None:
        if debug:
            print(f"[DBG] {nu_name}: empty mt_super or pair prune failed", file=sys.stderr)
        return {}

    r_nu = relative_rates_vec(nu_sub, ref_sub, min_ref_len=min_ref_len)
    L = min(r_mt.size, r_nu.size)
    if L < 3:
        if debug:
            print(f"[DBG] {nu_name}: effective edges <3, skip", file=sys.stderr)
        return {}

    x, y = r_mt[:L], r_nu[:L]
    r_p, p_p, r_s, p_s = cor_two(x, y)
    out = dict(gene=nu_name, n_branches=L,
               r_pearson=r_p, p_pearson=p_p,
               r_spearman=r_s, p_spearman=p_s)

    if do_jackknife:
        jk_vals = []
        if L >= 5:
            for i in range(L):
                xi = np.delete(x, i); yi = np.delete(y, i)
                if xi.size >= 3:
                    rr, _ = stats.pearsonr(xi, yi); jk_vals.append(rr)
        out["jk_mean"] = float(np.mean(jk_vals)) if jk_vals else np.nan
        out["jk_sd"]   = float(np.std(jk_vals))  if jk_vals else np.nan
        out["jk_min"]  = float(np.min(jk_vals))  if jk_vals else np.nan
        out["jk_max"]  = float(np.max(jk_vals))  if jk_vals else np.nan
        out["jk_n"]    = int(len(jk_vals))

    if do_perm:
        rng = np.random.default_rng(rng_seed)
        null_rs = []
        for _ in range(perm_n):
            y_perm = rng.permutation(y)
            rr, _ = stats.pearsonr(x, y_perm)
            null_rs.append(rr)
        null_rs = np.array(null_rs)
        more_extreme = np.sum(np.abs(null_rs) >= abs(r_p))
        out["perm_p_two_sided"] = float((more_extreme + 1) / (perm_n + 1))
        out["perm_mean_r"] = float(np.mean(null_rs))
        out["perm_sd_r"]   = float(np.std(null_rs))

    return out

# ============ 主程序 ============
def main():
    ap = argparse.ArgumentParser(description="ERC — mt-supergene vs nu genes (Nature Aging style) — restrict to PAML genes")
    ap.add_argument("--ref-tree", required=True)
    ap.add_argument("--mt-dir",   required=True)
    ap.add_argument("--nu-dir",   required=True)
    ap.add_argument("--out-dir",  required=True)
    ap.add_argument("--paml-tsv", required=True, help="codeml_sites_summary.merged.tsv")
    ap.add_argument("--roles",    nargs="*", default=None, help="仅保留这些 role（如 subunit mt assembly-factor 等）")
    ap.add_argument("--models",   nargs="*", default=None, help="仅保留这些 model（如 M0 M1a M2a M7 M8）")
    ap.add_argument("--require-omega", action="store_true", help="仅保留 omega 非 NA 的基因")
    ap.add_argument("--uppercase-genes", action="store_true", help="将 PAML 基因名转为大写后匹配树文件名")
    ap.add_argument("--min-ref-len", type=float, default=1e-6)
    ap.add_argument("--hits-n-branches", type=int, default=15)
    ap.add_argument("--hits-fdr",        type=float, default=0.10)
    ap.add_argument("--hits-r",          type=float, default=0.30)
    ap.add_argument("--jackknife", action="store_true")
    ap.add_argument("--perm",      action="store_true")
    ap.add_argument("--perm-n",    type=int, default=1000)
    ap.add_argument("--seed",      type=int, default=1)
    ap.add_argument("--debug",     action="store_true")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    out_all  = os.path.join(args.out_dir, "erc_mt_vs_nu_results.csv")
    out_hits = os.path.join(args.out_dir, "erc_hits_nu_top.csv")
    out_skips = os.path.join(args.out_dir, "skipped_genes.txt")

    # 读参考树 & mt/nu 树
    ref_full = read_tree(args.ref_tree)
    mt_trees = load_fixed_trees_relaxed(args.mt_dir, is_mt=True, debug=args.debug)
    if len(mt_trees) < 5:
        sys.exit(f"[ERROR] 可用 mtPCG 树太少：{list(mt_trees.keys())}")
    nu_trees_all = load_fixed_trees_relaxed(args.nu_dir, is_mt=False, debug=args.debug)

    # 读 PAML 白名单
    wl = load_paml_genes(args.paml_tsv, roles=args.roles, models=args.models,
                         require_omega=args.require_omega,
                         uppercase=args.uppercase_genes, debug=args.debug)

    # 根据白名单过滤 nu 树；允许树名大小写对齐
    # 树的 key 原样是文件基名（未强制大写），我们做一个大小写映射方便匹配
    nu_by_upper = {k.upper(): (k, v) for k, v in nu_trees_all.items()}
    nu_trees = {}
    skipped = []
    for g in wl:
        gU = g.upper()
        if gU in nu_by_upper:
            orig_name, tr = nu_by_upper[gU]
            nu_trees[orig_name] = tr
        else:
            skipped.append(g)

    if args.debug:
        print(f"[DEBUG] nu trees total: {len(nu_trees_all)}; whitelist: {len(wl)}; matched: {len(nu_trees)}; skipped: {len(skipped)}", file=sys.stderr)

    if skipped:
        with open(out_skips, "w") as fh:
            for s in skipped:
                fh.write(f"{s}\n")
        print(f"[INFO] 未匹配到树文件的 PAML 基因已写出: {out_skips}  (n={len(skipped)})")

    if not nu_trees:
        sys.exit("[ERROR] 白名单与树文件无交集。请检查基因命名与大小写/下划线。")

    # 逐基因计算 ERC
    rows = []
    for gn, nu_tr in nu_trees.items():
        res = erc_for_one_nu(
            gn, nu_tr, ref_full, mt_trees,
            min_ref_len=args.min_ref_len,
            do_jackknife=args.jackknife,
            do_perm=args.perm,
            perm_n=args.perm_n,
            rng_seed=args.seed,
            debug=args.debug
        )
        if res:
            rows.append(res)
        elif args.debug:
            print(f"[DBG] {gn}: no result", file=sys.stderr)

    if not rows:
        sys.exit("[ERROR] 没有任何 ERC 结果（可能 n_branches 太小或 mt-super 为空）。")

    # 组装、FDR、排序、输出
    df = pd.DataFrame(rows).replace([np.inf, -np.inf], np.nan)
    df["padj_pearson"]  = multipletests(df["p_pearson"].fillna(1.0).values,  method="fdr_bh")[1]
    df["padj_spearman"] = multipletests(df["p_spearman"].fillna(1.0).values, method="fdr_bh")[1]
    df["abs_r"] = df["r_pearson"].abs()
    df = df.sort_values(by=["padj_pearson", "abs_r"], ascending=[True, False], na_position="last")
    df.to_csv(out_all, index=False)
    print(f"[OK] 写出全表: {out_all}  (n={len(df)})")

    hits = df[
        (df["n_branches"] >= args.hits_n_branches) &
        (df["padj_pearson"] <= args.hits_fdr) &
        (df["abs_r"] >= args.hits_r)
    ].copy()
    hits = hits.sort_values(by=["padj_pearson", "abs_r"], ascending=[True, False], na_position="last")
    hits.to_csv(out_hits, index=False)
    print(f"[OK] 写出候选: {out_hits}  (n={len(hits)})")

if __name__ == "__main__":
    main()


python erc_mt_vs_nu_py.py \
  --ref-tree /mnt/spareHD_2/oxphos_gene_tree/species_astral.tre \
  --mt-dir   /mnt/spareHD_2/mt_gene_tree/pruned_trees_te \
  --nu-dir   /mnt/spareHD_2/oxphos_gene_tree/pruned_trees_te \
  --paml-tsv /mnt/spareHD_2/oxphos_codeml_ready/09_codeml_sites_models/codeml_sites_summary.merged.tsv \
  --roles subunit --models M0 --require-omega --uppercase-genes \
  --out-dir  /mnt/spareHD_2/oxphos_gene_tree/_erc_outputs_py_subunit \
  --jackknife --perm --perm-n 1000 --debug





#consensus gene tree 
cd /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus

#!/usr/bin/env bash
set -Eeuo pipefail

CF_DIR="/mnt/spareHD_2/mt_gene_tree/counts_top2"      # 你的 mt .cf
OUT_DIR="/mnt/spareHD_2/mt_gene_tree/unconstrained"   # 输出目录
THREADS=8
MODEL="GTR+P"   # PoMo；如稀疏可改 GTR+P

mkdir -p "$OUT_DIR"
shopt -s nullglob
for CF in "$CF_DIR"/*.cf; do
  GENE=$(basename "$CF" .cf)
  echo "🌿 [mt] $GENE ..."
  # 自由搜索拓扑 + 分支长度（不使用 -te）
  iqtree2 -s "$CF" \
          -m "$MODEL" \
          -nt "$THREADS" \
          -blmin 1e-12 -blmax 100 \
          -pre "$OUT_DIR/${GENE}" \
          -quiet --safe \
  && echo "✅ [mt] $GENE done" \
  || echo "❌ [mt] $GENE failed"
done

echo "🎯 mt 树输出：$OUT_DIR/*.treefile"



OUT_PER="_mt_gene_align_nt/per_gene_trees_dna"
mkdir -p "$OUT_PER"

while read g; do
  aln="${g}_aligned_clean.fasta"            # 你明确要用这个文件名
  if [[ ! -s "$aln" ]]; then
    echo "✖ Missing $aln，跳过 $g"; continue
  fi

  echo "🌿 IQ-TREE for $g ..."
  iqtree2 -s "$aln" \
          -st DNA \
          -m GTR+F+I+G4 \
          -bb 1000 --runs 5 \
          -nt 8 --safe \
          -pre "${OUT_PER}/${g}.nt_GTRF_IG4" -quiet \
  || { echo "❌ $g failed"; continue; }

  echo "✅ ${g} done → ${OUT_PER}/${g}.nt_GTRF_IG4.treefile"
done < mt_pcg.list


# 先列出 13 个 clean 对齐文件
: > nt_align_clean.list
while read g; do
  echo "${g}_aligned_clean.fasta" >> nt_align_clean.list
done < mt_pcg.list

# 严格取交集样本拼接（不写分区文件）
python3 concat_mtpcg_nt_strict_nopart.py \
  --list nt_align_clean.list \
  --out  mt_concat_nt.strict.clean.nopart.fasta

# 用单一模型建树（无分区）
iqtree2 -s mt_concat_nt.strict.clean.nopart.fasta \
        -st DNA -m GTR+F+G4 \
        -bb 1000 --runs 5 \
        -nt 20 -mem 12G \
        -pre tmtPCG_nt_STRICT_CLEAN_NO_PART_GTRF_G4




mnt/spareHD_2/oxphos_gene_tree/species_astral.nonum.tre


pomo te

python3 erc_matrix_mt13_vs_nu72.py \
  --ref     /mnt/spareHD_2/oxphos_gene_tree/species_astral.tre \
  --mt_glob "/mnt/spareHD_2/mt_gene_tree/_fixed_only_mt13/*_fixed.treefile" \
  --nu_glob "/mnt/spareHD_2/oxphos_gene_tree/_fixed_only_nu72/*_fixed.treefile" \
  --out_r   erc_r_mt13_vs_nu72.tsv \
  --out_p   erc_p_mt13_vs_nu72.tsv
[OK] wrote erc_r_mt13_vs_nu72.tsv and erc_p_mt13_vs_nu72.tsv




