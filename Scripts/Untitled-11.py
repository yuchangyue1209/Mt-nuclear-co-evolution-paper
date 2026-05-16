#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, sys, glob, csv, math, copy, re
from collections import defaultdict
from itertools import combinations
from typing import Dict, List, Tuple, Set

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
import numpy as np
from scipy import stats
from numpy.linalg import LinAlgError

EPS = 1e-8  # 极小正值，修复 0/缺失分支长度 & 协方差抖动

# ----------------------- 名字规范化 & 映射 -----------------------

def load_name_map(path: str) -> Dict[str, str]:
    mp = {}
    if not path:
        return mp
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t ]+", line)
            if len(parts) < 2:
                continue
            old, new = parts[0], parts[1]
            mp[old] = new
    return mp

def build_norm_fn(strip_square: bool, drop_prefix_number: bool, upper: bool, name_map: Dict[str,str]):
    def norm_name(s: str) -> str:
        if s is None:
            return None
        raw = s
        s = s.strip().strip("'").strip('"')
        if strip_square and len(s) >= 2 and s[0] == "[" and s[-1] == "]":
            s = s[1:-1]
        if drop_prefix_number:
            s = re.sub(r"^\d+[_-]?", "", s)
        if upper:
            s = s.upper()
        # 应用映射：原始键优先，其次规范化后键
        if raw in name_map:
            s = name_map[raw]
        elif s in name_map:
            s = name_map[s]
        return s
    return norm_name

# ----------------------- Tree 工具 -----------------------

def ensure_positive_branch_lengths(tr: Tree, eps: float = EPS) -> None:
    """将 None 或 <=0 的 branch_length 置为 eps。"""
    for cl in tr.find_clades(order="postorder"):
        if cl.branch_length is None or cl.branch_length <= 0:
            cl.branch_length = eps

def binarize_polytomies(tr: Tree, eps: float = EPS) -> None:
    """
    将多分叉转为二叉：对 >2 子节点的内部节点，迭代合并后两个子树为一个新父节点，
    给新内部边长设置为 eps。
    """
    def _binarize(clade: Clade):
        # 深度优先
        for ch in list(clade.clades):
            _binarize(ch)
        while len(clade.clades) > 2:
            a = clade.clades.pop()
            b = clade.clades.pop()
            new_parent = Clade(branch_length=eps)
            new_parent.clades = [a, b]
            clade.clades.append(new_parent)
    _binarize(tr.root)

def reroot_by_outgroup_or_midpoint(tr: Tree, outgroup: str) -> None:
    """优先用 outgroup 重根；否则 midpoint。异常时跳过。"""
    tips = {t.name: t for t in tr.get_terminals() if t.name}
    try:
        if outgroup and outgroup in tips:
            tr.root_with_outgroup(tips[outgroup])
        else:
            tr.root_at_midpoint()
    except Exception:
        # fallback：随便选一个 tip 重根（保证有根）
        try:
            any_tip = next(iter(tips.values()))
            tr.root_with_outgroup(any_tip)
        except Exception:
            pass

def prune_to_species(tr: Tree, keep: Set[str]) -> Tree:
    """复制并裁剪到 keep 物种集合（只保留这些名字的 tips）。"""
    tr2 = copy.deepcopy(tr)
    # 先删除不在 keep 的 tips
    for tip in list(tr2.get_terminals()):
        if tip.name not in keep:
            tr2.prune(tip, preserve_branch_length=True)
    return tr2

def normalize_tree_tip_names(tr: Tree, norm_name_fn) -> None:
    for tip in tr.get_terminals():
        if tip.name:
            tip.name = norm_name_fn(tip.name)

# ----------------------- 路径读取 -----------------------

def iter_tree_files(tree_dir: str) -> List[str]:
    pats = ["**/*.tre","**/*.tree","**/*.nwk","**/*.newick","**/*.txt",
            "**/*.treefile","**/*.contree","**/*"]
    skip_ext = {".fasta",".fa",".fas",".aln",".phy",".csv",".tsv",
                ".xlsx",".xls",".pdf",".png",".svg",".ckp",".gz",".log",".ckp.gz"}
    out, seen = [], set()
    for p in pats:
        for path in glob.glob(os.path.join(tree_dir, p), recursive=True):
            if os.path.isdir(path): 
                continue
            ext = os.path.splitext(path)[1].lower()
            if ext in skip_ext:
                continue
            if path not in seen:
                seen.add(path); out.append(path)
    return sorted(out)

def read_newick_any(path: str) -> List[Tree]:
    trees = []
    try:
        trees = list(Phylo.parse(path, "newick"))
    except Exception:
        trees = []
    if not trees:
        try:
            trees = [Phylo.read(path, "newick")]
        except Exception:
            return []
    return trees

# ----------------------- RTT 向量（根到叶距离） -----------------------

def root_by_outgroup_or_midpoint_and_norm(tr: Tree, outgroup: str, norm_name_fn) -> Tree:
    tr = copy.deepcopy(tr)
    normalize_tree_tip_names(tr, norm_name_fn)
    reroot_by_outgroup_or_midpoint(tr, outgroup)
    ensure_positive_branch_lengths(tr, EPS)
    return tr

def tree_root_to_tip_vector(tr: Tree) -> Dict[str, float]:
    depths = tr.depths()
    out = {}
    for tip in tr.get_terminals():
        if tip in depths and tip.name:
            out[tip.name] = float(depths[tip])
    return out

def transform_lengths(v: Dict[str, float], kind: str) -> Dict[str, float]:
    if kind == "none":
        return v
    out = {}
    if kind == "sqrt":
        for k, x in v.items():
            out[k] = math.sqrt(x) if x >= 0 else float("nan")
    elif kind == "log1p":
        for k, x in v.items():
            out[k] = math.log1p(x) if x >= 0 else float("nan")
    else:
        raise ValueError(f"Unknown transform: {kind}")
    return out

def load_gene_vectors(tree_dir: str, ref_species: Set[str], outgroup: str,
                      transform: str, norm_name_fn) -> Dict[str, Dict[str,float]]:
    out = {}
    if not (tree_dir and os.path.isdir(tree_dir)):
        print(f"[WARN] Dir not found or empty: {tree_dir}")
        return out
    files = iter_tree_files(tree_dir)
    print(f"[INFO] Scanning {tree_dir} … found {len(files)} candidate files")
    n_files_parsed = n_trees_ok = n_with_overlap = 0

    for path in files:
        trees = read_newick_any(path)
        if not trees:
            continue
        n_files_parsed += 1
        for idx, tr in enumerate(trees):
            try:
                tr = root_by_outgroup_or_midpoint_and_norm(tr, outgroup, norm_name_fn)
                # 不强制二叉化 gene 树；RTT 不需要
                vec = tree_root_to_tip_vector(tr)
                n_trees_ok += 1
            except Exception:
                continue
            vec = {sp: val for sp, val in vec.items() if sp in ref_species and np.isfinite(val)}
            if len(vec) == 0:
                continue
            n_with_overlap += 1
            vec = transform_lengths(vec, transform)
            base = os.path.splitext(os.path.basename(path))[0]
            gene = base if len(trees) == 1 else f"{base}__tree{idx+1}"
            out[gene] = vec

    print(f"[INFO] Parsed files: {n_files_parsed}, trees ok: {n_trees_ok}, with ref-overlap: {n_with_overlap}, genes kept: {len(out)}")
    return out

# ----------------------- PIC：独立对比 -----------------------

def repair_tree_for_pic(tr: Tree) -> Tree:
    """复制、修复后用于 PIC：正分支、二叉化、多点重根兜底。"""
    tr2 = copy.deepcopy(tr)
    ensure_positive_branch_lengths(tr2, EPS)
    # 二叉化（重要：否则无法逐节点对比）
    binarize_polytomies(tr2, EPS)
    # 参考树已在外部重根；这里不再变更根
    return tr2

def pic_contrasts(tr: Tree, values: Dict[str, float]) -> List[float]:
    """
    返回所有内部节点的对比列表。
    Felsenstein (1985) 实现：
      对叶子：x=value, v=0；
      对内部：v1=child1.branch+child1.v, v2=child2.branch+child2.v
             contrast = (x1-x2)/sqrt(v1+v2)
             x = (x1/v1 + x2/v2) / (1/v1 + 1/v2)
             v = 1/(1/v1 + 1/v2)
    """
    tr = repair_tree_for_pic(tr)
    # 每个 tip 必须有值
    for tip in tr.get_terminals():
        if tip.name not in values:
            raise ValueError(f"Missing value for species: {tip.name}")

    # 在节点属性上存值
    for tip in tr.get_terminals():
        tip._pic_x = float(values[tip.name])
        tip._pic_v = 0.0

    contrasts = []

    # 后序遍历，保证孩子先计算
    def postorder(cl: Clade):
        if len(cl.clades) == 0:
            return
        # 二叉已保证
        for ch in cl.clades:
            postorder(ch)
        if len(cl.clades) != 2:
            # 理论上不会发生；防御
            return
        c1, c2 = cl.clades
        bl1 = c1.branch_length if c1.branch_length is not None else EPS
        bl2 = c2.branch_length if c2.branch_length is not None else EPS
        v1 = bl1 + getattr(c1, "_pic_v", 0.0)
        v2 = bl2 + getattr(c2, "_pic_v", 0.0)
        # 避免除零
        v1 = max(v1, EPS)
        v2 = max(v2, EPS)
        x1 = getattr(c1, "_pic_x")
        x2 = getattr(c2, "_pic_x")
        # 该内部节点的对比
        denom = math.sqrt(v1 + v2)
        if denom <= 0:
            denom = math.sqrt(max(v1 + v2, EPS))
        contrasts.append((x1 - x2) / denom)
        # 向上回传
        w1 = 1.0 / v1
        w2 = 1.0 / v2
        cl._pic_x = (x1 * w1 + x2 * w2) / (w1 + w2)
        cl._pic_v = 1.0 / (w1 + w2)

    postorder(tr.root)
    return [float(c) for c in contrasts if np.isfinite(c)]

def pic_pearson_correlation(subtree: Tree, vec_x: Dict[str,float], vec_y: Dict[str,float]) -> Tuple[float,float]:
    """在子树上计算两基因的 PIC 相关（Pearson）。"""
    # 两个向量都要有所有 tips 的值
    tips = [t.name for t in subtree.get_terminals()]
    for sp in tips:
        if sp not in vec_x or sp not in vec_y:
            raise ValueError("PIC inputs must share identical species set.")
    cx = pic_contrasts(subtree, vec_x)
    cy = pic_contrasts(subtree, vec_y)
    m = min(len(cx), len(cy))
    if m < 3:
        return (float("nan"), float("nan"))
    r, p = stats.pearsonr(np.array(cx[:m]), np.array(cy[:m]))
    return (float(r), float(p))

# ----------------------- PGLS（Brownian 协方差） -----------------------

def tip_depths(tr: Tree) -> Dict[str, float]:
    d = tr.depths()
    out = {}
    for tip in tr.get_terminals():
        if tip in d and tip.name:
            out[tip.name] = float(d[tip])
    return out

def mrca_height(tr: Tree, name1: str, name2: str) -> float:
    # Biopython 的 common_ancestor 可用名字列表
    try:
        anc = tr.common_ancestor({"name": name1}, {"name": name2})
    except Exception:
        # 兜底：用 path 到根找最后公共节点
        path1 = tr.get_path(name1)
        path2 = tr.get_path(name2)
        common = tr.root
        for a, b in zip(path1, path2):
            if a is b:
                common = a
            else:
                break
        anc = common
    # 高度 = 根到该节点的路径长度
    total = 0.0
    # 从根走到 anc，利用 depths 更快：这里直接用 depths
    return float(tr.depths().get(anc, 0.0))

def brownian_covariance(tr: Tree, species: List[str]) -> np.ndarray:
    """C[i,j] = height(MRCA(i,j))，对角线 Var = height(tip)。"""
    ensure_positive_branch_lengths(tr, EPS)
    # 预计算 depth
    depth_map = tr.depths()
    # mrca() 对大量对儿较慢，这里用 Biopython 的 common_ancestor，必要时可缓存
    n = len(species)
    C = np.zeros((n, n), dtype=float)
    # 对角
    for i, sp in enumerate(species):
        tip = next(t for t in tr.get_terminals() if t.name == sp)
        C[i, i] = float(depth_map.get(tip, 0.0))
    # 非对角
    for i in range(n):
        for j in range(i+1, n):
            try:
                anc = tr.common_ancestor({"name": species[i]}, {"name": species[j]})
                h = float(depth_map.get(anc, 0.0))
            except Exception:
                h = 0.0
            C[i, j] = C[j, i] = h
    # 保证正定（最小抖动）
    # 若需要更强健，可做最近正定矩阵投影；这里先加极小对角
    C += np.eye(n) * EPS
    return C

def pgls_fit(subtree: Tree, vec_x: Dict[str,float], vec_y: Dict[str,float]) -> Tuple[float,float,float,float]:
    """
    在子树上做 y~1+x 的 GLS，协方差=Brownian C。
    返回：beta1, t, p, beta0
    """
    species = [t.name for t in subtree.get_terminals()]
    X = np.array([vec_x[sp] for sp in species], dtype=float)
    Y = np.array([vec_y[sp] for sp in species], dtype=float)

    # 去常量保护
    if np.allclose(X, X[0]) or np.allclose(Y, Y[0]):
        return (float("nan"), float("nan"), float("nan"), float("nan"))

    C = brownian_covariance(subtree, species)
    # Cholesky 白化
    try:
        L = np.linalg.cholesky(C)
    except LinAlgError:
        # 再加一点抖动后重试
        jitter = EPS
        for _ in range(5):
            try:
                L = np.linalg.cholesky(C + np.eye(len(C)) * jitter)
                break
            except LinAlgError:
                jitter *= 10
        else:
            return (float("nan"), float("nan"), float("nan"), float("nan"))

    Linv = np.linalg.inv(L)
    y_t = Linv @ Y
    X_t = Linv @ np.column_stack([np.ones_like(X), X])  # [intercept, x]

    # OLS on whitened
    XtX = X_t.T @ X_t
    try:
        XtX_inv = np.linalg.inv(XtX)
    except LinAlgError:
        return (float("nan"), float("nan"), float("nan"), float("nan"))
    beta = XtX_inv @ (X_t.T @ y_t)  # [b0, b1]
    yhat = X_t @ beta
    resid = y_t - yhat
    n, p = len(Y), 2
    dof = max(n - p, 1)
    sigma2 = float(resid.T @ resid) / dof
    cov_beta = XtX_inv * sigma2
    se_b1 = math.sqrt(max(cov_beta[1,1], EPS))
    t_b1 = float(beta[1] / se_b1)
    p_b1 = 2.0 * stats.t.sf(abs(t_b1), df=dof)
    return (float(beta[1]), t_b1, float(p_b1), float(beta[0]))

# ----------------------- ERC 主逻辑 -----------------------

def pairwise_all_metrics(ref_tree: Tree,
                         vec1: Dict[str,float], vec2: Dict[str,float],
                         min_species: int) -> Tuple[int, float,float, float,float, float,float,float,float]:
    shared = sorted(set(vec1).intersection(vec2))
    n = len(shared)
    if n < min_species:
        return (n, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    x = np.array([vec1[s] for s in shared], dtype=float)
    y = np.array([vec2[s] for s in shared], dtype=float)

    # Raw
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        r_raw, p_raw = (np.nan, np.nan)
        rho_raw, ps_raw = (np.nan, np.nan)
    else:
        r_raw, p_raw = stats.pearsonr(x, y)
        rho_raw, ps_raw = stats.spearmanr(x, y)

    # 子树（PIC & PGLS）
    sub = prune_to_species(ref_tree, set(shared))
    sub = repair_tree_for_pic(sub)

    # PIC
    try:
        r_pic, p_pic = pic_pearson_correlation(sub, {s: vec1[s] for s in shared},
                                                    {s: vec2[s] for s in shared})
    except Exception:
        r_pic, p_pic = (np.nan, np.nan)

    # PGLS
    try:
        b1, tval, pval, b0 = pgls_fit(sub, {s: vec1[s] for s in shared},
                                           {s: vec2[s] for s in shared})
    except Exception:
        b1 = tval = pval = b0 = np.nan

    return (n, float(r_raw), float(p_raw),
               float(rho_raw), float(ps_raw),
               float(r_pic), float(p_pic),
               float(b1), float(pval))

# ----------------------- 主程序 -----------------------

def main():
    ap = argparse.ArgumentParser(
        description="Compute ERC (Raw Pearson/Spearman, PIC-Pearson, PGLS) from gene trees using root-to-tip vectors. Robust tree repair to avoid NaNs."
    )
    ap.add_argument("--ref-tree", required=True, help="Reference species tree (Newick).")
    ap.add_argument("--mt-dir", required=False, default=None, help="Directory of mt gene trees.")
    ap.add_argument("--nu-dir", required=False, default=None, help="Directory of nu gene trees.")
    ap.add_argument("--out-prefix", required=True, help="Output file prefix.")
    ap.add_argument("--outgroup", default="17_SAY", help="Outgroup tip name for rerooting (default: 17_SAY).")
    ap.add_argument("--transform", choices=["none","sqrt","log1p"], default="sqrt", help="Transform branch lengths before ERC.")
    ap.add_argument("--min-species", type=int, default=15, help="Minimum shared species for a gene pair.")
    ap.add_argument("--between-only", action="store_true", help="Only compute mt–nu pairs (ignore mt–mt and nu–nu).")
    # 命名
    ap.add_argument("--name-map", default=None, help="TSV two columns: old_name  new_name.")
    ap.add_argument("--upper", action="store_true", help="Upper-case all tip names.")
    ap.add_argument("--strip-square", action="store_true", help="Strip surrounding [ ... ] in names.")
    ap.add_argument("--drop-prefix-number", action="store_true", help="Drop leading digits + optional '_'/'-' (e.g., '17_SAY'->'SAY').")

    args = ap.parse_args()

    name_map = load_name_map(args.name_map)
    norm_name_fn = build_norm_fn(
        strip_square=args.strip_square,
        drop_prefix_number=args.drop_prefix_number,
        upper=args.upper,
        name_map=name_map
    )

    # 读 & 规范化参考树
    ref = Phylo.read(args.ref_tree, "newick")
    normalize_tree_tip_names(ref, norm_name_fn)
    ensure_positive_branch_lengths(ref, EPS)
    reroot_by_outgroup_or_midpoint(ref, args.outgroup)
    ensure_positive_branch_lengths(ref, EPS)
    binarize_polytomies(ref, EPS)  # PIC 需要二叉；对参考树做一次全局处理

    ref_species = {t.name for t in ref.get_terminals() if t.name}
    if args.outgroup and args.outgroup not in ref_species:
        print(f"[WARN] Outgroup '{args.outgroup}' not found in normalized reference species set; will still try to outgroup-root gene trees if their tips contain it.")

    # 载入基因 RTT 向量
    gene_vecs = {}
    tag = {}
    if args.mt_dir:
        mt = load_gene_vectors(args.mt_dir, ref_species, args.outgroup, args.transform, norm_name_fn)
        for g,v in mt.items(): gene_vecs[g]=v; tag[g]="mt"
    if args.nu_dir:
        nu = load_gene_vectors(args.nu_dir, ref_species, args.outgroup, args.transform, norm_name_fn)
        for g,v in nu.items(): gene_vecs[g]=v; tag[g]="nu"

    genes = sorted(gene_vecs.keys())
    print(f"[INFO] Total genes loaded: {len(genes)}")
    if len(genes) < 2:
        sys.exit("ERROR: fewer than 2 genes loaded.")

    # 导出 gene×species 矩阵（QC）
    species_union = sorted(set().union(*[set(v.keys()) for v in gene_vecs.values()]))
    mat_path = args.out_prefix + ".gene_by_species_matrix.tsv"
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    with open(mat_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(["gene"] + species_union)
        for g in genes:
            row = [g] + [gene_vecs[g].get(sp, "") for sp in species_union]
            w.writerow(row)

    # 配对计算
    out_rows = []
    for g1, g2 in combinations(genes, 2):
        if args.between_only:
            if not ((tag.get(g1)=="mt" and tag.get(g2)=="nu") or (tag.get(g1)=="nu" and tag.get(g2)=="mt")):
                continue
        (n, r_raw, p_raw, rho_raw, ps_raw, r_pic, p_pic, b1, p_b1) = \
            pairwise_all_metrics(ref, gene_vecs[g1], gene_vecs[g2], args.min_species)
        if n >= args.min_species:
            out_rows.append([g1, g2, tag.get(g1,""), tag.get(g2,""), n,
                             r_raw, p_raw, rho_raw, ps_raw, r_pic, p_pic, b1, p_b1])

    out_path = args.out_prefix + ".erc_pairs.csv"
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene1","gene2","tag1","tag2","n_shared",
                    "raw_pearson_r","raw_pearson_p",
                    "raw_spearman_rho","raw_spearman_p",
                    "pic_pearson_r","pic_pearson_p",
                    "pgls_beta1","pgls_p"])
        w.writerows(out_rows)

    kept = sum(1 for r in out_rows if np.isfinite(r[5]))
    print(f"[INFO] Pairs kept (n_shared >= {args.min_species}): {kept}")
    print(f"[OK] Outputs:\n  - {out_path}\n  - {mat_path}")

if __name__ == "__main__":
    main()


python erc_from_trees2.py \
  --ref-tree /mnt/spareHD_2/oxphos_gene_tree/species_astral.tre \
  --mt-dir   /mnt/spareHD_2/mt_gene_tree/pruned_trees_te \
  --nu-dir   /mnt/spareHD_2/oxphos_gene_tree/pruned_trees_te \
  --between-only \
  --out-prefix /mnt/spareHD_2/erc_out/te_mt_nu_pic_pgls_fix2 \
  --outgroup 17_SAY \
  --transform sqrt \
  --min-species 20 \
  --upper --drop-prefix-number




#fdr
#!/usr/bin/env python3
import pandas as pd, numpy as np

IN  = "/mnt/spareHD_2/erc_out/te_mt_nu_pic_pgls_fix2.erc_pairs.csv"
OUT = "/mnt/spareHD_2/erc_out/te_mt_nu_triage_from_fix2"  # 输出前缀

df = pd.read_csv(IN)

# 只看 mt–nu 配对
df = df[(df["tag1"]=="mt") & (df["tag2"]=="nu")].copy()

# --- BH-FDR 计算 ---
def bh(p):
    x = np.array(p, dtype=float)
    x[np.isnan(x)] = 1.0
    order = np.argsort(x)
    ranked = np.empty_like(x, dtype=float)
    m = len(x)
    ranked[order] = x[order] * m / (np.arange(m)+1)
    return np.minimum.accumulate(ranked[::-1])[::-1]

# 与 CSV 列名一一对应
df["FDR_raw"]  = bh(df["raw_pearson_p"].values)   if "raw_pearson_p"  in df else np.nan
df["FDR_pic"]  = bh(df["pic_pearson_p"].values)   if "pic_pearson_p"  in df else np.nan
df["FDR_pgls"] = bh(df["pgls_p"].values)          if "pgls_p"         in df else np.nan

# --- 三角判据打标签 ---
def label(row, alpha=0.05, near=0.05):
    r_raw = row.get("raw_pearson_r", np.nan)
    r_pic = row.get("pic_pearson_r", np.nan)
    b1    = row.get("pgls_beta1",    np.nan)

    FDR_raw  = row.get("FDR_raw",  np.nan)
    FDR_pic  = row.get("FDR_pic",  np.nan)
    FDR_pgls = row.get("FDR_pgls", np.nan)

    # 缺失时的安全默认
    if np.isnan(r_pic): r_pic = 0.0
    if np.isnan(b1):   b1    = 0.0

    any_sig = any([
        (FDR_raw  < alpha) if not np.isnan(FDR_raw)  else False,
        (FDR_pic  < alpha) if not np.isnan(FDR_pic)  else False,
        (FDR_pgls < alpha) if not np.isnan(FDR_pgls) else False,
    ])

    # 1) Raw<0 且 PIC≈0 且 PGLS≈0 且都不显著 → 结构性负号
    if (r_raw < 0) and (abs(r_pic) < near) and (abs(b1) < near) and (not any_sig):
        return "likely_structural_negative"

    # 2) Raw<0、PIC<0、PGLS<0 且 ≥1 显著 → 反向耦合候选
    if (r_raw < 0) and (r_pic < 0) and (b1 < 0) and any_sig:
        return "candidate_negative_coupling"

    # 3) Raw<0，但 PIC>0 或 PGLS>0 → 结构性差异，校正后为正
    if (r_raw < 0) and ((r_pic > 0) or (b1 > 0)):
        return "negative_raw_but_corrected_positive"

    return "other"

df["triage"] = df.apply(label, axis=1)

# --- 导出 ---
df.to_csv(OUT+"_all.csv", index=False)
df[df["triage"]=="candidate_negative_coupling"].to_csv(OUT+"_neg_coupling.csv", index=False)
df[df["triage"]=="likely_structural_negative"].to_csv(OUT+"_likely_structural.csv", index=False)
df[df["triage"]=="negative_raw_but_corrected_positive"].to_csv(OUT+"_neg_raw_pos_corrected.csv", index=False)

print("Done. Outputs:")
print(" ", OUT+"_all.csv")
print(" ", OUT+"_neg_coupling.csv")
print(" ", OUT+"_likely_structural.csv")
print(" ", OUT+"_neg_raw_pos_corrected.csv")


