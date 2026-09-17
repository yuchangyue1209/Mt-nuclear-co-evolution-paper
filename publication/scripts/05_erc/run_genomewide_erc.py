#!/usr/bin/env python3
"""Genome-wide mitonuclear ERC from fixed-topology codeml M0 outputs.

The script requires nuclear and mitochondrial genes to have been fitted on the
same 27-population topology. Branches are matched by PAML's parent..child labels.
Relative rates are residuals from a through-origin regression of each gene's
branch lengths on the genome-wide mean branch-length vector. This avoids the
zero-median problem caused by sparse within-species codeml branch changes. The
primary mitochondrial predictor is the median standardized residual across all
13 mt protein-coding genes.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


MT13 = {
    "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
}

MT_COMPLEX = {
    "ND1": "Complex I", "ND2": "Complex I", "ND3": "Complex I",
    "ND4": "Complex I", "ND4L": "Complex I", "ND5": "Complex I",
    "ND6": "Complex I", "CYTB": "Complex III",
    "COX1": "Complex IV", "COX2": "Complex IV", "COX3": "Complex IV",
    "ATP6": "Complex V", "ATP8": "Complex V",
}

BRANCH_HEADER = re.compile(r"^\s*branch\s+t\s+N\s+S\s+dN/dS", re.I)
BRANCH_ROW = re.compile(
    r"^\s*(\d+\.\.\d+)\s+"
    r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+"
    r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)"
)


def canonical_mt(name: str) -> str:
    x = name.upper().replace("-", "").replace("_", "")
    aliases = {"COI": "COX1", "COXI": "COX1", "CO1": "COX1",
               "COII": "COX2", "COXII": "COX2", "CO2": "COX2",
               "COIII": "COX3", "COXIII": "COX3", "CO3": "COX3"}
    return aliases.get(x, x)


def read_gene_list(path: Path) -> list[str]:
    genes = [x.strip() for x in path.read_text().splitlines() if x.strip()]
    if len(genes) != len(set(genes)):
        raise ValueError(f"Duplicate genes in {path}")
    return genes


def parse_mlc(path: Path, metric: str) -> dict[str, float]:
    """Return branch_label -> t, dN, or dS from the M0 branch table."""
    values: dict[str, float] = {}
    in_table = False
    with path.open(errors="replace") as handle:
        for line in handle:
            if BRANCH_HEADER.search(line):
                in_table = True
                continue
            if not in_table:
                continue
            m = BRANCH_ROW.search(line)
            if m:
                branch = m.group(1)
                t, dn, ds = float(m.group(2)), float(m.group(6)), float(m.group(7))
                values[branch] = {"t": t, "dN": dn, "dS": ds}[metric]
            elif values and not line.strip():
                break
    if not values:
        raise ValueError(f"No codeml branch table parsed from {path}")
    return values


def locate_mlc(root: Path, gene: str, suffix: str) -> Path:
    candidates = [
        root / gene / suffix,
        root / gene.upper() / suffix,
        root / gene.lower() / suffix,
        root / f"{gene}_M0" / suffix,
        root / f"{gene.upper()}_M0" / suffix,
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f"Cannot locate {suffix} for {gene} under {root}")


def load_matrix(root: Path, genes: list[str], metric: str, suffix: str) -> pd.DataFrame:
    rows = {}
    failures = []
    for i, gene in enumerate(genes, 1):
        try:
            rows[gene] = parse_mlc(locate_mlc(root, gene, suffix), metric)
        except Exception as exc:
            failures.append((gene, str(exc)))
        if i % 1000 == 0:
            print(f"[load] {i}/{len(genes)}", flush=True)
    if failures:
        preview = "\n".join(f"  {g}: {e}" for g, e in failures[:10])
        raise RuntimeError(f"Failed to load {len(failures)} genes:\n{preview}")
    matrix = pd.DataFrame.from_dict(rows, orient="index").sort_index(axis=1)
    if matrix.isna().any().any():
        bad = int(matrix.isna().sum().sum())
        raise ValueError(f"Branch matrix contains {bad} missing cells; topology/labels differ")
    return matrix


def bh_adjust(pvalues: pd.Series) -> np.ndarray:
    p = pd.to_numeric(pvalues, errors="coerce").fillna(1.0).to_numpy(float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, float)
    out[order] = np.minimum(ranked, 1.0)
    return out


def regression_relative_rates(raw: pd.DataFrame, background: pd.Series) -> pd.DataFrame:
    """Return standardized residual branch rates for every gene.

    For gene g, fit t_g,b = beta_g * background_b through the origin. Residuals
    remove both the branch-wide opportunity for change and the gene's overall
    rate. Each residual vector is then centered and RMS-scaled; scaling does not
    change correlations but makes the mt composite give genes comparable weight.
    """
    bg = background.to_numpy(float)
    denom = float(np.dot(bg, bg))
    if not np.isfinite(denom) or denom <= 0:
        raise ValueError("Genome-wide branch background has no positive information")
    out = np.empty(raw.shape, dtype=float)
    for i, (_, row) in enumerate(raw.iterrows()):
        y = row.to_numpy(float)
        beta = max(0.0, float(np.dot(y, bg) / denom))
        resid = y - beta * bg
        resid -= np.mean(resid)
        scale = float(np.sqrt(np.mean(resid * resid)))
        out[i, :] = resid / scale if scale > 0 else 0.0
    return pd.DataFrame(out, index=raw.index, columns=raw.columns)


def safe_cor(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float]:
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 6 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return math.nan, math.nan
    result = stats.pearsonr(x, y) if method == "pearson" else stats.spearmanr(x, y)
    return float(result.statistic), float(result.pvalue)


def leave_one_branch_out(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float, float]:
    vals = []
    for i in range(len(x)):
        r, _ = safe_cor(np.delete(x, i), np.delete(y, i), method)
        if np.isfinite(r):
            vals.append(r)
    if not vals:
        return math.nan, math.nan, math.nan
    return float(np.min(vals)), float(np.median(vals)), float(np.max(vals))


def correlate_predictor(name: str, predictor: pd.Series, nuclear_rer: pd.DataFrame,
                        method: str, do_loo: bool) -> pd.DataFrame:
    y = predictor.to_numpy(float)
    rows = []
    for gene, row in nuclear_rer.iterrows():
        x = row.to_numpy(float)
        r, p = safe_cor(y, x, method)
        record = {"mt_predictor": name, "nuclear_gene": gene,
                  "n_branches": int(np.isfinite(x * y).sum()), "r": r, "p": p}
        if do_loo:
            lo, med, hi = leave_one_branch_out(y, x, method)
            record.update({"loo_min_r": lo, "loo_median_r": med, "loo_max_r": hi})
        rows.append(record)
    out = pd.DataFrame(rows)
    out["padj"] = bh_adjust(out["p"])
    return out


def add_metadata(results: pd.DataFrame, metadata_path: Path | None) -> pd.DataFrame:
    if metadata_path is None:
        return results
    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)
    candidates = ["gene_id", "stickleback_gene_id", "stickleback_gene",
                  "ensembl_gene_id", "Gene stable ID", "gene"]
    id_column = next((x for x in candidates if x in meta.columns), None)
    if id_column is None:
        raise ValueError(f"Cannot identify metadata gene-ID column: {list(meta.columns)}")
    print(f"[metadata] Using gene-ID column: {id_column}", flush=True)
    meta = meta.drop_duplicates(id_column)
    merged = results.merge(meta, how="left", left_on="nuclear_gene", right_on=id_column)
    print(f"[metadata] Matched {merged[id_column].notna().sum()}/{len(merged)} ERC rows", flush=True)
    return merged


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nuclear-results", type=Path, required=True,
                    help="Directory containing one <gene>/mlc per nuclear gene")
    ap.add_argument("--nuclear-genes", type=Path, required=True,
                    help="One variable nuclear gene ID per line")
    ap.add_argument("--mt-results", type=Path, required=True,
                    help="Directory containing one <mtgene>/mlc per mtPCG")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, default=None,
                    help="Optional gene metadata/classification TSV with gene_id")
    ap.add_argument("--metric", choices=["t", "dN", "dS"], default="t")
    ap.add_argument("--method", choices=["pearson", "spearman"], default="spearman")
    ap.add_argument("--mlc-name", default="mlc")
    ap.add_argument("--expected-nuclear", type=int, default=15241)
    ap.add_argument("--leave-one-branch-out", action="store_true")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    nuclear_genes = read_gene_list(args.nuclear_genes)
    if len(nuclear_genes) != args.expected_nuclear:
        raise SystemExit(f"Expected {args.expected_nuclear} variable nuclear genes; found {len(nuclear_genes)}")

    print("[1/6] Loading nuclear codeml branch lengths", flush=True)
    nuclear_raw = load_matrix(args.nuclear_results, nuclear_genes, args.metric, args.mlc_name)
    # Include zeros: the mean across all variable genes estimates the total
    # genome-wide opportunity/rate for each branch. The median is inappropriate
    # here because most individual branches have zero change in >50% of genes.
    background = nuclear_raw.mean(axis=0)
    if (background <= 0).any():
        bad = background[background <= 0].index.tolist()
        raise SystemExit(f"Non-positive genome-wide background for branches: {bad}")

    print("[2/6] Loading the 13 mitochondrial PCGs", flush=True)
    mt_raw = load_matrix(args.mt_results, sorted(MT13), args.metric, args.mlc_name)
    got = {canonical_mt(x) for x in mt_raw.index}
    if got != MT13:
        raise SystemExit(f"mtPCG set mismatch; missing={sorted(MT13-got)}, extra={sorted(got-MT13)}")
    if list(mt_raw.columns) != list(nuclear_raw.columns):
        raise SystemExit("Nuclear and mitochondrial PAML branch labels do not match exactly")

    print("[3/6] Calculating centered relative evolutionary rates", flush=True)
    nuclear_rer = regression_relative_rates(nuclear_raw, background)
    mt_rer = regression_relative_rates(mt_raw, background)

    predictors = {"mtPCG_composite": mt_rer.median(axis=0)}
    for gene in sorted(MT13):
        predictors[gene] = mt_rer.loc[gene]
    for complex_name in sorted(set(MT_COMPLEX.values())):
        members = [g for g, c in MT_COMPLEX.items() if c == complex_name]
        predictors[complex_name.replace(" ", "_")] = mt_rer.loc[members].median(axis=0)

    print("[4/6] Calculating genome-wide ERC", flush=True)
    results = pd.concat([
        correlate_predictor(name, pred, nuclear_rer, args.method, args.leave_one_branch_out)
        for name, pred in predictors.items()
    ], ignore_index=True)
    results = add_metadata(results, args.metadata)
    results["metric"] = args.metric
    results["correlation_method"] = args.method
    results = results.sort_values(["mt_predictor", "padj", "r"], ascending=[True, True, False])

    print("[5/6] Writing result and QC tables", flush=True)
    nuclear_raw.to_csv(args.out_dir / "nuclear_branch_lengths.tsv.gz", sep="\t", compression="gzip")
    nuclear_rer.to_csv(args.out_dir / "nuclear_relative_rates.tsv.gz", sep="\t", compression="gzip")
    mt_raw.to_csv(args.out_dir / "mt13_branch_lengths.tsv", sep="\t")
    mt_rer.to_csv(args.out_dir / "mt13_relative_rates.tsv", sep="\t")
    background.rename("genomewide_mean_branch_length").to_csv(
        args.out_dir / "genomewide_branch_background.tsv", sep="\t", header=True)
    results.to_csv(args.out_dir / "erc_all_predictors.tsv.gz", sep="\t", index=False, compression="gzip")
    results[results.mt_predictor == "mtPCG_composite"].to_csv(
        args.out_dir / "erc_mtPCG_composite.tsv", sep="\t", index=False)
    results[(results.mt_predictor == "mtPCG_composite") & (results.padj < 0.05)].to_csv(
        args.out_dir / "erc_mtPCG_composite_FDR05.tsv", sep="\t", index=False)

    with (args.out_dir / "run_summary.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["item", "value"])
        w.writerow(["nuclear_genes", nuclear_raw.shape[0]])
        w.writerow(["mt_genes", mt_raw.shape[0]])
        w.writerow(["branches", nuclear_raw.shape[1]])
        w.writerow(["predictors", len(predictors)])
        w.writerow(["metric", args.metric])
        w.writerow(["correlation_method", args.method])

    print("[6/6] Complete", flush=True)
    print(f"[OK] Results: {args.out_dir}")


if __name__ == "__main__":
    main()
