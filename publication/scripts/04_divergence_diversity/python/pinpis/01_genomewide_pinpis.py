#!/usr/bin/env python3
"""Genome-wide pooled piN/piS from a PoPoolation2 sync file.

Main corrections relative to the original targeted-gene implementation:
1. Complement alleles for CDS features on the negative strand.
2. Normalize piN and piS by nonsynonymous and synonymous site opportunities.
3. Read the large sync once and construct/release one chromosome index at a time.
4. Use an explicit sync-column map; excluded samples are never analyzed.
"""

import argparse
import gzip
import math
import re
import sys
from collections import defaultdict

import pysam

CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}
BASES = ("A", "C", "G", "T")
COMP = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}


def op(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def attrs(text):
    out = dict(re.findall(r'(\S+)\s+"([^"]+)"', text))
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.strip().split("=", 1)
            out.setdefault(key, value)
    return out


def read_column_map(path):
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip().split("\t")
        ix = {x: i for i, x in enumerate(header)}
        required = {"sync_column", "population", "include", "n_individuals"}
        if not required.issubset(ix):
            raise ValueError(f"Column map requires {sorted(required)}")
        for line in handle:
            f = line.rstrip().split("\t")
            rows.append({
                "column": int(f[ix["sync_column"]]),
                "population": f[ix["population"]],
                "include": f[ix["include"]].lower() == "yes",
                "n": int(f[ix["n_individuals"]]),
                "max_dp": int(f[ix["max_depth"]]) if "max_depth" in ix and f[ix["max_depth"]] else 0,
            })
    rows.sort(key=lambda x: x["column"])
    if [x["column"] for x in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("sync_column must be consecutive and start at 1")
    return rows


def read_models(path):
    tmp = defaultdict(lambda: {"strand": None, "chrom": None, "blocks": []})
    with op(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 9 or f[2] != "CDS":
                continue
            a = attrs(f[8])
            gene, tx = a.get("gene_id"), a.get("transcript_id")
            if not gene or not tx:
                continue
            key = (gene, tx)
            tmp[key]["strand"] = f[6]
            tmp[key]["chrom"] = f[0]
            tmp[key]["blocks"].append((int(f[3]), int(f[4])))
    models = {}
    for (gene, tx), x in tmp.items():
        if gene in models:
            raise ValueError(f"Multiple selected transcripts for {gene}")
        x.update(gene=gene, transcript=tx)
        models[gene] = x
    return models


def read_classes(path):
    out = {}
    with open(path) as handle:
        h = handle.readline().rstrip().split("\t")
        ix = {x: i for i, x in enumerate(h)}
        gi = ix.get("gene_id", 0)
        ci = next((ix[x] for x in ("primary_class", "Kuster_class", "class") if x in ix), None)
        si = next((ix[x] for x in ("symbol", "gene_symbol", "new_symbol", "gene_name") if x in ix), None)
        if ci is None:
            raise ValueError("Cannot find classification column")
        for line in handle:
            f = line.rstrip().split("\t")
            gene = f[gi]
            out[gene] = (f[ci], f[si] if si is not None and si < len(f) and f[si] else gene)
    return out


def coding_positions(model):
    blocks = sorted(model["blocks"], reverse=model["strand"] == "-")
    positions = []
    for start, end in blocks:
        positions.extend(range(start, end + 1) if model["strand"] == "+" else range(end, start - 1, -1))
    return positions


def aa(codon):
    return CODE.get(codon) if len(codon) == 3 and set(codon) <= set(BASES) else None


def opportunities(codon, pos):
    refaa = aa(codon)
    if refaa in (None, "*"):
        return None
    n = s = 0.0
    for alt in BASES:
        if alt == codon[pos]:
            continue
        c = list(codon); c[pos] = alt
        altaa = aa("".join(c))
        if altaa in (None, "*"):
            continue
        if altaa == refaa:
            s += 1.0 / 3.0
        else:
            n += 1.0 / 3.0
    return n, s


def build_chrom_index(chrom, genes, models, fasta):
    seq = fasta.fetch(chrom).upper()
    posmap = defaultdict(list)
    lengths = {}
    for gene in genes:
        model = models[gene]
        positions = coding_positions(model)
        lengths[gene] = len(positions)
        if len(positions) % 3:
            print(f"[skip non-triplet] {gene}\t{len(positions)}", file=sys.stderr)
            continue
        for start in range(0, len(positions), 3):
            gp = positions[start:start + 3]
            genomic = [seq[p - 1] for p in gp]
            codon = "".join(genomic if model["strand"] == "+" else [COMP.get(x, "N") for x in genomic])
            if aa(codon) in (None, "*"):
                continue
            for cp, p in enumerate(gp):
                opp = opportunities(codon, cp)
                if opp:
                    posmap[p].append((gene, codon, cp, model["strand"], opp[0], opp[1]))
    print(f"[index] {chrom}: genes={len(genes)} positions={len(posmap)}", file=sys.stderr)
    return posmap, lengths


def counts(field, minus):
    try:
        x = [int(v) for v in field.split(":")[:4]]
    except (ValueError, IndexError):
        return None
    genomic = {"A":x[0], "T":x[1], "C":x[2], "G":x[3]}
    return ({"A":genomic["T"], "C":genomic["G"], "G":genomic["C"], "T":genomic["A"]}
            if minus else {b:genomic[b] for b in BASES})


def diversity(cnt, codon, pos, correction):
    cov = sum(cnt.values())
    freq = {b: cnt[b] / cov for b in BASES}
    pn = ps = 0.0
    for i, b1 in enumerate(BASES):
        for b2 in BASES[i + 1:]:
            if not freq[b1] or not freq[b2]:
                continue
            c1 = list(codon); c1[pos] = b1
            c2 = list(codon); c2[pos] = b2
            a1, a2 = aa("".join(c1)), aa("".join(c2))
            if a1 in (None, "*") or a2 in (None, "*"):
                continue
            value = 2 * freq[b1] * freq[b2] * correction
            if a1 == a2: ps += value
            else: pn += value
    return pn, ps


def fmt(x):
    return "NA" if x is None or not math.isfinite(x) else f"{x:.10g}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sync", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--column-map", required=True)
    p.add_argument("--classification", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-coverage", type=int, default=10)
    p.add_argument("--min-callable-fraction", type=float, default=0.70)
    p.add_argument("--only-chromosome")
    args = p.parse_args()

    columns = read_column_map(args.column_map)
    selected = [(i, x) for i, x in enumerate(columns) if x["include"]]
    models = read_models(args.gtf)
    classes = read_classes(args.classification)
    fasta = pysam.FastaFile(args.reference)
    chrom_genes = defaultdict(list)
    for gene, model in models.items():
        chrom_genes[model["chrom"]].append(gene)

    acc = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0]))
    cds_lengths = {}
    current_chrom = None
    target_chromosome_seen = False
    posmap = None
    completed_chromosomes = set()
    rows = targets = 0

    with op(args.sync) as handle:
        for line in handle:
            rows += 1
            f = line.rstrip().split("\t")
            if len(f) != 3 + len(columns):
                raise ValueError(f"sync row {rows}: expected {len(columns)} samples, observed {len(f)-3}")
            chrom = f[0]
            if args.only_chromosome:
                if chrom == args.only_chromosome:
                    target_chromosome_seen = True
                elif target_chromosome_seen:
                    break
                else:
                    continue
            if chrom != current_chrom:
                if chrom in completed_chromosomes:
                    raise ValueError(f"Sync is not chromosome-grouped: {chrom} occurs again")
                if current_chrom is not None:
                    completed_chromosomes.add(current_chrom)
                current_chrom = chrom
                if chrom in chrom_genes:
                    posmap, lens = build_chrom_index(chrom, chrom_genes[chrom], models, fasta)
                    cds_lengths.update(lens)
                else:
                    posmap = None
            if not posmap:
                continue
            try: pos = int(f[1])
            except ValueError: continue
            records = posmap.get(pos)
            if not records:
                continue
            targets += 1
            for gene, codon, cp, strand, nopp, sopp in records:
                for field_index, meta in selected:
                    cnt = counts(f[3 + field_index], strand == "-")
                    if cnt is None:
                        continue
                    dp = sum(cnt.values())
                    if dp < args.min_coverage or (meta["max_dp"] and dp > meta["max_dp"]):
                        continue
                    corr = (2 * meta["n"]) / (2 * meta["n"] - 1.0)
                    pn, ps = diversity(cnt, codon, cp, corr)
                    a = acc[gene][meta["population"]]
                    a[0] += pn; a[1] += ps; a[2] += nopp; a[3] += sopp; a[4] += 1
            if rows % 10_000_000 == 0:
                print(f"[stream] rows={rows:,} target_rows={targets:,}", file=sys.stderr)

    header = ["gene_id","symbol","Kuster_class","population","piN","piS","piN_piS",
              "piN_sum","piS_sum","N_opportunities","S_opportunities","callable_sites",
              "CDS_length","callable_fraction","pinpis_eligible","ratio_status"]
    with open(args.output, "w") as out:
        out.write("\t".join(header) + "\n")
        for gene in sorted(models):
            cls, symbol = classes.get(gene, ("unclassified", gene))
            length = cds_lengths.get(gene, len(coding_positions(models[gene])))
            for _, meta in selected:
                pn_sum, ps_sum, no, so, callable_n = acc[gene][meta["population"]]
                pn = pn_sum / no if no > 0 else None
                ps = ps_sum / so if so > 0 else None
                if pn is None or ps is None: ratio, status = None, "missing_opportunities"
                elif ps > 0: ratio, status = pn / ps, "finite"
                elif pn == 0: ratio, status = None, "piN0_piS0"
                else: ratio, status = None, "piNpositive_piS0"
                frac = callable_n / length if length else 0
                out.write("\t".join(map(str, [gene,symbol,cls,meta["population"],fmt(pn),fmt(ps),fmt(ratio),
                    fmt(pn_sum),fmt(ps_sum),fmt(no),fmt(so),callable_n,length,fmt(frac),
                    "yes" if frac >= args.min_callable_fraction else "no",status])) + "\n")
    print(f"[complete] rows={rows:,} targets={targets:,} output={args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
