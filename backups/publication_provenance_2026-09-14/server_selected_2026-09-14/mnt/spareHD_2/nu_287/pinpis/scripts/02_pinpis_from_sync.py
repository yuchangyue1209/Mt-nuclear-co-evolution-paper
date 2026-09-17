#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

# -----------------------------
# Genetic code (standard + mitochondrial vertebrate)
# -----------------------------
STD_CODE = {
    # U(T) first
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    # C first
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    # A first
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    # G first
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

# Vertebrate mitochondrial (translation table 2): AGA/AGG=Stop, ATA=M, TGA=W
MITO2_CODE = dict(STD_CODE)
MITO2_CODE.update({
    "ATA":"M",
    "TGA":"W",
    "AGA":"*",
    "AGG":"*",
})

DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


# -----------------------------
# FASTA fetcher (pysam)
# -----------------------------
class FastaFetcher:
    def __init__(self, fasta_path: str):
        try:
            import pysam  # type: ignore
        except Exception as e:
            raise RuntimeError(
                "ERROR: pysam is required for random access FASTA fetch.\n"
                "Install in your env, e.g.:\n"
                "  conda install -c bioconda pysam\n"
            ) from e
        self.pysam = pysam
        self.fasta_path = fasta_path
        self.fa = pysam.FastaFile(fasta_path)

    def fetch_base_1based(self, chrom: str, pos1: int) -> str:
        # pysam uses 0-based half-open
        return self.fa.fetch(chrom, pos1 - 1, pos1).upper()

    def fetch_seq_1based(self, chrom: str, start1: int, end1: int) -> str:
        # inclusive end
        return self.fa.fetch(chrom, start1 - 1, end1).upper()


# -----------------------------
# Parse pool sizes
# -----------------------------
def read_pool_sizes(pool_tsv: str):
    """
    Expect:
      Population <tab> n_individuals
    Returns dict pop->n_individuals (int)
    """
    d = {}
    with open(pool_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.lower().startswith("population"):
                continue
            parts = re.split(r"\t+", line)
            if len(parts) < 2:
                continue
            pop = parts[0].strip()
            try:
                n = int(float(parts[1]))
            except:
                continue
            d[pop] = n
    if not d:
        raise ValueError(f"Empty/invalid pool sizes file: {pool_tsv}")
    return d


# -----------------------------
# Parse your params format:
#   GENE aars2
#   STRAND +
#   BLOCK chrIX 2045.. 2045.. phase
# -----------------------------
def parse_params_blocks(params_path):
    """
    Returns:
      gene (str|None)
      strand ('+'|'-')
      blocks: list of (chrom, start, end, phase)   (1-based inclusive)
    """
    gene = None
    strand = None
    blocks = []

    with open(params_path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if not parts:
                continue
            key = parts[0]

            if key == "GENE" and len(parts) >= 2:
                gene = parts[1]
                continue

            if key == "STRAND" and len(parts) >= 2:
                strand = parts[1]
                continue

            if key == "BLOCK" and len(parts) >= 5:
                chrom = parts[1]
                start = int(parts[2])
                end = int(parts[3])
                phase = int(parts[4])
                blocks.append((chrom, start, end, phase))
                continue

            # tolerate STRAND=+ or STRAND:+
            if key.startswith("STRAND") and len(parts) == 1:
                if "=" in key:
                    strand = key.split("=", 1)[1].strip()
                elif ":" in key:
                    strand = key.split(":", 1)[1].strip()

    if strand is None:
        strand = "+"

    if not blocks:
        raise ValueError(f"Could not parse CDS blocks from params: {params_path}")

    # Keep as listed; but sort by genomic coordinate for determinism
    blocks = sorted(blocks, key=lambda x: (x[0], x[1], x[2]))
    return gene, strand, blocks


# -----------------------------
# Build CDS position map using blocks + phase
# We create an ordered list of genomic positions in coding order,
# trimming phase bases from the START of each CDS feature (GFF-like).
# -----------------------------
def build_cds_position_map(blocks, strand):
    """
    blocks: list of (chrom, start, end, phase), 1-based inclusive
    strand: '+'|'-'
    Returns:
      pos2cdsidx: dict[(chrom,pos1)] -> cds_index (0-based in CDS)
      codon2gpos: dict[codon_index] -> list of 3 genomic positions [(chrom,pos1),...]
    """

    # For '-' strand, coding order is reverse on genome; also phase trimming applies to the
    # 5' end of the CDS on the coding strand:
    #   '+' : trim from start => start += phase
    #   '-' : trim from end   => end   -= phase
    if strand not in ("+", "-"):
        strand = "+"

    ordered = []

    if strand == "+":
        for chrom, start, end, phase in blocks:
            s = start + phase
            e = end
            if s > e:
                continue
            for p in range(s, e + 1):
                ordered.append((chrom, p))
    else:
        # Reverse blocks in coding order: higher->lower by coordinate within each chrom
        blocks_rev = sorted(blocks, key=lambda x: (x[0], x[1], x[2]), reverse=True)
        for chrom, start, end, phase in blocks_rev:
            s = start
            e = end - phase
            if s > e:
                continue
            # coding order is from e down to s
            for p in range(e, s - 1, -1):
                ordered.append((chrom, p))

    pos2cdsidx = {}
    for i, gpos in enumerate(ordered):
        pos2cdsidx[gpos] = i

    codon2gpos = {}
    ncod = len(ordered) // 3
    for ci in range(ncod):
        codon2gpos[ci] = [ordered[3 * ci], ordered[3 * ci + 1], ordered[3 * ci + 2]]

    return pos2cdsidx, codon2gpos


# -----------------------------
# Parse sync counts
# sync columns typical:
#   chr  pos  ref  A:T:C:G:N:del   A:T:C:G:N:del  ...
# We'll use A,C,G,T only; ignore N/del.
# -----------------------------
BASES = ["A", "C", "G", "T"]
SYNC_ORDER = ["A", "T", "C", "G", "N", "DEL"]  # common PoPoolation order


def parse_sync_counts_field(field: str):
    parts = field.split(":")
    if len(parts) < 4:
        return None
    # Expect 6, but allow longer
    try:
        nums = [int(x) for x in parts[:6]]
    except:
        return None

    # Map to A,C,G,T using known order A:T:C:G
    a = nums[0]
    t = nums[1]
    c = nums[2]
    g = nums[3]
    return {"A": a, "C": c, "G": g, "T": t}


def translate_codon(codon: str, code_map):
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3 or re.search(r"[^ACGT]", codon):
        return "X"
    return code_map.get(codon, "X")


# -----------------------------
# π decomposition into syn/nonsyn at a site
# We approximate amino acid state for each allele at this position using:
#   codon = ref_codon with this base replaced (other 2 positions fixed as ref)
# Then classify allele pairs as syn/nonsyn based on aa equality.
# -----------------------------
def calc_site_pin_pis(freqs, ref_codon, pos_in_codon, code_map):
    """
    freqs: dict base->p  (A/C/G/T)
    ref_codon: 3-mer in coding orientation
    pos_in_codon: 0,1,2 within codon
    Returns (piN_raw, piS_raw) where raw = sum_{i<j} 2 p_i p_j classified
    """
    aa_of = {}
    for b in BASES:
        cod = list(ref_codon)
        cod[pos_in_codon] = b
        aa_of[b] = translate_codon("".join(cod), code_map)

    # pairwise contributions
    piN = 0.0
    piS = 0.0
    for i in range(4):
        bi = BASES[i]
        pi = freqs.get(bi, 0.0)
        if pi <= 0:
            continue
        for j in range(i + 1, 4):
            bj = BASES[j]
            pj = freqs.get(bj, 0.0)
            if pj <= 0:
                continue
            contrib = 2.0 * pi * pj
            if aa_of[bi] == aa_of[bj]:
                piS += contrib
            else:
                piN += contrib
    return piN, piS


def main():
    ap = argparse.ArgumentParser(
        description="Compute piN, piS, piN/piS from PoPoolation sync + custom .sync.params (BLOCK format)."
    )
    ap.add_argument("--sync-dir", required=True, help="Directory containing *.sync (per gene).")
    ap.add_argument("--pool-sizes", required=True, help="TSV: Population\\t n_individuals")
    ap.add_argument("--ref-fasta", required=True, help="Reference FASTA used for mapping (must have .fai).")
    ap.add_argument("--out-tsv", required=True, help="Output TSV (long format).")
    ap.add_argument("--params-dir", default=None, help="Directory containing *.sync.params. Default=sync-dir")
    ap.add_argument("--params-suffix", default=".params", help="Suffix for params files.")
    ap.add_argument("--sync-glob", default="*.sync", help="Glob pattern for sync files.")
    ap.add_argument("--min-cov", type=int, default=0, help="Minimum coverage (A+C+G+T) to include a site. 0=off")
    ap.add_argument("--genetic-code", type=int, default=1, choices=[1, 2],
                    help="Translation table: 1=standard nuclear, 2=vertebrate mitochondrial.")
    ap.add_argument("--gene-from-filename", action="store_true",
                    help="Use filename stem as GENE name (ignore GENE in params).")
    ap.add_argument("--skip-missing-params", action="store_true",
                    help="Skip sync if corresponding params not found.")
    args = ap.parse_args()

    pool_sizes = read_pool_sizes(args.pool_sizes)

    code_map = STD_CODE if args.genetic_code == 1 else MITO2_CODE
    fetcher = FastaFetcher(args.ref_fasta)

    sync_dir = args.sync_dir
    params_dir = args.params_dir or sync_dir

    sync_paths = sorted(glob.glob(os.path.join(sync_dir, args.sync_glob)))
    if not sync_paths:
        raise FileNotFoundError(f"No sync files found: {os.path.join(sync_dir, args.sync_glob)}")

    # Ensure output dir exists
    out_dir = os.path.dirname(os.path.abspath(args.out_tsv))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Output header
    with open(args.out_tsv, "w") as out:
        out.write("\t".join([
            "gene", "population",
            "piN", "piS", "piN_piS",
            "sites_used",
            "ref_fasta", "genetic_code"
        ]) + "\n")

        for sync_path in sync_paths:
            fn = os.path.basename(sync_path)
            stem = fn
            if stem.endswith(".sync"):
                stem = stem[:-5]

            params_path = os.path.join(params_dir, fn + args.params_suffix)  # aars2.sync + .sync.params
            if not os.path.exists(params_path):
                # also allow: <stem>.sync.params
                alt = os.path.join(params_dir, stem + args.params_suffix)
                if os.path.exists(alt):
                    params_path = alt

            if not os.path.exists(params_path):
                msg = f"[WARN] missing params for {fn}: {params_path}"
                if args.skip_missing_params:
                    print(msg + " -> skip", file=sys.stderr)
                    continue
                else:
                    raise FileNotFoundError(msg)

            gene_in_params, strand, blocks = parse_params_blocks(params_path)
            gene = stem if args.gene_from_filename else (gene_in_params or stem)

            # Build CDS maps
            pos2cdsidx, codon2gpos = build_cds_position_map(blocks, strand)

            # Initialize accumulators per population (in sync column order)
            # We'll read populations from the sync file columns.
            # But your sync files (per gene) are created consistently across 27 pops.
            pop_names = None
            acc_piN = []
            acc_piS = []
            acc_sites = []

            # Parse sync stream
            with open(sync_path, "r") as f:
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 4:
                        continue

                    chrom = parts[0]
                    try:
                        pos1 = int(parts[1])
                    except:
                        continue
                    ref_base = parts[2].upper()

                    # Initialize populations by counting sample columns
                    if pop_names is None:
                        nsamp = len(parts) - 3
                        # your population order is the bamlist order; we don't have names here.
                        # So we label them as S1..Sn unless user provides POPLIST elsewhere.
                        # BUT: your downstream expects real pop names -> you should pass a poplist
                        # by keeping sync column order consistent with pool_sizes.tsv.
                        #
                        # Here we assume sync column order == sorted(pool_sizes keys) is NOT safe.
                        # Therefore: require a sidecar file <sync>.pops if exists, else try infer from params header.
                        side = sync_path + ".pops"
                        if os.path.exists(side):
                            with open(side, "r") as pf:
                                pops = [x.strip() for x in pf if x.strip()]
                            if len(pops) != nsamp:
                                raise ValueError(f"{side} has {len(pops)} pops but sync has {nsamp} samples")
                            pop_names = pops
                        else:
                            # Fallback: if the sync was generated from your bamlist, order is your 27 pops
                            # (FG, LG, ...). We'll attempt to infer from pool_sizes.tsv header order by using
                            # the exact list you used earlier (common in your pipeline).
                            # If mismatch happens, you MUST create <sync>.pops.
                            default_order = [
                                "THE","JOE","BEA","MUC","PYE","AMO","SAY","GOS","ROB","FG",
                                "BOOT","ECHO","FRED","LAW","PACH","RS","SC","LB","CH","LG",
                                "SR","SL","TL","WB","WT","WK","SWA"
                            ]
                            if len(default_order) != nsamp:
                                raise ValueError(
                                    f"Sync has {nsamp} samples but no {side}.pops; "
                                    f"default list has {len(default_order)}. "
                                    f"Create {side} with one population per line (in sync column order)."
                                )
                            pop_names = default_order

                        acc_piN = [0.0] * len(pop_names)
                        acc_piS = [0.0] * len(pop_names)
                        acc_sites = [0] * len(pop_names)

                        # Validate pool sizes availability
                        missing = [p for p in pop_names if p not in pool_sizes]
                        if missing:
                            raise ValueError(f"Missing pool sizes for populations: {missing}")

                    gpos = (chrom, pos1)
                    if gpos not in pos2cdsidx:
                        continue

                    cds_idx = pos2cdsidx[gpos]
                    codon_index = cds_idx // 3
                    pos_in_codon = cds_idx % 3
                    if codon_index not in codon2gpos:
                        continue
                    codon_gpos = codon2gpos[codon_index]  # list of 3 (chrom,pos)

                    # Fetch reference codon in coding orientation
                    try:
                        b0 = fetcher.fetch_base_1based(codon_gpos[0][0], codon_gpos[0][1])
                        b1 = fetcher.fetch_base_1based(codon_gpos[1][0], codon_gpos[1][1])
                        b2 = fetcher.fetch_base_1based(codon_gpos[2][0], codon_gpos[2][1])
                    except Exception as e:
                        raise RuntimeError(f"FASTA fetch failed at {codon_gpos}: {e}")

                    ref_codon = (b0 + b1 + b2).upper()
                    if strand == "-":
                        ref_codon = revcomp(ref_codon)

                    # per-pop sample fields
                    sample_fields = parts[3:]
                    if len(sample_fields) != len(pop_names):
                        raise ValueError(f"{sync_path}: sample columns changed within file")

                    for i, field in enumerate(sample_fields):
                        counts = parse_sync_counts_field(field)
                        if counts is None:
                            continue
                        cov = counts["A"] + counts["C"] + counts["G"] + counts["T"]
                        if args.min_cov and cov < args.min_cov:
                            continue
                        if cov <= 0:
                            continue

                        freqs = {b: counts[b] / cov for b in BASES}

                        # raw syn/nonsyn pairwise diversity at this site
                        piN_raw, piS_raw = calc_site_pin_pis(freqs, ref_codon, pos_in_codon, code_map)

                        # finite sample correction using 2N chromosomes (diploid pools)
                        nchrom = 2 * int(pool_sizes[pop_names[i]])
                        if nchrom <= 1:
                            continue
                        corr = nchrom / (nchrom - 1.0)

                        acc_piN[i] += corr * piN_raw
                        acc_piS[i] += corr * piS_raw
                        acc_sites[i] += 1

            # Write gene results
            for i, pop in enumerate(pop_names):
                sites = acc_sites[i]
                if sites == 0:
                    piN = "NA"
                    piS = "NA"
                    ratio = "NA"
                else:
                    piN_val = acc_piN[i] / sites
                    piS_val = acc_piS[i] / sites
                    piN = f"{piN_val:.6g}"
                    piS = f"{piS_val:.6g}"
                    ratio = "NA" if piS_val == 0 else f"{(piN_val/piS_val):.6g}"

                out.write("\t".join([
                    gene, pop, str(piN), str(piS), str(ratio),
                    str(sites),
                    args.ref_fasta, str(args.genetic_code)
                ]) + "\n")

    print(f"[OK] wrote: {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
