#piNpiS for 287+13 genes


#00) scripts/00_make_pool_sizes.sh
#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/spareHD_2/nu_287/pinpis"
IN="${ROOT}/meta/pool_sizes.raw.tsv"
OUT="${ROOT}/meta/pool_sizes.tsv"

mkdir -p "${ROOT}/meta"

awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1 {next}                     # skip header (pool_ID Population Region n_individuals)
  NF<4 {next}
  $2=="" || $4=="" {next}
  {
    pop=$2; n=$4;
    gsub(/^[ \t]+|[ \t]+$/,"",pop);
    gsub(/^[ \t]+|[ \t]+$/,"",n);
    if(pop=="Population") next;
    if(n ~ /^[0-9]+$/) print pop, n;
  }' "$IN" \
| sort -k1,1 \
| awk 'BEGIN{OFS="\t"} NR==1{print "Population","n_individuals"} {print}' \
> "$OUT"

echo "[OK] wrote: $OUT"
echo "[OK] lines: $(wc -l < "$OUT")"

01) scripts/01_make_pop_order_nuclear.sh

从 pool_sizes.tsv 生成固定 population 顺序（后面 merge/画图用）

#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/spareHD_2/nu_287/pinpis"
POOL="${ROOT}/meta/pool_sizes.tsv"
OUT="${ROOT}/meta/pop_order.tsv"

awk -F'\t' 'NR==1{next} NF>=1{print $1}' "$POOL" > "$OUT"
echo "[OK] wrote: $OUT ($(wc -l < "$OUT") pops)"

# 01) scripts/01_make_nuclear_params_from_phase.py
#!/usr/bin/env python3
import argparse, os, glob, re
from collections import defaultdict

def read_phase(phase_tsv):
    """
    phase.tsv columns (expected):
      chr start end strand phase tx gene
    """
    blocks = defaultdict(lambda: {"strand": None, "rows": []})
    with open(phase_tsv) as f:
        for line in f:
            line=line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts=line.split("\t")
            if len(parts) < 7:
                continue
            chrom, start, end, strand, phase, tx, gene = parts[:7]
            gene = gene.strip()
            if gene == "":
                continue
            d = blocks[gene]
            if d["strand"] is None:
                d["strand"] = strand
            d["rows"].append((chrom, int(start), int(end), int(phase)))
    # sort blocks
    for g in blocks:
        blocks[g]["rows"].sort(key=lambda x:(x[0], x[1], x[2]))
    return blocks

def sync_gene_from_filename(sync_path):
    # aars2.ENSG...sync  -> gene=aars2
    base=os.path.basename(sync_path)
    return base.split(".")[0]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--phase-tsv", required=True)
    ap.add_argument("--sync-dir", required=True)
    ap.add_argument("--out-dir", default=None, help="Default=sync-dir")
    ap.add_argument("--suffix", default=".sync.params", help="Default=.sync.params")
    ap.add_argument("--skip-genes", default="", help="Comma-separated gene symbols to skip (case-insensitive)")
    args=ap.parse_args()

    out_dir = args.out_dir or args.sync_dir
    os.makedirs(out_dir, exist_ok=True)

    skip=set([x.strip().lower() for x in args.skip_genes.split(",") if x.strip()])

    phase_blocks = read_phase(args.phase_tsv)

    wrote=0
    skipped=0
    for sp in sorted(glob.glob(os.path.join(args.sync_dir, "*.sync"))):
        gene = sync_gene_from_filename(sp)
        gl = gene.lower()
        if gl in skip:
            skipped += 1
            continue
        if gene not in phase_blocks:
            skipped += 1
            continue

        strand = phase_blocks[gene]["strand"] or "+"
        rows = phase_blocks[gene]["rows"]
        if not rows:
            skipped += 1
            continue

        outp = os.path.join(out_dir, os.path.basename(sp) + args.suffix)
        with open(outp, "w") as o:
            o.write("# pinpis params (auto-generated)\n")
            o.write(f"GENE\t{gene}\n")
            o.write(f"STRAND\t{strand}\n")
            o.write("CDS_BLOCKS\tchr\tstart\tend\tphase\n")
            for chrom, s, e, ph in rows:
                o.write(f"BLOCK\t{chrom}\t{s}\t{e}\t{ph}\n")
        wrote += 1

    print(f"[OK] wrote params for {wrote} genes -> {out_dir}")
    if skipped:
        print(f"[WARN] skipped {skipped} sync files (no match / empty blocks / in skip list)")

if __name__ == "__main__":
    main()


#01) scripts/01_make_mt_params_from_gff.py

#!/usr/bin/env python3
import argparse, os, re

PCG = set(["ND1","ND2","ND3","ND4","ND4L","ND5","ND6","COX1","COX2","COX3","CYTB","ATP6","ATP8"])

def parse_attrs(attr):
    d={}
    for kv in attr.split(";"):
        if "=" in kv:
            k,v=kv.split("=",1)
            d[k]=v
    return d

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="mt genes.gff (GFF3)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--contig", default=None, help="If set, replace contig name to this (e.g. MH205729.1)")
    args=ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # collect CDS features for PCG
    cds_by_gene = {}
    strand_by_gene = {}

    with open(args.gff) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts=line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, src, feat, start, end, score, strand, phase, attr = parts
            if feat != "CDS":
                continue
            attrs=parse_attrs(attr)
            gene = attrs.get("gene", attrs.get("Name",""))
            if gene not in PCG:
                continue
            chrom2 = args.contig if args.contig else chrom
            cds_by_gene.setdefault(gene, []).append((chrom2, int(start), int(end), 0))
            strand_by_gene[gene] = strand

    wrote=0
    for gene, rows in sorted(cds_by_gene.items()):
        rows.sort(key=lambda x:(x[0], x[1], x[2]))
        outp=os.path.join(args.out_dir, f"{gene}.sync.params")
        with open(outp,"w") as o:
            o.write("# pinpis params (auto-generated)\n")
            o.write(f"GENE\t{gene}\n")
            o.write(f"STRAND\t{strand_by_gene.get(gene,'+')}\n")
            o.write("CDS_BLOCKS\tchr\tstart\tend\tphase\n")
            for chrom,s,e,ph in rows:
                o.write(f"BLOCK\t{chrom}\t{s}\t{e}\t{ph}\n")
        wrote += 1

    print(f"[OK] wrote mt params: {wrote} genes -> {args.out_dir}")

if __name__=="__main__":
    main()

# 02) scripts/02_pinpis_from_sync.py
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



#03) scripts/03_run_pinpis_nuclear_all.sh （忽略 si gene）


#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/spareHD_2/nu_287/pinpis"

SYNC_DIR="/mnt/spareHD_2/nu_287/sync"
POOL="${ROOT}/meta/pool_sizes.tsv"
REF="/work/cyu/stickleback_nuclear_only.fa"

OUTDIR="${ROOT}/results"
OUT="${OUTDIR}/pinpis_nuclear.tsv"

mkdir -p "$OUTDIR" "${ROOT}/logs"

# ---- ignore this gene (file) ----
# file name contains: si_dkey-31b16.7.ENSG....sync
IGNORE_REGEX='^si_dkey-31b16\.7\.'

python "${ROOT}/scripts/02_pinpis_from_sync.py" \
  --sync-dir   "$SYNC_DIR" \
  --sync-glob  "*.sync" \
  --pool-sizes "$POOL" \
  --ref-fasta  "$REF" \
  --genetic-code 1 \
  --params-dir "$SYNC_DIR" \
  --params-suffix ".params" \
  --skip-missing-params \
  --out-tsv    "$OUT" \
  > "${ROOT}/logs/pinpis_nuclear.log" 2>&1

# 过滤掉 si 那个 gene（双保险：有些会写 gene=si:dkey-31b16.7）
tmp="${OUT}.tmp"
awk -F'\t' -v re="$IGNORE_REGEX" '
  NR==1{print; next}
  $1 ~ /^si:dkey-31b16\.7$/ {next}
  $0 ~ re {next}
  {print}
' "$OUT" > "$tmp" && mv "$tmp" "$OUT"

echo "[OK] nuclear πNπS -> $OUT"


解释：你现在 nuclear 的 .params 已经写在 sync dir 里了，所以 params-dir 就指向 sync dir。
--skip-missing-params 保证缺一个也不会整个停掉。

#03) scripts/03_run_pinpis_mt13.sh


#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/spareHD_2/nu_287/pinpis"

SYNC_DIR="/work/cyu/poolseq/PPalign_output/ann_mt_pergene_sync"
POOL="${ROOT}/meta/pool_sizes.tsv"
REF="/work/cyu/sequence.fasta"
PARAMS_DIR="${ROOT}/params_mt"

OUTDIR="${ROOT}/results"
OUT="${OUTDIR}/pinpis_mt13.tsv"

mkdir -p "$OUTDIR" "${ROOT}/logs"

python "${ROOT}/scripts/02_pinpis_from_sync.py" \
  --sync-dir   "$SYNC_DIR" \
  --pool-sizes "$POOL" \
  --ref-fasta  "$REF" \
  --genetic-code 2 \
  --params-dir "$PARAMS_DIR" \
  --params-suffix ".params" \
  --skip-missing-params \
  --out-tsv    "$OUT" \
  > "${ROOT}/logs/pinpis_mt13.log" 2>&1

echo "[OK] mt13 πNπS -> $OUT"


⚠️ 注意：这要求你的 mt params 文件名是 ND1.sync.params 这种。
所以你要用我上面那个 01_make_mt_params_from_gff.py 生成后，再把它们改名成 .params 后缀匹配，最省事的做法是：

生成时就输出 .sync.params，再软链一个 .params：

mkdir -p /mnt/spareHD_2/nu_287/pinpis/params_mt
python /mnt/spareHD_2/nu_287/pinpis/scripts/01_make_mt_params_from_gff.py \
  --gff /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff \
  --out-dir /mnt/spareHD_2/nu_287/pinpis/params_mt \
  --contig MH205729.1

# 给每个 *.sync.params 建一个 *.sync.params? 现在脚本输出 ND1.sync.params
# 我们要让 pinpis 用 ".params"：ND1.sync + ".params" => ND1.sync.params （已经匹配）
# 所以这里不用改名，直接 OK。
ls -lh /mnt/spareHD_2/nu_287/pinpis/params_mt | head

#04) scripts/04_merge_pinpis_tables.py

把 nuclear + mt13 合并成一个总表（long format）

保存为：scripts/04_merge_pinpis_tables.py

#!/usr/bin/env python3
import argparse, pandas as pd

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--nuclear", required=True)
    ap.add_argument("--mt", required=True)
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    nu=pd.read_csv(args.nuclear, sep="\t")
    mt=pd.read_csv(args.mt, sep="\t")

    all_df=pd.concat([nu, mt], ignore_index=True)

    # normalize columns
    # expected columns: gene population piN piS piN_piS sites_used ref_fasta genetic_code
    # ensure consistent dtypes
    for c in ["piN","piS","piN_piS"]:
        if c in all_df.columns:
            all_df[c]=pd.to_numeric(all_df[c], errors="coerce")
    if "sites_used" in all_df.columns:
        all_df["sites_used"]=pd.to_numeric(all_df["sites_used"], errors="coerce").astype("Int64")

    # sort
    all_df = all_df.sort_values(["gene","population"])

    all_df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote: {args.out} rows={len(all_df)}")

if __name__=="__main__":
    main()

#05) scripts/05_run_all.sh （一键跑完）


#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/spareHD_2/nu_287/pinpis"
mkdir -p "${ROOT}/meta" "${ROOT}/results" "${ROOT}/logs" "${ROOT}/params_mt"

echo "[00] pool_sizes.tsv"
bash "${ROOT}/scripts/00_make_pool_sizes.sh"
bash "${ROOT}/scripts/01_make_pop_order_nuclear.sh"

echo "[01] nuclear params from phase.tsv -> sync-dir/*.params"
python "${ROOT}/scripts/01_make_nuclear_params_from_phase.py" \
  --phase-tsv /work/cyu/oxphos_from_ref_no_biomart/06_igv/dnds_annotations_final/oxphos_assembly.CDS.final.phase.tsv \
  --sync-dir  /mnt/spareHD_2/nu_287/sync \
  --out-dir   /mnt/spareHD_2/nu_287/sync \
  --suffix    ".params" \
  --skip-genes "si:dkey-31b16.7,si_dkey-31b16.7"

echo "[01] mt params from genes.gff -> params_mt/"
python "${ROOT}/scripts/01_make_mt_params_from_gff.py" \
  --gff /home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff \
  --out-dir "${ROOT}/params_mt" \
  --contig MH205729.1

echo "[03] run nuclear pinpis (ignore si gene)"
bash "${ROOT}/scripts/03_run_pinpis_nuclear_all.sh"

echo "[03] run mt13 pinpis"
bash "${ROOT}/scripts/03_run_pinpis_mt13.sh"

echo "[04] merge tables"
python "${ROOT}/scripts/04_merge_pinpis_tables.py" \
  --nuclear "${ROOT}/results/pinpis_nuclear.tsv" \
  --mt      "${ROOT}/results/pinpis_mt13.tsv" \
  --out     "${ROOT}/results/pinpis_all_287plus13.tsv"

echo "[DONE] ${ROOT}/results/pinpis_all_287plus13.tsv"



