#!/usr/bin/env python3

import os
import re
import pandas as pd
from collections import defaultdict

# ============================================================
# Paths
# ============================================================
TABLES3   = "/work/cyu/TableS3.csv"
ALIGN_DIR = "/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
OUTDIR    = "/work/cyu/provean_inputs_OXPHOS72"

os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# 1. Read all OXPHOS subunit genes from TableS3
# ============================================================
tab = pd.read_csv(TABLES3)

CANDIDATE_GENES = (
    tab[
        (tab["model"] == "M0") &
        (tab["role"] == "subunit")
    ]["gene"]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
    .tolist()
)

print("Number of subunit genes from TableS3:", len(CANDIDATE_GENES))

# ============================================================
# 2. Population sets
# ============================================================
AK_FW = {"FG","LG","SR","SL","TL","WB","WT","WK","LB"}
BC_FW = {"SWA","THE","JOE","BEA","MUC","PYE","BOOT","ECHO","LAW","GOS","ROB"}

MARINE_REF = {
    "AK": "RS",
    "BC": "SAY"
}

ALL_POPS = AK_FW | BC_FW | {"RS","SAY","AMO","PACH","FRED","SC","CH"}

BAD_AA = {"-", "X", "N", "*", "?"}

# ============================================================
# Helper functions
# ============================================================
def read_fasta(path):
    seqs = {}
    name = None
    chunks = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        seqs[name] = "".join(chunks)

    return seqs


def write_fasta(path, name, seq):
    with open(path, "w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")


def infer_pop(sample_name):
    parts = re.split(r"[_\.\-]+", sample_name)

    for p in parts:
        if p in ALL_POPS:
            return p

    for p in sorted(ALL_POPS, key=len, reverse=True):
        if re.search(rf"(^|[_\.\-]){p}([_\.\-]|$)", sample_name):
            return p

    return None


def aligned_to_ungapped_position(ref_aln):
    pos_map = {}
    pos = 0

    for i, aa in enumerate(ref_aln):
        if aa != "-":
            pos += 1
            pos_map[i] = pos
        else:
            pos_map[i] = None

    return pos_map


def ungap(seq):
    return seq.replace("-", "")


def get_variants(ref_aln, query_aln):
    pos_map = aligned_to_ungapped_position(ref_aln)
    variants = []

    for i, (r, q) in enumerate(zip(ref_aln, query_aln)):
        if pos_map[i] is None:
            continue
        if r in BAD_AA or q in BAD_AA:
            continue
        if r == q:
            continue

        pos = pos_map[i]
        variants.append(f"{r}{pos}{q}")

    return variants

# ============================================================
# 3. Main
# ============================================================
summary_rows = []
missing_genes = []

for gene in CANDIDATE_GENES:
    aln_file = os.path.join(ALIGN_DIR, gene, f"{gene}.pep.aln.faa")

    if not os.path.exists(aln_file):
        missing_genes.append(gene)
        print(f"[skip] missing alignment: {gene}")
        continue

    seqs = read_fasta(aln_file)

    pop_to_sample = {}
    for sample in seqs:
        pop = infer_pop(sample)
        if pop:
            pop_to_sample[pop] = sample

    for region, marine_pop in MARINE_REF.items():
        if marine_pop not in pop_to_sample:
            print(f"[skip] {gene} {region}: missing marine ref {marine_pop}")
            continue

        fw_pops = AK_FW if region == "AK" else BC_FW

        ref_sample = pop_to_sample[marine_pop]
        ref_aln = seqs[ref_sample]
        ref_protein = ungap(ref_aln)

        out_region = os.path.join(OUTDIR, f"{gene}_{region}")
        os.makedirs(out_region, exist_ok=True)

        ref_fa = os.path.join(out_region, f"{gene}_{region}_ref_{marine_pop}.faa")
        write_fasta(ref_fa, f"{gene}_{marine_pop}", ref_protein)

        variant_to_pops = defaultdict(list)

        for pop in sorted(fw_pops):
            if pop not in pop_to_sample:
                continue

            sample = pop_to_sample[pop]
            q_aln = seqs[sample]

            variants = get_variants(ref_aln, q_aln)

            for v in variants:
                variant_to_pops[v].append(pop)

            summary_rows.append({
                "gene": gene,
                "region": region,
                "marine_ref": marine_pop,
                "pop": pop,
                "sample": sample,
                "n_variants": len(variants),
                "variants": ",".join(variants)
            })

        var_file = os.path.join(out_region, f"{gene}_{region}_variants.txt")
        with open(var_file, "w") as out:
            for v in sorted(variant_to_pops):
                out.write(v + "\n")

        map_file = os.path.join(out_region, f"{gene}_{region}_variant_population_map.tsv")
        with open(map_file, "w") as out:
            out.write("gene\tregion\tmarine_ref\tvariant\tpopulations\tn_pops\n")
            for v in sorted(variant_to_pops):
                pops = sorted(set(variant_to_pops[v]))
                out.write(
                    f"{gene}\t{region}\t{marine_pop}\t{v}\t{','.join(pops)}\t{len(pops)}\n"
                )

        print(f"[ok] {gene} {region}")

# ============================================================
# 4. Write summary
# ============================================================
summary_file = os.path.join(OUTDIR, "OXPHOS72_AA_substitution_summary.tsv")

with open(summary_file, "w") as out:
    out.write("gene\tregion\tmarine_ref\tpop\tsample\tn_variants\tvariants\n")
    for r in summary_rows:
        out.write(
            f"{r['gene']}\t{r['region']}\t{r['marine_ref']}\t{r['pop']}\t"
            f"{r['sample']}\t{r['n_variants']}\t{r['variants']}\n"
        )

missing_file = os.path.join(OUTDIR, "missing_alignment_genes.txt")
with open(missing_file, "w") as out:
    for g in missing_genes:
        out.write(g + "\n")

print("\nDONE")
print("Summary:", summary_file)
print("Missing alignments:", missing_file)
print("Output:", OUTDIR)
