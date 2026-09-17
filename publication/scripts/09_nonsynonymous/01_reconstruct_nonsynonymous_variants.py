#!/usr/bin/env python3

import csv
import gzip
import itertools
import math
import os
import random
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

try:
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "ERROR: Biopython is unavailable.\n"
        "Run: conda activate poolseq_env; conda install -c conda-forge biopython"
    )

# ============================================================
# Configuration
# ============================================================

SYNC_DIR = Path("/path/to/data/nu_287/sync")
BAMLIST = Path("/path/to/data/oxphos_gene_tree/bamlist_nuclear.txt")

GTF = Path("/path/to/workspace/OXPHOS_structural_validation/287_recalled_nonsyn/canonical_289_loci_from_phase.gtf")

OUTDIR = Path(
    "/path/to/workspace/OXPHOS_structural_validation/287_recalled_nonsyn"
)

MIN_DEPTH = 20
MIN_FRESH_POPS = 3
HIGH_EFFECT_DELTA = 0.50
PERM_THRESHOLD = 0.05

AK_FRESH = ["FG", "LG", "SR", "SL", "TL", "WB", "WT", "WK"]
BC_FRESH = [
    "SWA", "THE", "JOE", "BEA", "MUC", "PYE",
    "AMO", "BOOT", "ECHO", "LAW", "GOS", "ROB"
]

MARINE = {
    "AK": "RS",
    "BC": "SAY",
}

FRESH = {
    "AK": AK_FRESH,
    "BC": BC_FRESH,
}

KEEP_POPS = set(AK_FRESH + BC_FRESH + ["RS", "SAY"])

BASE_ORDER = ["A", "T", "C", "G"]
COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N",
}

random.seed(287)

OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Helpers
# ============================================================

def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_pop(path):
    name = os.path.basename(path.strip())

    name = re.sub(r"_subset\.bam$", "", name)
    name = re.sub(r"_mtDNA\.bam$", "", name)
    name = re.sub(r"\.bam$", "", name)

    match = re.match(r"^(?:\d+_)?([A-Za-z]+)(?:_S\d+)?$", name)
    if match:
        return match.group(1).upper()

    return name.upper()


def parse_attributes(text):
    attrs = {}

    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue

        match = re.match(r'(\S+)\s+"([^"]+)"', item)
        if match:
            attrs[match.group(1)] = match.group(2)
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()

    return attrs


def parse_sync_counts(cell):
    fields = cell.split(":")

    if len(fields) < 4:
        return None

    try:
        # Correct PoPoolation2 order: A:T:C:G:N:del
        A, T, C, G = map(int, fields[:4])
    except ValueError:
        return None

    return {
        "A": A,
        "T": T,
        "C": C,
        "G": G,
    }


def transcript_base(genomic_base, strand):
    genomic_base = genomic_base.upper()

    if strand == "+":
        return genomic_base

    return COMPLEMENT.get(genomic_base, "N")


def genomic_base(transcript_nt, strand):
    transcript_nt = transcript_nt.upper()

    if strand == "+":
        return transcript_nt

    return COMPLEMENT.get(transcript_nt, "N")


def translate_codon(codon):
    codon = codon.upper()

    if len(codon) != 3 or any(base not in "ATCG" for base in codon):
        return "X"

    return str(Seq(codon).translate())


def mean(values):
    return sum(values) / len(values) if values else math.nan


def sample_sd(values):
    if len(values) < 2:
        return math.nan
    return statistics.stdev(values)


def median(values):
    return statistics.median(values) if values else math.nan


def sign_flip_permutation(values):
    """
    Two-sided sign-flip permutation test for consistent ΔAF across
    freshwater populations. Exact for <=14 populations.
    """

    values = [float(x) for x in values if math.isfinite(float(x))]

    if len(values) < 2:
        return math.nan, 0

    observed = abs(mean(values))
    n = len(values)

    if n <= 14:
        permuted = []

        for signs in itertools.product([-1, 1], repeat=n):
            statistic = abs(
                mean([value * sign for value, sign in zip(values, signs)])
            )
            permuted.append(statistic)

        p = sum(x >= observed - 1e-15 for x in permuted) / len(permuted)
        return p, len(permuted)

    nperm = 9999
    extreme = 0

    for _ in range(nperm):
        statistic = abs(
            mean([
                value * random.choice([-1, 1])
                for value in values
            ])
        )

        if statistic >= observed - 1e-15:
            extreme += 1

    return (extreme + 1) / (nperm + 1), nperm


# ============================================================
# Read sample order
# ============================================================

if not BAMLIST.exists():
    sys.exit(f"ERROR: bamlist not found: {BAMLIST}")

bam_entries = [
    line.strip()
    for line in open(BAMLIST)
    if line.strip()
]

all_pops = [normalize_pop(x) for x in bam_entries]

if len(all_pops) != len(set(all_pops)):
    duplicates = sorted({
        pop for pop in all_pops
        if all_pops.count(pop) > 1
    })
    sys.exit(f"ERROR: duplicated normalized populations: {duplicates}")

missing_pops = sorted(KEEP_POPS - set(all_pops))

if missing_pops:
    sys.exit(
        "ERROR: required populations absent from bamlist: "
        + ",".join(missing_pops)
    )

pop_to_column = {
    pop: index
    for index, pop in enumerate(all_pops)
}

print("[Sample order]")
print("number of samples:", len(all_pops))
print(",".join(all_pops))
print("kept:", ",".join(sorted(KEEP_POPS)))

# ============================================================
# Read sync manifest
# ============================================================

sync_records = []

pattern = re.compile(
    r"^(.+)\.(ENSGACG\d+)\.sync(?:\.gz)?$",
    flags=re.IGNORECASE,
)

for sync_file in sorted(SYNC_DIR.glob("*.sync")):
    match = pattern.match(sync_file.name)

    if not match:
        print(
            f"WARNING: filename not recognized: {sync_file}",
            file=sys.stderr,
        )
        continue

    gene = match.group(1).lower()
    gene_id = match.group(2).upper()

    sync_records.append({
        "gene": gene,
        "gene_id": gene_id,
        "sync_file": sync_file,
    })

print("\n[Sync]")
print("recognized sync files:", len(sync_records))

manifest_file = OUTDIR / "all_sync_manifest.tsv"

with open(manifest_file, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=["gene", "gene_id", "sync_file"],
    )
    writer.writeheader()

    for record in sync_records:
        writer.writerow({
            "gene": record["gene"],
            "gene_id": record["gene_id"],
            "sync_file": str(record["sync_file"]),
        })

# ============================================================
# Read canonical CDS annotations
# ============================================================

if not GTF.exists():
    sys.exit(
        f"ERROR: GTF not found: {GTF}\n"
        "Edit GTF near the top of this script."
    )

cds_by_transcript = defaultdict(list)
transcript_meta = {}

with open_text(GTF) as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 9 or fields[2] != "CDS":
            continue

        chromosome = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        phase = fields[7]
        attrs = parse_attributes(fields[8])

        gene_id = attrs.get("gene_id", "").upper()
        transcript_id = attrs.get("transcript_id", "")
        gene_name = (
            attrs.get("gene_name")
            or attrs.get("gene")
            or attrs.get("Name")
            or gene_id
        ).lower()

        if not gene_id or not transcript_id:
            continue

        cds_by_transcript[transcript_id].append({
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "strand": strand,
            "phase": phase,
        })

        transcript_meta[transcript_id] = {
            "gene_id": gene_id,
            "gene": gene_name,
            "strand": strand,
            "chromosome": chromosome,
        }

# Select longest CDS transcript for each gene ID.
# The final GTF should already contain canonical transcripts, but this
# protects against accidental duplicate transcripts.

transcripts_by_gene = defaultdict(list)

for transcript_id, segments in cds_by_transcript.items():
    gene_id = transcript_meta[transcript_id]["gene_id"]
    cds_length = sum(
        segment["end"] - segment["start"] + 1
        for segment in segments
    )

    transcripts_by_gene[gene_id].append(
        (cds_length, transcript_id)
    )

canonical_transcript = {}

for gene_id, candidates in transcripts_by_gene.items():
    candidates.sort(reverse=True)
    canonical_transcript[gene_id] = candidates[0][1]

# ============================================================
# Construct transcript-oriented CDS coordinate maps
# ============================================================

cds_coordinate_map = {}
annotation_checks = []

for record in sync_records:
    gene = record["gene"]
    gene_id = record["gene_id"]

    transcript_id = canonical_transcript.get(gene_id)

    if transcript_id is None:
        annotation_checks.append({
            "gene": gene,
            "gene_id": gene_id,
            "transcript": "NA",
            "status": "NO_CDS_TRANSCRIPT",
            "cds_length": "NA",
            "cds_mod3": "NA",
            "strand": "NA",
        })
        continue

    segments = cds_by_transcript[transcript_id]
    strand = transcript_meta[transcript_id]["strand"]

    if strand == "+":
        segments = sorted(
            segments,
            key=lambda x: (x["start"], x["end"]),
        )
    else:
        segments = sorted(
            segments,
            key=lambda x: (x["start"], x["end"]),
            reverse=True,
        )

    coordinates = []

    for segment in segments:
        chromosome = segment["chromosome"]

        if strand == "+":
            positions = range(segment["start"], segment["end"] + 1)
        else:
            positions = range(segment["end"], segment["start"] - 1, -1)

        for position in positions:
            coordinates.append((chromosome, position))

    cds_length = len(coordinates)

    status = "PASS"

    if cds_length % 3 != 0:
        status = "CDS_LENGTH_NOT_MOD3"

    cds_coordinate_map[gene_id] = {
        "gene": gene,
        "gene_id": gene_id,
        "transcript": transcript_id,
        "strand": strand,
        "coordinates": coordinates,
    }

    annotation_checks.append({
        "gene": gene,
        "gene_id": gene_id,
        "transcript": transcript_id,
        "status": status,
        "cds_length": cds_length,
        "cds_mod3": cds_length % 3,
        "strand": strand,
    })

annotation_file = OUTDIR / "canonical_CDS_annotation_check.tsv"

with open(annotation_file, "w", newline="") as handle:
    fields = [
        "gene", "gene_id", "transcript", "status",
        "cds_length", "cds_mod3", "strand",
    ]
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(annotation_checks)

# ============================================================
# Main analysis
# ============================================================

variant_population_rows = []
variant_values = defaultdict(dict)
consensus_proteins = []
gene_qc = []

for number, record in enumerate(sync_records, start=1):
    gene = record["gene"]
    gene_id = record["gene_id"]
    sync_file = record["sync_file"]

    print(
        f"[{number}/{len(sync_records)}] "
        f"{gene} {gene_id}"
    )

    cds_info = cds_coordinate_map.get(gene_id)

    if cds_info is None:
        print("  SKIP: no CDS annotation")
        continue

    transcript_id = cds_info["transcript"]
    strand = cds_info["strand"]
    coordinates = cds_info["coordinates"]

    sync_data = {}

    with open_text(sync_file) as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split()

            if len(fields) < 3 + len(all_pops):
                continue

            chromosome = fields[0]
            position = int(fields[1])
            reference = fields[2].upper()
            cells = fields[3:3 + len(all_pops)]

            sync_data[(chromosome, position)] = {
                "reference": reference,
                "cells": cells,
            }

    reference_cds = []
    missing_reference_positions = 0

    for chromosome, position in coordinates:
        site = sync_data.get((chromosome, position))

        if site is None:
            reference_cds.append("N")
            missing_reference_positions += 1
            continue

        reference_cds.append(
            transcript_base(site["reference"], strand)
        )

    reference_cds = "".join(reference_cds)

    # Population consensus CDS/proteins
    for pop in sorted(KEEP_POPS):
        sample_index = pop_to_column[pop]
        population_cds = []

        for transcript_nt_index, (chromosome, position) in enumerate(coordinates):
            site = sync_data.get((chromosome, position))

            if site is None:
                population_cds.append("N")
                continue

            counts = parse_sync_counts(site["cells"][sample_index])

            if counts is None:
                population_cds.append("N")
                continue

            depth = sum(counts.values())

            if depth < MIN_DEPTH:
                population_cds.append("N")
                continue

            major_genomic_base = max(
                BASE_ORDER,
                key=lambda base: counts[base],
            )

            population_cds.append(
                transcript_base(major_genomic_base, strand)
            )

        population_cds = "".join(population_cds)

        protein_parts = []

        for codon_start in range(0, len(population_cds) - 2, 3):
            codon = population_cds[codon_start:codon_start + 3]
            protein_parts.append(translate_codon(codon))

        protein = "".join(protein_parts)

        consensus_proteins.append({
            "gene": gene,
            "gene_id": gene_id,
            "transcript": transcript_id,
            "population": pop,
            "protein": protein,
        })

    # Variant-specific amino-acid frequencies
    for cds_index, (chromosome, position) in enumerate(coordinates):
        site = sync_data.get((chromosome, position))

        if site is None:
            continue

        codon_number = cds_index // 3 + 1
        codon_offset = cds_index % 3
        codon_start = (codon_number - 1) * 3
        ref_codon = reference_cds[codon_start:codon_start + 3]

        if len(ref_codon) != 3 or "N" in ref_codon:
            continue

        ref_aa = translate_codon(ref_codon)

        if ref_aa == "X":
            continue

        genomic_ref = site["reference"].upper()

        pooled_counts = {base: 0 for base in BASE_ORDER}

        for pop in KEEP_POPS:
            sample_index = pop_to_column[pop]
            counts = parse_sync_counts(site["cells"][sample_index])

            if counts is None:
                continue

            for base in BASE_ORDER:
                pooled_counts[base] += counts[base]

        for genomic_alt in BASE_ORDER:
            if genomic_alt == genomic_ref:
                continue

            if pooled_counts[genomic_alt] == 0:
                continue

            transcript_alt = transcript_base(genomic_alt, strand)

            alt_codon_list = list(ref_codon)
            alt_codon_list[codon_offset] = transcript_alt
            alt_codon = "".join(alt_codon_list)
            alt_aa = translate_codon(alt_codon)

            if alt_aa in {"X", ref_aa}:
                continue

            variant_id = (
                f"{gene}:{chromosome}:{position}:"
                f"{genomic_ref}>{genomic_alt}:"
                f"{ref_aa}{codon_number}{alt_aa}"
            )

            variant_key = (
                gene,
                gene_id,
                transcript_id,
                chromosome,
                position,
                genomic_ref,
                genomic_alt,
                codon_number,
                ref_codon,
                alt_codon,
                ref_aa,
                alt_aa,
                strand,
                variant_id,
            )

            for pop in KEEP_POPS:
                sample_index = pop_to_column[pop]
                counts = parse_sync_counts(site["cells"][sample_index])

                if counts is None:
                    continue

                depth = sum(counts.values())

                if depth < MIN_DEPTH:
                    continue

                alt_count = counts[genomic_alt]
                alt_frequency = alt_count / depth

                variant_values[variant_key][pop] = {
                    "alt_count": alt_count,
                    "depth": depth,
                    "alt_aa_frequency": alt_frequency,
                }

    gene_qc.append({
        "gene": gene,
        "gene_id": gene_id,
        "transcript": transcript_id,
        "strand": strand,
        "cds_length": len(coordinates),
        "sync_positions": len(sync_data),
        "missing_CDS_positions_in_sync": missing_reference_positions,
        "reference_X_codons": sum(
            translate_codon(reference_cds[i:i+3]) == "X"
            for i in range(0, len(reference_cds) - 2, 3)
        ),
    })

# ============================================================
# Regional marine-freshwater comparisons
# ============================================================

region_summary_rows = []

for variant_key, population_data in variant_values.items():
    (
        gene,
        gene_id,
        transcript_id,
        chromosome,
        position,
        genomic_ref,
        genomic_alt,
        codon_number,
        ref_codon,
        alt_codon,
        ref_aa,
        alt_aa,
        strand,
        variant_id,
    ) = variant_key

    for region in ["AK", "BC"]:
        marine_pop = MARINE[region]
        marine_data = population_data.get(marine_pop)

        if marine_data is None:
            continue

        marine_af = marine_data["alt_aa_frequency"]
        deltas = []
        freshwater_afs = []
        populations = []

        for pop in FRESH[region]:
            fresh_data = population_data.get(pop)

            if fresh_data is None:
                continue

            fresh_af = fresh_data["alt_aa_frequency"]
            delta = fresh_af - marine_af

            deltas.append(delta)
            freshwater_afs.append(fresh_af)
            populations.append(pop)

            variant_population_rows.append({
                "variant_id": variant_id,
                "gene": gene,
                "gene_id": gene_id,
                "transcript": transcript_id,
                "region": region,
                "marine_population": marine_pop,
                "freshwater_population": pop,
                "chromosome": chromosome,
                "position": position,
                "strand": strand,
                "genomic_ref": genomic_ref,
                "genomic_alt": genomic_alt,
                "codon_position": codon_number,
                "ref_codon": ref_codon,
                "alt_codon": alt_codon,
                "ref_AA": ref_aa,
                "alt_AA": alt_aa,
                "marine_altAA_frequency": marine_af,
                "freshwater_altAA_frequency": fresh_af,
                "deltaAA_frequency": delta,
                "abs_deltaAA_frequency": abs(delta),
                "marine_depth": marine_data["depth"],
                "freshwater_depth": fresh_data["depth"],
            })

        if len(deltas) < MIN_FRESH_POPS:
            continue

        p_perm, n_perm = sign_flip_permutation(deltas)

        n_positive = sum(delta > 0 for delta in deltas)
        n_negative = sum(delta < 0 for delta in deltas)

        median_delta = median(deltas)
        direction = (
            "positive" if median_delta > 0
            else "negative" if median_delta < 0
            else "zero"
        )

        region_summary_rows.append({
            "variant_id": variant_id,
            "gene": gene,
            "gene_id": gene_id,
            "transcript": transcript_id,
            "region": region,
            "marine_population": marine_pop,
            "chromosome": chromosome,
            "position": position,
            "strand": strand,
            "genomic_ref": genomic_ref,
            "genomic_alt": genomic_alt,
            "codon_position": codon_number,
            "ref_codon": ref_codon,
            "alt_codon": alt_codon,
            "ref_AA": ref_aa,
            "alt_AA": alt_aa,
            "n_freshwater_populations": len(populations),
            "freshwater_populations": ",".join(populations),
            "marine_altAA_frequency": marine_af,
            "mean_freshwater_altAA_frequency": mean(freshwater_afs),
            "median_freshwater_altAA_frequency": median(freshwater_afs),
            "mean_deltaAA_frequency": mean(deltas),
            "median_deltaAA_frequency": median_delta,
            "sd_deltaAA_frequency": sample_sd(deltas),
            "max_abs_deltaAA_frequency": max(abs(x) for x in deltas),
            "n_abs_delta_ge_0.5": sum(
                abs(x) >= HIGH_EFFECT_DELTA for x in deltas
            ),
            "n_positive": n_positive,
            "n_negative": n_negative,
            "direction": direction,
            "p_permutation": p_perm,
            "n_permutations": n_perm,
        })

# ============================================================
# Cross-region overlap
# ============================================================

summary_by_aa = defaultdict(dict)

for row in region_summary_rows:
    key = (
        row["gene"],
        row["codon_position"],
        row["ref_AA"],
        row["alt_AA"],
    )

    summary_by_aa[key][row["region"]] = row

overlap_rows = []
candidate_rows = []

for key, regions in summary_by_aa.items():
    if "AK" not in regions or "BC" not in regions:
        continue

    ak = regions["AK"]
    bc = regions["BC"]

    ak_delta = float(ak["median_deltaAA_frequency"])
    bc_delta = float(bc["median_deltaAA_frequency"])

    concordant = (
        ak_delta != 0
        and bc_delta != 0
        and (ak_delta > 0) == (bc_delta > 0)
    )

    high_effect_both = (
        abs(ak_delta) >= HIGH_EFFECT_DELTA
        and abs(bc_delta) >= HIGH_EFFECT_DELTA
    )

    perm_both = (
        float(ak["p_permutation"]) < PERM_THRESHOLD
        and float(bc["p_permutation"]) < PERM_THRESHOLD
    )

    output = {
        "gene": key[0],
        "codon_position": key[1],
        "ref_AA": key[2],
        "alt_AA": key[3],
        "AK_variant_id": ak["variant_id"],
        "BC_variant_id": bc["variant_id"],
        "AK_marine_AF": ak["marine_altAA_frequency"],
        "BC_marine_AF": bc["marine_altAA_frequency"],
        "AK_median_freshwater_AF":
            ak["median_freshwater_altAA_frequency"],
        "BC_median_freshwater_AF":
            bc["median_freshwater_altAA_frequency"],
        "AK_median_delta":
            ak["median_deltaAA_frequency"],
        "BC_median_delta":
            bc["median_deltaAA_frequency"],
        "AK_p_permutation":
            ak["p_permutation"],
        "BC_p_permutation":
            bc["p_permutation"],
        "direction_concordant": concordant,
        "high_effect_both": high_effect_both,
        "permutation_P_below_0.05_both": perm_both,
    }

    overlap_rows.append(output)

    if concordant and high_effect_both and perm_both:
        candidate_rows.append(output)

# ============================================================
# Write outputs
# ============================================================

def write_table(path, rows):
    if not rows:
        Path(path).write_text("")
        return

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


write_table(
    OUTDIR / "gene_processing_QC.tsv",
    gene_qc,
)

write_table(
    OUTDIR / "nonsynonymous_variant_population_frequencies.tsv",
    variant_population_rows,
)

write_table(
    OUTDIR / "nonsynonymous_variant_region_summary.tsv",
    region_summary_rows,
)

write_table(
    OUTDIR / "AK_BC_nonsynonymous_overlap.tsv",
    overlap_rows,
)

write_table(
    OUTDIR / "AK_BC_high_effect_permutation_candidates.tsv",
    candidate_rows,
)

# Consensus proteins
consensus_fasta = OUTDIR / "population_consensus_proteins_287.fasta"

with open(consensus_fasta, "w") as handle:
    for row in consensus_proteins:
        header = (
            f">{row['gene']}|{row['gene_id']}|"
            f"{row['transcript']}|{row['population']}"
        )

        handle.write(header + "\n")
        protein = row["protein"]

        for start in range(0, len(protein), 60):
            handle.write(protein[start:start + 60] + "\n")

print("\n[Written]")
for path in [
    manifest_file,
    annotation_file,
    OUTDIR / "gene_processing_QC.tsv",
    OUTDIR / "nonsynonymous_variant_population_frequencies.tsv",
    OUTDIR / "nonsynonymous_variant_region_summary.tsv",
    OUTDIR / "AK_BC_nonsynonymous_overlap.tsv",
    OUTDIR / "AK_BC_high_effect_permutation_candidates.tsv",
    consensus_fasta,
]:
    print(path)

print("\n[Counts]")
print("sync genes:", len(sync_records))
print("regional nonsynonymous summaries:", len(region_summary_rows))
print("AK–BC overlaps:", len(overlap_rows))
print("strict high-effect candidates:", len(candidate_rows))
