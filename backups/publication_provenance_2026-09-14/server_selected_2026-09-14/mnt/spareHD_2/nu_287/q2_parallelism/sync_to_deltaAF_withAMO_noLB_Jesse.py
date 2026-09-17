#!/usr/bin/env python3

import os
import re
import gzip
import sys


# ============================================================
# Inputs and outputs
# ============================================================

SYNC_DIR = "/mnt/spareHD_2/nu_287/sync"

BAMLIST = (
    "/mnt/spareHD_2/oxphos_gene_tree/"
    "bamlist_nuclear.txt"
)

OUT_AF = (
    "/mnt/spareHD_2/nu_287/q2_parallelism/"
    "af_long.withAMO_noLB_Jesse.tsv.gz"
)

OUT_META = (
    "/mnt/spareHD_2/nu_287/q2_parallelism/"
    "snp_meta.withAMO_noLB_Jesse.tsv.gz"
)


# ============================================================
# Parameters
# ============================================================

DROP = {"LB"}

# Final minimum depth used throughout the pipeline
MIN_DEPTH = 20


# ============================================================
# Helper functions
# ============================================================

def open_any(path):
    """
    Open a plain-text or gzipped file for reading.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def norm_pop(value):
    """
    Convert a BAM path or BAM filename into the population code.

    Examples:
      1_FG_S1_subset.bam -> FG
      FG_subset.bam      -> FG
      25_RS.bam          -> RS
    """

    value = os.path.basename(value)

    value = re.sub(
        r"_subset\.bam$",
        "",
        value
    )

    value = re.sub(
        r"\.bam$",
        "",
        value
    )

    value = re.sub(
        r"^(?:\d+_)?([A-Za-z]+)(?:_S\d+)?$",
        r"\1",
        value
    )

    return value.upper()


def gene_from_fn(filename):
    """
    Derive the gene name from a per-gene sync filename.

    Example:
      ndufa1.sync.gz -> ndufa1
    """

    filename = re.sub(
        r"\.sync(\.gz)?$",
        "",
        filename
    )

    return filename.split(".")[0].lower()


def parse_sync_cell(cell):
    """
    Parse one PoPoolation2 sync cell.

    PoPoolation2 sync order:
      A:T:C:G:N:deletion
    """

    values = cell.split(":")

    if len(values) < 4:
        raise ValueError(
            f"Malformed sync cell: {cell}"
        )

    A, T, C, G = map(int, values[:4])

    return A, T, C, G


# ============================================================
# Read population order from BAM list
# ============================================================

with open(BAMLIST, "r") as handle:
    pops_raw = [
        line.strip()
        for line in handle
        if line.strip()
    ]

pops = [
    norm_pop(value)
    for value in pops_raw
]

keep_idx = [
    index
    for index, pop in enumerate(pops)
    if pop not in DROP
]

keep_pops = [
    pops[index]
    for index in keep_idx
]


print(
    "[info] all populations:",
    ",".join(pops),
    file=sys.stderr
)

print(
    "[info] kept populations:",
    ",".join(keep_pops),
    file=sys.stderr
)

print(
    "[info] dropped populations:",
    ",".join(sorted(DROP)),
    file=sys.stderr
)

print(
    f"[info] minimum per-population depth: {MIN_DEPTH}",
    file=sys.stderr
)


# ============================================================
# Process all per-gene sync files
# ============================================================

with gzip.open(OUT_AF, "wt") as out_af, \
     gzip.open(OUT_META, "wt") as out_meta:

    out_af.write(
        "chr\tpos\tgene\tpop\taf\tdepth\tfocal_allele\n"
    )

    out_meta.write(
        "chr\tpos\tgene\tfocal_allele\n"
    )

    sync_files = sorted(
        filename
        for filename in os.listdir(SYNC_DIR)
        if filename.endswith(".sync")
        or filename.endswith(".sync.gz")
    )

    if len(sync_files) == 0:
        raise RuntimeError(
            f"No .sync or .sync.gz files found in {SYNC_DIR}"
        )

    for file_index, filename in enumerate(
        sync_files,
        start=1
    ):

        gene = gene_from_fn(filename)

        path = os.path.join(
            SYNC_DIR,
            filename
        )

        print(
            f"[progress] {file_index}/{len(sync_files)} "
            f"{filename}",
            file=sys.stderr
        )

        with open_any(path) as handle:

            for line_number, line in enumerate(
                handle,
                start=1
            ):

                if not line.strip():
                    continue

                fields = line.rstrip("\n").split()

                expected_fields = 3 + len(pops)

                if len(fields) < expected_fields:

                    print(
                        "[warning] skipped short row: "
                        f"{filename}:{line_number}; "
                        f"found {len(fields)} columns, "
                        f"expected at least {expected_fields}",
                        file=sys.stderr
                    )

                    continue

                chrom = fields[0]
                pos = fields[1]

                cells = fields[
                    3:(3 + len(pops))
                ]

                total_counts = {
                    "A": 0,
                    "T": 0,
                    "C": 0,
                    "G": 0
                }

                per_pop = []

                for index in keep_idx:

                    pop = pops[index]

                    A, T, C, G = parse_sync_cell(
                        cells[index]
                    )

                    depth = A + T + C + G

                    total_counts["A"] += A
                    total_counts["T"] += T
                    total_counts["C"] += C
                    total_counts["G"] += G

                    per_pop.append(
                        (
                            pop,
                            A,
                            T,
                            C,
                            G,
                            depth
                        )
                    )

                # ------------------------------------------------
                # Define one globally consistent focal allele
                # ------------------------------------------------
                #
                # The focal allele is the allele with the largest
                # total read count across all retained populations.
                #
                # This ensures that AF and deltaAF always refer to
                # the same allele in every population.
                #
                # This is a pooled-read-count definition, not an
                # unweighted mean AF across populations.
                # ------------------------------------------------

                focal = max(
                    total_counts.items(),
                    key=lambda item: item[1]
                )[0]

                if total_counts[focal] == 0:
                    continue

                out_meta.write(
                    f"{chrom}\t"
                    f"{pos}\t"
                    f"{gene}\t"
                    f"{focal}\n"
                )

                for pop, A, T, C, G, depth in per_pop:

                    if depth < MIN_DEPTH:
                        continue

                    focal_count = {
                        "A": A,
                        "T": T,
                        "C": C,
                        "G": G
                    }[focal]

                    af = focal_count / depth

                    out_af.write(
                        f"{chrom}\t"
                        f"{pos}\t"
                        f"{gene}\t"
                        f"{pop}\t"
                        f"{af:.8f}\t"
                        f"{depth}\t"
                        f"{focal}\n"
                    )


print(
    "[OK] wrote:",
    OUT_AF,
    file=sys.stderr
)

print(
    "[OK] wrote:",
    OUT_META,
    file=sys.stderr
)
