#!/usr/bin/env python3

import csv
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path("/mnt/spareHD_2/genomewide_codeml_kuster")

FASTA_DIR = (
    ROOT
    / "05_gene_alignments"
    / "complete_codon_fasta"
)

SUMMARY = (
    ROOT
    / "05_gene_alignments"
    / "qa"
    / "gene_alignment_summary.with_class.tsv"
)

MAIN_LIST = (
    ROOT
    / "05_gene_alignments"
    / "qa"
    / "codeml_main_genes.txt"
)

TREE = Path(
    "/work/cyu/Kuster2026_genomewide_reanalysis/"
    "genomewide_erc/00_inputs/"
    "erc_master_chr21_27_unrooted_topology.nwk"
)

OUT = ROOT / "06_codeml"
PHYLIP_DIR = OUT / "phylip"
PILOT_DIR = OUT / "pilot"
QA_DIR = OUT / "qa"

EXPECTED_SAMPLES = 27
EXPECTED_MAIN_GENES = 17965

PILOT_QUOTA = {
    "direct_n-mt": 5,
    "indirect_n-mt": 5,
    "non-n-mt": 10
}


def read_fasta(path):
    records = []
    name = None
    parts = []

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))

                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())

        if name is not None:
            records.append((name, "".join(parts)))

    return records


def write_phylip(path, records):
    lengths = {len(sequence) for _, sequence in records}

    if len(lengths) != 1:
        raise ValueError(f"Unequal alignment lengths: {path}")

    length = lengths.pop()
    temporary = Path(str(path) + ".tmp")

    with open(temporary, "w") as handle:
        handle.write(f"{len(records)} {length}\n")

        for name, sequence in records:
            handle.write(f"{name}  {sequence}\n")

    temporary.replace(path)


def alignment_variation(records):
    sequences = [sequence for _, sequence in records]
    length = len(sequences[0])

    variable_nt = 0
    variable_codons = 0

    for position in range(length):
        bases = {sequence[position] for sequence in sequences}

        if len(bases) > 1:
            variable_nt += 1

    for start in range(0, length, 3):
        codons = {
            sequence[start:start + 3]
            for sequence in sequences
        }

        if len(codons) > 1:
            variable_codons += 1

    return variable_nt, variable_codons


def write_ctl(path):
    text = """seqfile = alignment.phy
treefile = tree.nwk
outfile = mlc

noisy = 0
verbose = 0
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 0
NSsites = 0
icode = 0

fix_kappa = 0
kappa = 2

fix_omega = 0
omega = 0.2

fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8

getSE = 0
RateAncestor = 0
Small_Diff = 5e-7
cleandata = 1
method = 0
"""

    with open(path, "w") as handle:
        handle.write(text)


def main():
    for directory in (OUT, PHYLIP_DIR, PILOT_DIR, QA_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    for required in (SUMMARY, MAIN_LIST, TREE):
        if not required.is_file():
            raise FileNotFoundError(f"Missing input: {required}")

    with open(MAIN_LIST) as handle:
        main_genes = [
            line.strip()
            for line in handle
            if line.strip()
        ]

    if len(main_genes) != EXPECTED_MAIN_GENES:
        raise ValueError(
            f"Expected {EXPECTED_MAIN_GENES} main genes; "
            f"found {len(main_genes)}"
        )

    with open(SUMMARY) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        metadata = {
            row["gene_id"]: row
            for row in reader
        }

    tree_text = TREE.read_text().strip()

    tree_tips = set(re.findall(
        r"(?<=[(,])([A-Za-z0-9_]+)(?=[:),;])",
        tree_text
    ))

    qa_file = QA_DIR / "codeml_input_QA.tsv"
    eligible_by_class = defaultdict(list)

    with open(qa_file, "w", newline="") as qa_handle:
        writer = csv.writer(
            qa_handle,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "gene_id",
            "Kuster_class",
            "n_samples",
            "alignment_nt",
            "alignment_codons",
            "variable_nt_sites",
            "variable_codon_sites",
            "codeml_input_status"
        ])

        for index, gene_id in enumerate(main_genes, start=1):
            fasta = (
                FASTA_DIR
                / f"{gene_id}.complete_codons.fa"
            )

            if not fasta.is_file():
                raise FileNotFoundError(
                    f"Missing alignment: {fasta}"
                )

            records = read_fasta(fasta)

            if len(records) != EXPECTED_SAMPLES:
                raise ValueError(
                    f"{gene_id}: expected 27 sequences, "
                    f"found {len(records)}"
                )

            sample_names = {name for name, _ in records}

            if tree_tips and sample_names != tree_tips:
                missing_tree = sorted(sample_names - tree_tips)
                missing_alignment = sorted(tree_tips - sample_names)

                raise ValueError(
                    f"{gene_id}: tree/alignment tip mismatch; "
                    f"alignment_only={missing_tree}; "
                    f"tree_only={missing_alignment}"
                )

            lengths = {len(sequence) for _, sequence in records}

            if len(lengths) != 1:
                raise ValueError(
                    f"{gene_id}: unequal sequence lengths"
                )

            alignment_nt = lengths.pop()

            if alignment_nt % 3 != 0:
                raise ValueError(
                    f"{gene_id}: non-triplet alignment"
                )

            variable_nt, variable_codons = (
                alignment_variation(records)
            )

            status = (
                "variable"
                if variable_codons > 0
                else "invariant"
            )

            phylip = PHYLIP_DIR / f"{gene_id}.phy"
            write_phylip(phylip, records)

            gene_class = metadata[gene_id]["Kuster_class"]

            writer.writerow([
                gene_id,
                gene_class,
                len(records),
                alignment_nt,
                alignment_nt // 3,
                variable_nt,
                variable_codons,
                status
            ])

            if variable_codons > 0:
                eligible_by_class[gene_class].append((
                    gene_id,
                    variable_codons,
                    alignment_nt // 3
                ))

            if index % 500 == 0:
                print(
                    f"[progress] {index}/"
                    f"{EXPECTED_MAIN_GENES}",
                    flush=True
                )

    pilot_genes = []

    for gene_class, quota in PILOT_QUOTA.items():
        candidates = eligible_by_class[gene_class]

        # Avoid testing only extreme genes: prioritize alignments
        # with at least 200 codons and moderate variation.
        candidates.sort(
            key=lambda item: (
                item[2] < 200,
                abs(item[1] - 20),
                item[0]
            )
        )

        selected = candidates[:quota]

        if len(selected) != quota:
            raise ValueError(
                f"Insufficient pilot candidates for {gene_class}"
            )

        for gene_id, variable_codons, alignment_codons in selected:
            pilot_genes.append((
                gene_id,
                gene_class,
                variable_codons,
                alignment_codons
            ))

    pilot_manifest = PILOT_DIR / "pilot_manifest.tsv"

    with open(pilot_manifest, "w", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "gene_id",
            "Kuster_class",
            "variable_codon_sites",
            "alignment_codons"
        ])

        for (
            gene_id,
            gene_class,
            variable_codons,
            alignment_codons
        ) in pilot_genes:
            gene_dir = PILOT_DIR / gene_id
            gene_dir.mkdir(exist_ok=True)

            shutil.copy2(
                PHYLIP_DIR / f"{gene_id}.phy",
                gene_dir / "alignment.phy"
            )

            shutil.copy2(
                TREE,
                gene_dir / "tree.nwk"
            )

            write_ctl(gene_dir / "codeml.ctl")

            writer.writerow([
                gene_id,
                gene_class,
                variable_codons,
                alignment_codons
            ])

    marker = OUT / "STEP06_PREPARATION_COMPLETE.txt"

    with open(marker, "w") as handle:
        handle.write(
            f"main_genes\t{len(main_genes)}\n"
        )
        handle.write(
            f"phylip_files\t"
            f"{len(list(PHYLIP_DIR.glob('*.phy')))}\n"
        )
        handle.write(
            f"pilot_genes\t{len(pilot_genes)}\n"
        )

    print()
    print("===== Step 06 preparation complete =====")
    print(f"Main genes: {len(main_genes)}")
    print(
        "PHYLIP files: "
        f"{len(list(PHYLIP_DIR.glob('*.phy')))}"
    )
    print(f"Pilot genes: {len(pilot_genes)}")
    print(f"QA: {qa_file}")
    print(f"Pilot manifest: {pilot_manifest}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
