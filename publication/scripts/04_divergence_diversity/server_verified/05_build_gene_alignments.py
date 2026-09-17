#!/usr/bin/env python3

import csv
import gzip
import re
import sys
from itertools import zip_longest
from pathlib import Path


ROOT = Path("/path/to/data/genomewide_codeml_kuster")

CDS_DIR = ROOT / "04_consensus_cds" / "canonical_cds"
GTF = (
    ROOT
    / "00_targets_codeml_ready"
    / "stickleback_20347_unique_codeml.gtf"
)

OUT = ROOT / "05_gene_alignments"
MASKED_DIR = OUT / "masked_fasta"
STRICT_DIR = OUT / "complete_codon_fasta"
QA_DIR = OUT / "qa"

EXPECTED_SAMPLES = 27
EXPECTED_GENES = 20347

# This flag is descriptive only; no gene is deleted in Step 05.
MIN_COMPLETE_CODONS = 100
MIN_COMPLETE_FRACTION = 0.50

VALID_BASES = set("ACGT")
STOP_CODONS = {"TAA", "TAG", "TGA"}

transcript_re = re.compile(r'transcript_id "([^"]+)"')
gene_re = re.compile(r'gene_id "([^"]+)"')


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def fasta_iterator(path):
    """Stream one FASTA record at a time."""
    with open_text(path) as handle:
        header = None
        sequence_parts = []

        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:].split()[0]
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence encountered before FASTA header: {path}"
                    )

                sequence_parts.append(line)

        if header is not None:
            yield header, "".join(sequence_parts)


def read_gtf_transcript_map(gtf_path):
    transcript_to_gene = {}
    gene_to_transcript = {}

    with open(gtf_path, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                raise ValueError(
                    f"Invalid GTF row at line {line_number}: "
                    f"expected 9 fields, observed {len(fields)}"
                )

            if fields[2] != "transcript":
                continue

            attributes = fields[8]

            transcript_match = transcript_re.search(attributes)
            gene_match = gene_re.search(attributes)

            if transcript_match is None or gene_match is None:
                raise ValueError(
                    f"Missing transcript_id or gene_id at "
                    f"GTF line {line_number}"
                )

            transcript_id = transcript_match.group(1)
            gene_id = gene_match.group(1)

            if transcript_id in transcript_to_gene:
                raise ValueError(
                    f"Duplicated transcript ID in GTF: {transcript_id}"
                )

            if gene_id in gene_to_transcript:
                raise ValueError(
                    f"More than one selected transcript for gene: {gene_id}"
                )

            transcript_to_gene[transcript_id] = gene_id
            gene_to_transcript[gene_id] = transcript_id

    return transcript_to_gene, gene_to_transcript


def write_fasta(path, records):
    temporary = Path(str(path) + ".tmp")

    with open(temporary, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")

            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    temporary.replace(path)


def main():
    for directory in (OUT, MASKED_DIR, STRICT_DIR, QA_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    if not GTF.is_file():
        raise FileNotFoundError(f"GTF not found: {GTF}")

    cds_files = sorted(CDS_DIR.glob("*.canonical_cds.fa"))

    if len(cds_files) != EXPECTED_SAMPLES:
        raise ValueError(
            f"Expected {EXPECTED_SAMPLES} CDS FASTA files, "
            f"found {len(cds_files)}"
        )

    samples = [
        path.name.replace(".canonical_cds.fa", "")
        for path in cds_files
    ]

    if len(set(samples)) != EXPECTED_SAMPLES:
        raise ValueError("Duplicated sample names detected")

    print("===== Step 05: gene-wise codon alignments =====")
    print(f"[input] CDS files: {len(cds_files)}")
    print(f"[input] GTF: {GTF}")
    print(f"[samples] {','.join(samples)}")
    print(f"[output] {OUT}")

    transcript_to_gene, gene_to_transcript = (
        read_gtf_transcript_map(GTF)
    )

    print(
        f"[GTF] Selected transcripts: "
        f"{len(transcript_to_gene)}"
    )
    print(
        f"[GTF] Selected genes: "
        f"{len(gene_to_transcript)}"
    )

    if len(transcript_to_gene) != EXPECTED_GENES:
        raise ValueError(
            f"Expected {EXPECTED_GENES} GTF transcripts, "
            f"found {len(transcript_to_gene)}"
        )

    iterators = [
        fasta_iterator(path)
        for path in cds_files
    ]

    summary_path = QA_DIR / "gene_alignment_summary.tsv"
    manifest_path = OUT / "gene_alignment_manifest.tsv"
    sample_summary_path = QA_DIR / "sample_masking_summary.tsv"

    sample_ambiguous_codons = {
        sample: 0 for sample in samples
    }
    sample_internal_stops = {
        sample: 0 for sample in samples
    }
    sample_terminal_stops = {
        sample: 0 for sample in samples
    }
    sample_total_codons = {
        sample: 0 for sample in samples
    }

    observed_transcripts = set()
    genes_processed = 0
    recommended_genes = 0

    with (
        open(summary_path, "w", newline="") as summary_handle,
        open(manifest_path, "w", newline="") as manifest_handle
    ):
        summary_writer = csv.writer(
            summary_handle,
            delimiter="\t",
            lineterminator="\n"
        )

        manifest_writer = csv.writer(
            manifest_handle,
            delimiter="\t",
            lineterminator="\n"
        )

        summary_writer.writerow([
            "gene_id",
            "transcript_id",
            "n_samples",
            "cds_length_nt",
            "total_codons",
            "complete_codons",
            "removed_codons",
            "complete_codon_fraction",
            "samples_with_ambiguous_codons",
            "ambiguous_sample_codons",
            "internal_stop_sample_codons",
            "terminal_stop_sample_codons",
            "codeml_recommended"
        ])

        manifest_writer.writerow([
            "gene_id",
            "transcript_id",
            "masked_alignment",
            "complete_codon_alignment"
        ])

        grouped_records = zip_longest(
            *iterators,
            fillvalue=None
        )

        for records in grouped_records:
            if any(record is None for record in records):
                raise ValueError(
                    "CDS FASTA files contain unequal record counts"
                )

            transcript_ids = [
                record[0] for record in records
            ]

            if len(set(transcript_ids)) != 1:
                details = ", ".join(
                    f"{sample}={transcript_id}"
                    for sample, transcript_id
                    in zip(samples, transcript_ids)
                )

                raise ValueError(
                    "CDS record order differs among samples: "
                    + details
                )

            transcript_id = transcript_ids[0]

            if transcript_id in observed_transcripts:
                raise ValueError(
                    f"Duplicated CDS transcript: {transcript_id}"
                )

            observed_transcripts.add(transcript_id)

            if transcript_id not in transcript_to_gene:
                raise ValueError(
                    f"Transcript absent from GTF mapping: "
                    f"{transcript_id}"
                )

            gene_id = transcript_to_gene[transcript_id]

            sequences = [
                record[1].upper()
                for record in records
            ]

            lengths = {len(sequence) for sequence in sequences}

            if len(lengths) != 1:
                length_details = ", ".join(
                    f"{sample}={len(sequence)}"
                    for sample, sequence
                    in zip(samples, sequences)
                )

                raise ValueError(
                    f"Unequal sequence lengths for {gene_id}: "
                    f"{length_details}"
                )

            cds_length = lengths.pop()

            if cds_length % 3 != 0:
                raise ValueError(
                    f"Non-triplet CDS length for {gene_id}: "
                    f"{cds_length}"
                )

            total_codons = cds_length // 3

            masked_codons_by_sample = [
                [] for _ in samples
            ]

            strict_codons_by_sample = [
                [] for _ in samples
            ]

            ambiguous_samples_for_gene = set()
            ambiguous_sample_codons = 0
            internal_stop_sample_codons = 0
            terminal_stop_sample_codons = 0
            complete_codons = 0

            for codon_index in range(total_codons):
                start = codon_index * 3

                original_codons = [
                    sequence[start:start + 3]
                    for sequence in sequences
                ]

                masked_column = []
                complete_column = True

                for sample_index, (sample, codon) in enumerate(
                    zip(samples, original_codons)
                ):
                    sample_total_codons[sample] += 1

                    contains_ambiguity = (
                        len(codon) != 3
                        or any(base not in VALID_BASES for base in codon)
                    )

                    is_stop = codon in STOP_CODONS
                    is_terminal = codon_index == total_codons - 1

                    if contains_ambiguity:
                        masked_codon = "NNN"
                        complete_column = False
                        ambiguous_samples_for_gene.add(sample)
                        ambiguous_sample_codons += 1
                        sample_ambiguous_codons[sample] += 1

                    elif is_stop:
                        masked_codon = "NNN"
                        complete_column = False

                        if is_terminal:
                            terminal_stop_sample_codons += 1
                            sample_terminal_stops[sample] += 1
                        else:
                            internal_stop_sample_codons += 1
                            sample_internal_stops[sample] += 1
                    else:
                        masked_codon = codon

                    masked_codons_by_sample[sample_index].append(
                        masked_codon
                    )
                    masked_column.append(masked_codon)

                if complete_column:
                    complete_codons += 1

                    for sample_index, codon in enumerate(masked_column):
                        strict_codons_by_sample[sample_index].append(
                            codon
                        )

            removed_codons = total_codons - complete_codons

            if total_codons > 0:
                complete_fraction = complete_codons / total_codons
            else:
                complete_fraction = 0.0

            recommended = (
                complete_codons >= MIN_COMPLETE_CODONS
                and complete_fraction >= MIN_COMPLETE_FRACTION
            )

            if recommended:
                recommended_genes += 1

            masked_records = [
                (sample, "".join(codons))
                for sample, codons
                in zip(samples, masked_codons_by_sample)
            ]

            strict_records = [
                (sample, "".join(codons))
                for sample, codons
                in zip(samples, strict_codons_by_sample)
            ]

            masked_path = MASKED_DIR / f"{gene_id}.masked.fa"
            strict_path = STRICT_DIR / f"{gene_id}.complete_codons.fa"

            write_fasta(masked_path, masked_records)
            write_fasta(strict_path, strict_records)

            summary_writer.writerow([
                gene_id,
                transcript_id,
                len(samples),
                cds_length,
                total_codons,
                complete_codons,
                removed_codons,
                f"{complete_fraction:.8f}",
                len(ambiguous_samples_for_gene),
                ambiguous_sample_codons,
                internal_stop_sample_codons,
                terminal_stop_sample_codons,
                "yes" if recommended else "no"
            ])

            manifest_writer.writerow([
                gene_id,
                transcript_id,
                str(masked_path),
                str(strict_path)
            ])

            genes_processed += 1

            if genes_processed % 500 == 0:
                print(
                    f"[progress] {genes_processed}/"
                    f"{EXPECTED_GENES}",
                    flush=True
                )

    if genes_processed != EXPECTED_GENES:
        raise ValueError(
            f"Expected {EXPECTED_GENES} alignments, "
            f"created {genes_processed}"
        )

    missing_from_fasta = (
        set(transcript_to_gene) - observed_transcripts
    )

    extra_in_fasta = (
        observed_transcripts - set(transcript_to_gene)
    )

    if missing_from_fasta:
        raise ValueError(
            f"Transcripts missing from FASTAs: "
            f"{len(missing_from_fasta)}"
        )

    if extra_in_fasta:
        raise ValueError(
            f"Unexpected transcripts in FASTAs: "
            f"{len(extra_in_fasta)}"
        )

    with open(sample_summary_path, "w", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "sample",
            "total_codons",
            "ambiguous_codons",
            "ambiguous_codon_fraction",
            "internal_stop_codons",
            "terminal_stop_codons"
        ])

        for sample in samples:
            total = sample_total_codons[sample]
            ambiguous = sample_ambiguous_codons[sample]

            fraction = (
                ambiguous / total if total > 0 else 0.0
            )

            writer.writerow([
                sample,
                total,
                ambiguous,
                f"{fraction:.8f}",
                sample_internal_stops[sample],
                sample_terminal_stops[sample]
            ])

    marker = OUT / "STEP05_COMPLETE.txt"

    with open(marker, "w") as handle:
        handle.write(
            f"genes_processed\t{genes_processed}\n"
        )
        handle.write(
            f"samples\t{len(samples)}\n"
        )
        handle.write(
            f"recommended_genes\t{recommended_genes}\n"
        )

    print()
    print("===== Final audit =====")
    print(f"Samples: {len(samples)}")
    print(f"Genes processed: {genes_processed}")
    print(f"Masked alignments: {len(list(MASKED_DIR.glob('*.masked.fa')))}")
    print(
        "Complete-codon alignments: "
        f"{len(list(STRICT_DIR.glob('*.complete_codons.fa')))}"
    )
    print(f"Preliminary recommended genes: {recommended_genes}")
    print(f"[summary] {summary_path}")
    print(f"[sample QA] {sample_summary_path}")
    print(f"[manifest] {manifest_path}")
    print("[05] COMPLETE")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
