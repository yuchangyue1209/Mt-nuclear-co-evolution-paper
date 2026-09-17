#!/usr/bin/env python3
"""Extract N and S site counts from the first branch row in each codeml mlc."""

import csv
import re
from pathlib import Path

ROOT = Path("/path/to/data/genomewide_codeml_kuster")
RESULTS = ROOT / "07_codeml_genomewide" / "results"
RUN_LIST = ROOT / "07_codeml_genomewide" / "codeml_variable_genes.txt"
OUTPUT = ROOT / "07_codeml_genomewide" / "codeml_site_information.tsv"

EXPECTED = 15241
HEADER = re.compile(r"^\s*branch\s+t\s+N\s+S\s+dN/dS")
ROW = re.compile(
    r"^\s*\S+\s+[-+0-9.eE]+\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+"
)


def parse_mlc(path: Path):
    found_header = False
    with path.open(errors="replace") as handle:
        for line in handle:
            if not found_header:
                found_header = bool(HEADER.search(line))
                continue
            match = ROW.search(line)
            if match:
                return float(match.group(1)), float(match.group(2))
    return None, None


def main():
    genes = [x.strip() for x in RUN_LIST.read_text().splitlines() if x.strip()]
    if len(genes) != EXPECTED:
        raise SystemExit(f"Expected {EXPECTED} genes; found {len(genes)}")

    parsed = 0
    with OUTPUT.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "N_sites", "S_sites", "site_parse_status"])
        for index, gene in enumerate(genes, 1):
            n_sites, s_sites = parse_mlc(RESULTS / gene / "mlc")
            status = "parsed" if n_sites is not None else "failed"
            parsed += status == "parsed"
            writer.writerow([gene, n_sites or "NA", s_sites or "NA", status])
            if index % 1000 == 0:
                print(f"[progress] {index}/{len(genes)}", flush=True)

    if parsed != EXPECTED:
        raise SystemExit(f"Only {parsed}/{EXPECTED} files parsed")
    print(f"[OK] {parsed} genes -> {OUTPUT}")


if __name__ == "__main__":
    main()

