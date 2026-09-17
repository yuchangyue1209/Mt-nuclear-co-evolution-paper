#!/usr/bin/env python3
"""Remove SC, CH, LB, PACH, and FRED from codon alignments."""
import argparse
from pathlib import Path

EXCLUDED = {"SC", "CH", "LB", "PACH", "FRED"}

def read_fasta(path):
    records, name, sequence = [], None, []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: records.append((name, "".join(sequence)))
            name, sequence = line[1:].split()[0], []
        else: sequence.append(line)
    if name is not None: records.append((name, "".join(sequence)))
    return records

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True); p.add_argument("--output-dir", required=True)
    a = p.parse_args(); indir, outdir = Path(a.input_dir), Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    files = sorted(indir.glob("*.fa")) + sorted(indir.glob("*.fasta"))
    if not files: raise SystemExit("No FASTA alignments found")
    audit = [("alignment", "input_tips", "retained_tips")]
    for path in files:
        records = read_fasta(path)
        retained = [(n, s) for n, s in records if n.upper() not in EXCLUDED]
        if len(records) != 27 or len(retained) != 22:
            raise ValueError(f"{path.name}: expected 27 input and 22 retained tips")
        lengths = {len(s) for _, s in retained}
        if len(lengths) != 1 or next(iter(lengths)) % 3:
            raise ValueError(f"{path.name}: invalid codon alignment")
        with (outdir / path.name).open("w") as h:
            for name, sequence in retained: h.write(f">{name}\n{sequence}\n")
        audit.append((path.name, 27, 22))
    with (outdir / "alignment_filter_audit.tsv").open("w") as h:
        for row in audit: h.write("\t".join(map(str, row)) + "\n")

if __name__ == "__main__": main()
