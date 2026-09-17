#!/usr/bin/env python3

import os
import re
import math
import pandas as pd
from collections import Counter

BASE = "/work/cyu/provean_inputs_OXPHOS72"
ALIGN_DIR = "/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
INFILE = f"{BASE}/unique_AA_substitutions_by_population.tsv"
OUTFILE = f"{BASE}/AA_substitution_conservation_BLOSUM.tsv"

# BLOSUM62 matrix
BLOSUM62 = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""

def parse_blosum(mat):
    lines = [l.split() for l in mat.strip().splitlines()]
    header = lines[0]
    d = {}
    for row in lines[1:]:
        aa = row[0]
        for b, val in zip(header, row[1:]):
            d[(aa, b)] = int(val)
    return d

B = parse_blosum(BLOSUM62)

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
                if name:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        seqs[name] = "".join(chunks)
    return seqs

def parse_variant(v):
    m = re.match(r"^([A-Z])([0-9]+)([A-Z])$", v)
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)

def aligned_col_for_ref_pos(ref_seq, pos):
    count = 0
    for i, aa in enumerate(ref_seq):
        if aa != "-":
            count += 1
        if count == pos:
            return i
    return None

def shannon_entropy(chars):
    n = len(chars)
    if n == 0:
        return None
    c = Counter(chars)
    ent = 0
    for x in c.values():
        p = x / n
        ent -= p * math.log2(p)
    return ent

df = pd.read_csv(INFILE, sep="\t")
rows = []

for _, r in df.iterrows():
    gene = r["gene"]
    region = r["region"]
    variant = r["variant"]
    n_pops = r["n_pops"]
    populations = r["populations"]

    parsed = parse_variant(variant)
    if parsed is None:
        continue
    ref_aa, pos, alt_aa = parsed

    aln_path = os.path.join(ALIGN_DIR, gene, f"{gene}.pep.aln.faa")
    if not os.path.exists(aln_path):
        continue

    seqs = read_fasta(aln_path)

    marine = "RS" if region == "AK" else "SAY"
    ref_name = None
    for name in seqs:
        if marine in re.split(r"[_\\.\\-]+", name):
            ref_name = name
            break

    if ref_name is None:
        ref_name = list(seqs.keys())[0]

    ref_seq = seqs[ref_name]
    col = aligned_col_for_ref_pos(ref_seq, pos)
    if col is None:
        continue

    actual_ref = ref_seq[col]
    site_aas = [
        s[col] for s in seqs.values()
        if col < len(s) and s[col] not in {"-", "X", "?", "*"}
    ]

    total = len(site_aas)
    counts = Counter(site_aas)
    major_aa, major_count = counts.most_common(1)[0] if total > 0 else ("NA", 0)

    ref_freq = counts.get(ref_aa, 0) / total if total else None
    alt_freq = counts.get(alt_aa, 0) / total if total else None
    major_freq = major_count / total if total else None
    entropy = shannon_entropy(site_aas)

    blosum = B.get((ref_aa, alt_aa), None)

    # simple category
    if major_freq is not None and major_freq >= 0.90 and blosum is not None and blosum <= 0:
        impact = "high_conservation_radical"
    elif major_freq is not None and major_freq >= 0.75:
        impact = "conserved_site"
    elif blosum is not None and blosum <= -1:
        impact = "radical_change"
    else:
        impact = "weak_or_uncertain"

    rows.append({
        "gene": gene,
        "region": region,
        "variant": variant,
        "n_pops": n_pops,
        "populations": populations,
        "marine_ref": marine,
        "actual_ref_at_pos": actual_ref,
        "alignment_column_1based": col + 1,
        "n_sequences_scored": total,
        "major_aa": major_aa,
        "major_freq": round(major_freq, 3) if major_freq is not None else "NA",
        "ref_aa_freq": round(ref_freq, 3) if ref_freq is not None else "NA",
        "alt_aa_freq": round(alt_freq, 3) if alt_freq is not None else "NA",
        "entropy": round(entropy, 3) if entropy is not None else "NA",
        "BLOSUM62": blosum,
        "impact_proxy": impact
    })

out = pd.DataFrame(rows)
out = out.sort_values(["impact_proxy", "BLOSUM62", "major_freq"], ascending=[True, True, False])
out.to_csv(OUTFILE, sep="\t", index=False)

print(out.to_string(index=False))
print("\nSaved:", OUTFILE)
