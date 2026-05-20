#provean
#/work/cyu/make_provean_inputs.py
#!/usr/bin/env python3

import os
import re
from collections import defaultdict

# ============================================================
# Paths
# ============================================================
ALIGN_DIR = "/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
OUTDIR = "/work/cyu/provean_inputs_oxphos_candidates"
os.makedirs(OUTDIR, exist_ok=True)

# candidate genes; 可以改成你想看的
CANDIDATE_GENES = [
    "cox4i2",
    "cox7a2l",
    "uqcrc2b",
    "ndufb3",
    "atp5mc1",
    "atp5pd",
    "ndufs3",
    "uqcrc1",
    "ndufs7",
    "ndufa7",
    "ndufs4",
    "cox4i1"
]

AK_FW = {"FG","LG","SR","SL","TL","WB","WT","WK","LB"}
BC_FW = {"SWA","THE","JOE","BEA","MUC","PYE","BOOT","ECHO","LAW","GOS","ROB"}

MARINE_REF = {
    "AK": "RS",
    "BC": "SAY"
}

ALL_POPS = AK_FW | BC_FW | {"RS","SAY","AMO","PACH","FRED","SC","CH"}

BAD_AA = set(["-", "X", "N", "*", "?"])


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
    """
    Handles sample names like:
    1_FG
    1_FG_S1
    output_1_FG_S1
    25_RS
    17_SAY
    """
    parts = re.split(r"[_\.\-]+", sample_name)

    for p in parts:
        if p in ALL_POPS:
            return p

    # fallback: direct contains
    for p in sorted(ALL_POPS, key=len, reverse=True):
        if re.search(rf"(^|[_\.\-]){p}([_\.\-]|$)", sample_name):
            return p

    return None


def aligned_to_ungapped_position(ref_aln):
    """
    For each aligned index, return ungapped ref amino-acid position.
    If ref has gap at that alignment column, return None.
    Position is 1-based.
    """
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
    """
    Return variants like A123T relative to ungapped ref.
    Ignore:
    - ref gaps
    - query gaps
    - ambiguous residues
    - identical residues
    """
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
# Main
# ============================================================
summary_rows = []

for gene in CANDIDATE_GENES:
    aln_file = os.path.join(ALIGN_DIR, gene, f"{gene}.pep.aln.faa")

    if not os.path.exists(aln_file):
        print(f"[skip] missing alignment: {aln_file}")
        continue

    seqs = read_fasta(aln_file)

    # map pop -> sample header
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

        # PROVEAN reference protein fasta
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

        # unique substitutions for PROVEAN
        var_file = os.path.join(out_region, f"{gene}_{region}_variants.txt")
        with open(var_file, "w") as out:
            for v in sorted(variant_to_pops):
                out.write(v + "\n")

        # detailed map: variant -> populations
        map_file = os.path.join(out_region, f"{gene}_{region}_variant_population_map.tsv")
        with open(map_file, "w") as out:
            out.write("gene\tregion\tmarine_ref\tvariant\tpopulations\tn_pops\n")
            for v in sorted(variant_to_pops):
                pops = sorted(set(variant_to_pops[v]))
                out.write(
                    f"{gene}\t{region}\t{marine_pop}\t{v}\t{','.join(pops)}\t{len(pops)}\n"
                )

        print(f"[ok] {gene} {region}")
        print(f"  ref fasta : {ref_fa}")
        print(f"  variants  : {var_file}")
        print(f"  map       : {map_file}")

# write summary
summary_file = os.path.join(OUTDIR, "candidate_gene_AA_substitution_summary.tsv")
with open(summary_file, "w") as out:
    out.write("gene\tregion\tmarine_ref\tpop\tsample\tn_variants\tvariants\n")
    for r in summary_rows:
        out.write(
            f"{r['gene']}\t{r['region']}\t{r['marine_ref']}\t{r['pop']}\t"
            f"{r['sample']}\t{r['n_variants']}\t{r['variants']}\n"
        )

print("\nDONE")
print(f"Summary: {summary_file}")
print(f"Output : {OUTDIR}")


awk -F'\t' 'NR>1 {sum[$1]+=$6} END{for(g in sum) print g, sum[g]}' \
/work/cyu/provean_inputs_oxphos_candidates/candidate_gene_AA_substitution_summary.tsv | sort -k2,2nr

awk -F'\t' 'NR>1 && $6>0 {
  split($7,a,",")
  for(i in a){
    print $1"\t"$2"\t"a[i]"\t"$4
  }
}' /work/cyu/provean_inputs_oxphos_candidates/candidate_gene_AA_substitution_summary.tsv \
| sort \
| awk -F'\t' '
{
  key=$1 FS $2 FS $3
  pops[key]=pops[key]","$4
  n[key]++
}
END{
  print "gene\tregion\tvariant\tn_pops\tpopulations"
  for(k in n){
    p=pops[k]; sub(/^,/,"",p)
    print k"\t"n[k]"\t"p
  }
}' OFS='\t' \
| sort -k1,1 -k2,2 -k4,4nr \
> /work/cyu/provean_inputs_oxphos_candidates/unique_AA_substitutions_by_population.tsv

head -n 30 /work/cyu/provean_inputs_oxphos_candidates/unique_AA_substitutions_by_population.tsv



mkdir -p /work/cyu/provean_db
cd /work/cyu/provean_db
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz








#
cat score_AA_substitutions_conservation.py
#!/usr/bin/env python3

import os
import re
import math
import pandas as pd
from collections import Counter

BASE = "/work/cyu/provean_inputs_oxphos_candidates"
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



/work/cyu/provean_inputs_oxphos_candidates/AA_substitution_conservation_BLOSUM.tsv


#sift candidates
cd /work/cyu/provean_inputs_oxphos_candidates
cat run_sift4g_candidates.py
#!/usr/bin/env python3

import os
import re
import subprocess
import pandas as pd
from pathlib import Path

BASE = Path("/work/cyu/provean_inputs_oxphos_candidates")
DB = "/work/cyu/provean_db/custom_oxphos/custom_oxphos_clean_unique.faa"
SIFT4G = "/home/cyu/.conda/envs/provean_env/bin/sift4g"

INFILE = BASE / "unique_AA_substitutions_by_population.tsv"
OUTBASE = BASE / "sift4g_all_candidates"
OUTBASE.mkdir(exist_ok=True)

df = pd.read_csv(INFILE, sep="\t")

rows = []

for _, r in df.iterrows():
    gene = r["gene"]
    region = r["region"]
    variant = r["variant"]

    d = BASE / f"{gene}_{region}"
    ref_files = list(d.glob("*_ref_*.faa"))

    if not ref_files:
        print(f"[skip] missing ref fasta: {gene}_{region}")
        continue

    ref = ref_files[0]

    # fasta header determines output prefix and subst filename
    with open(ref) as f:
        header = next(line[1:].strip().split()[0] for line in f if line.startswith(">"))

    work = OUTBASE / f"{gene}_{region}_{variant}"
    subst_dir = work / "subst"
    outdir = work / "out"
    subst_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    subst_file = subst_dir / f"{header}.subst"
    subst_file.write_text(f"{variant}\n")

    cmd = [
        SIFT4G,
        "-q", str(ref),
        "-d", DB,
        "--subst", str(subst_dir),
        "--out", str(outdir),
        "-t", "4"
    ]

    print("[run]", gene, region, variant)
    subprocess.run(cmd, check=True)

    pred_file = outdir / f"{header}.SIFTprediction"
    if not pred_file.exists():
        print(f"[warn] missing prediction: {pred_file}")
        continue

    score = None
    pred = "NA"

    with open(pred_file) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            fields = line.split()

            # expected:
            # A32V TOLERATED 1.00 4.32 5 59

            if len(fields) >= 3 and fields[0] == variant:
                pred = fields[1].lower()
                score = float(fields[2])
                break
    rows.append({
        "gene": gene,
        "region": region,
        "variant": variant,
        "n_pops": r["n_pops"],
        "populations": r["populations"],
        "marine_ref": r.get("marine_ref", ""),
        "SIFT_score": score,
        "SIFT_prediction": pred,
        "prediction_file": str(pred_file)
    })

out = pd.DataFrame(rows)
outfile = BASE / "SIFT4G_candidate_substitution_results.tsv"
out.to_csv(outfile, sep="\t", index=False)

print("\nDONE")
print("Saved:", outfile)
print(out.to_string(index=False))












#72 aa
nano /work/cyu/make_AA_inputs_OXPHOS72.py
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

cd /work/cyu/provean_inputs_OXPHOS72

awk -F'\t' 'NR>1 && $6>0 {
  split($7,a,",")
  for(i in a){
    print $1"\t"$2"\t"a[i]"\t"$4
  }
}' OXPHOS72_AA_substitution_summary.tsv \
| sort \
| awk -F'\t' '
BEGIN{OFS="\t"; print "gene","region","variant","n_pops","populations"}
{
  key=$1 FS $2 FS $3
  pops[key]=(key in pops ? pops[key]","$4 : $4)
  n[key]++
}
END{
  for(k in n){
    print k,n[k],pops[k]
  }
}' \
| sort -t$'\t' -k1,1 -k2,2 -k4,4nr \
> unique_AA_substitutions_by_population.tsv

cd /work/cyu/provean_inputs_OXPHOS72

(
echo -e "gene\tregion\tvariant\tn_pops\tpopulations"
cat unique_AA_substitutions_by_population.tsv
) > tmp.tsv

mv tmp.tsv unique_AA_substitutions_by_population.tsv

python3 score_AA_substitutions_conservation_OXPHOS72.py


#sift72
cd /work/cyu/provean_inputs_OXPHOS72
python3 run_sift4g_OXPHOS72.py
cd /work/cyu/provean_inputs_OXPHOS72

cat > run_sift4g_OXPHOS72.py <<'PY'
#!/usr/bin/env python3

import subprocess
import pandas as pd
from pathlib import Path

BASE = Path("/work/cyu/provean_inputs_OXPHOS72")
DB = "/work/cyu/provean_db/custom_oxphos/custom_oxphos_clean_unique.faa"
SIFT4G = "/home/cyu/.conda/envs/provean_env/bin/sift4g"

INFILE = BASE / "unique_AA_substitutions_by_population.tsv"
OUTBASE = BASE / "sift4g_all_OXPHOS72"
OUTBASE.mkdir(exist_ok=True)

df = pd.read_csv(INFILE, sep="\t")
rows = []

for _, r in df.iterrows():
    gene = str(r["gene"])
    region = str(r["region"])
    variant = str(r["variant"])

    d = BASE / f"{gene}_{region}"
    ref_files = list(d.glob("*_ref_*.faa"))

    if not ref_files:
        print(f"[skip] missing ref fasta: {gene}_{region}")
        continue

    ref = ref_files[0]

    with open(ref) as f:
        header = next(line[1:].strip().split()[0] for line in f if line.startswith(">"))

    work = OUTBASE / f"{gene}_{region}_{variant}"
    subst_dir = work / "subst"
    outdir = work / "out"
    subst_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    subst_file = subst_dir / f"{header}.subst"
    subst_file.write_text(f"{variant}\n")

    cmd = [
        SIFT4G,
        "-q", str(ref),
        "-d", DB,
        "--subst", str(subst_dir),
        "--out", str(outdir),
        "-t", "4"
    ]

    print("[run]", gene, region, variant)
    subprocess.run(cmd, check=True)

    pred_file = outdir / f"{header}.SIFTprediction"

    score = None
    pred = "NA"

    if pred_file.exists():
        with open(pred_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                fields = line.split()

                # expected format:
                # A32V TOLERATED 1.00 4.32 5 59
                if len(fields) >= 3 and fields[0] == variant:
                    pred = fields[1].lower()
                    score = float(fields[2])
                    break
    else:
        print(f"[warn] missing prediction: {pred_file}")

    rows.append({
        "gene": gene,
        "region": region,
        "variant": variant,
        "n_pops": r["n_pops"],
        "populations": r["populations"],
        "SIFT_score": score,
        "SIFT_prediction": pred,
        "prediction_file": str(pred_file)
    })

out = pd.DataFrame(rows)
outfile = BASE / "SIFT4G_candidate_substitution_results.tsv"
out.to_csv(outfile, sep="\t", index=False)

print("\nDONE")
print("Saved:", outfile)
print(out.to_string(index=False))
PY

/work/cyu/provean_inputs_OXPHOS72/AA_substitution_conservation_BLOSUM.tsv
/work/cyu/provean_inputs_OXPHOS72/SIFT4G_candidate_substitution_results.tsv





#teleos aa sift 72genes
cd /work/cyu/sift_db

awk '
BEGIN{keep=0}
/^>/{
  keep = ($0 ~ /Danio|Oryzias|Takifugu|Gasterosteus|Salmo|Xiphophorus|Gadus/)
}
keep{print}
' uniref90.fasta > teleost_uniref90.fa


python3 run_sift4g_OXPHOS72_teleostDB.py


cd /work/cyu/provean_inputs_OXPHOS72


cp run_sift4g_OXPHOS72.py run_sift4g_OXPHOS72_teleostDB.py

perl -pi -e 's|DB = ".*"|DB = "/work/cyu/provean_db/teleost_uniref90.fa"|' run_sift4g_OXPHOS72_teleostDB.py

perl -pi -e 's|OUTBASE = BASE / "sift4g_all_OXPHOS72"|OUTBASE = BASE / "sift4g_all_OXPHOS72_teleostUniRef90"|' run_sift4g_OXPHOS72_teleostDB.py

perl -pi -e 's|SIFT4G_candidate_substitution_results.tsv|SIFT4G_candidate_substitution_results_teleostUniRef90.tsv|' run_sift4g_OXPHOS72_teleostDB.py





#blosum for 287
cat > /work/cyu/make_AA_substitution_all_categories.py <<'PY'
#!/usr/bin/env python3
import os, re, math
from collections import defaultdict, Counter
import pandas as pd

OUTDIR = "/work/cyu/AA_substitution_all_categories"
os.makedirs(OUTDIR, exist_ok=True)

NU_ALIGN_DIR = "/mnt/spareHD_2/oxphos_codeml_ready/06_gene_align_72"
MT_ALIGN_DIR = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13"
ANNOT_TABLE  = "/work/cyu/codeml_sites_summary.merged.tsv"

CATEGORY_CONFIG = {
    "mtPCG": {
        "align_dir": MT_ALIGN_DIR,
        "genes": ["ND1","ND2","ND3","ND4","ND4L","ND5","ND6","CYTB","COX1","COX2","COX3","ATP6","ATP8"],
        "genetic_code": "vertebrate_mito"
    },
    "nuOXPHOS_subunit": {
        "align_dir": NU_ALIGN_DIR,
        "genes_file": "/work/cyu/gene_lists/nuOXPHOS_subunit.list",
        "genetic_code": "standard"
    },
    "OXPHOS_assembly": {
        "align_dir": NU_ALIGN_DIR,
        "genes_file": "/work/cyu/gene_lists/assembly.list",
        "genetic_code": "standard"
    },
    "nmt_ARS": {
        "align_dir": NU_ALIGN_DIR,
        "genes_file": "/work/cyu/gene_lists/nmt_ARS.list",
        "genetic_code": "standard"
    },
    "nmt_ribo": {
        "align_dir": NU_ALIGN_DIR,
        "genes_file": "/work/cyu/gene_lists/nmt_ribo.list",
        "genetic_code": "standard"
    },
    "cyto_ARS": {
        "align_dir": NU_ALIGN_DIR,
        "genes_file": "/work/cyu/gene_lists/cyto_ARS.list",
        "genetic_code": "standard"
    },
    "cyto_ribo": {
        "align_dir": NU_ALIGN_DIR,
        "genes_file": "/work/cyu/gene_lists/cyto_ribo.list",
        "genetic_code": "standard"
    }
}

AK_FW = {"FG","LG","SR","SL","TL","WB","WT","WK","LB"}
BC_FW = {"SWA","THE","JOE","BEA","MUC","PYE","BOOT","ECHO","LAW","GOS","ROB"}
MARINE_REF = {"AK": "RS", "BC": "SAY"}
ALL_POPS = AK_FW | BC_FW | {"RS","SAY","AMO","PACH","FRED","SC","CH"}
BAD_AA = {"-", "X", "N", "*", "?"}

MT_COMPLEX = {
    "ND1":"CI", "ND2":"CI", "ND3":"CI", "ND4":"CI", "ND4L":"CI", "ND5":"CI", "ND6":"CI",
    "CYTB":"CIII",
    "COX1":"CIV", "COX2":"CIV", "COX3":"CIV",
    "ATP6":"CV", "ATP8":"CV"
}

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
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0
"""

def parse_blosum(mat):
    lines = [x.split() for x in mat.strip().splitlines()]
    header = lines[0]
    d = {}
    for row in lines[1:]:
        aa1 = row[0]
        for aa2, val in zip(header, row[1:]):
            d[(aa1, aa2)] = int(val)
    return d

BLOSUM = parse_blosum(BLOSUM62)

STD_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
MT_CODE = STD_CODE.copy()
MT_CODE.update({'TGA':'W', 'ATA':'M', 'AGA':'*', 'AGG':'*'})

def clean_role(x):
    x = str(x).strip()
    if x.lower() in {"nan", "", "na", "none"}:
        return "NA"
    if x == "subunit":
        return "nuOXPHOS_subunit"
    if x == "assembly_factor":
        return "OXPHOS_assembly"
    if x.lower() == "nmt-ars":
        return "nmt_ARS"
    if x.lower() == "nmt-ribo":
        return "nmt_ribo"
    if x.lower() == "cyto-ars":
        return "cyto_ARS"
    if x.lower() == "cyto-ribo":
        return "cyto_ribo"
    return x

def clean_complex(x):
    x = str(x).strip()
    if x.lower() in {"nan", "", "na", "none"}:
        return "NA"
    x = x.replace("Complex ", "")
    mapx = {
        "I": "CI", "II": "CII", "III": "CIII", "IV": "CIV", "V": "CV",
        "CI": "CI", "CII": "CII", "CIII": "CIII", "CIV": "CIV", "CV": "CV"
    }
    return mapx.get(x, x)

def load_annotation():
    gene2category = {}
    gene2complex = {}

    if os.path.exists(ANNOT_TABLE):
        ann = pd.read_csv(ANNOT_TABLE, sep="\t", dtype=str)
        ann.columns = [c.strip() for c in ann.columns]

        ann = ann[ann["model"].astype(str) == "M0"].copy()

        for _, row in ann.iterrows():
            gene = str(row["gene"]).strip().lower()
            role = clean_role(row.get("role", "NA"))
            comp = clean_complex(row.get("complex", "NA"))

            if gene and gene != "nan":
                gene2category[gene] = role
                gene2complex[gene] = comp

    for g, comp in MT_COMPLEX.items():
        gene2category[g.lower()] = "mtPCG"
        gene2complex[g.lower()] = comp

    return gene2category, gene2complex

gene2category, gene2complex = load_annotation()

def read_gene_list(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]

def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs

def is_dna_alignment(seqs):
    chars = set("".join(seqs.values()).upper())
    return chars.issubset(set("ACGTN-?"))

def translate_aligned_codon(seq, code_name):
    code = MT_CODE if code_name == "vertebrate_mito" else STD_CODE
    clean = seq.upper().replace("?", "N")
    aa = []
    for i in range(0, len(clean) - 2, 3):
        codon = clean[i:i+3]
        if codon == "---":
            aa.append("-")
        elif "-" in codon:
            aa.append("-")
        elif "N" in codon:
            aa.append("X")
        else:
            aa.append(code.get(codon, "X"))
    return "".join(aa)

def infer_pop(sample):
    parts = re.split(r"[_\.\-\|]+", sample)
    for p in parts:
        if p in ALL_POPS:
            return p
    for p in sorted(ALL_POPS, key=len, reverse=True):
        if re.search(rf"(^|[_\.\-\|]){p}([_\.\-\|]|$)", sample):
            return p
    return None

def get_alignment_path(align_dir, gene):
    names = [gene, gene.lower(), gene.upper()]
    suffixes = [
        ".pep.aln.faa", ".aln.faa", ".faa",
        ".codon.fas", ".codon.fa", ".fas", ".fasta"
    ]
    candidates = []
    for g in names:
        for s in suffixes:
            candidates.append(os.path.join(align_dir, g, f"{g}{s}"))
            candidates.append(os.path.join(align_dir, f"{g}{s}"))
    for p in candidates:
        if os.path.exists(p):
            return p
    return None

def aligned_to_ref_pos(ref_aln):
    pos_map, pos = {}, 0
    for i, aa in enumerate(ref_aln):
        if aa != "-":
            pos += 1
            pos_map[i] = pos
        else:
            pos_map[i] = None
    return pos_map

def get_variants(ref_aln, q_aln):
    pos_map = aligned_to_ref_pos(ref_aln)
    out = []
    for i, (r, q) in enumerate(zip(ref_aln, q_aln)):
        if pos_map[i] is None:
            continue
        if r in BAD_AA or q in BAD_AA:
            continue
        if r != q:
            out.append(f"{r}{pos_map[i]}{q}")
    return out

def shannon_entropy(chars):
    n = len(chars)
    c = Counter(chars)
    ent = 0
    for v in c.values():
        p = v / n
        ent -= p * math.log2(p)
    return ent

def parse_variant(v):
    m = re.match(r"^([A-Z])([0-9]+)([A-Z])$", v)
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)

def impact_class(major_freq, blosum, n_pops):
    recurrent = n_pops >= 2
    conserved = major_freq >= 0.90
    radical = blosum is not None and blosum <= -2
    nonconservative = blosum is not None and blosum <= 0

    if recurrent and conserved and radical:
        return "recurrent_conserved_radical"
    if recurrent and conserved and nonconservative:
        return "recurrent_conserved_nonconservative"
    if conserved and radical:
        return "conserved_radical"
    if recurrent and radical:
        return "recurrent_radical"
    if recurrent:
        return "recurrent"
    if radical:
        return "radical"
    if nonconservative:
        return "nonconservative"
    return "weak_or_uncertain"

all_rows = []
summary_rows = []

for config_category, cfg in CATEGORY_CONFIG.items():
    genes = cfg["genes"] if "genes" in cfg else read_gene_list(cfg["genes_file"])
    align_dir = cfg["align_dir"]
    code_name = cfg.get("genetic_code", "standard")

    print(f"\n[category] {config_category}: {len(genes)} genes")

    for gene in genes:
        gene_key = gene.lower()
        aln_path = get_alignment_path(align_dir, gene)

        if aln_path is None:
            print(f"[skip] missing alignment: {config_category} {gene}")
            continue

        seqs_raw = read_fasta(aln_path)

        if is_dna_alignment(seqs_raw):
            seqs = {k: translate_aligned_codon(v, code_name) for k, v in seqs_raw.items()}
            seq_type = "translated_codon_alignment"
        else:
            seqs = seqs_raw
            seq_type = "protein_alignment"

        category = gene2category.get(gene_key, config_category)
        complex_ = gene2complex.get(gene_key, "NA")

        pop_to_sample = {}
        for sample in seqs:
            pop = infer_pop(sample)
            if pop:
                pop_to_sample[pop] = sample

        for region, marine in MARINE_REF.items():
            if marine not in pop_to_sample:
                print(f"[skip] {config_category} {gene} {region}: missing marine ref {marine}")
                continue

            fw_pops = AK_FW if region == "AK" else BC_FW
            ref_sample = pop_to_sample[marine]
            ref_aln = seqs[ref_sample]
            pos_map = aligned_to_ref_pos(ref_aln)

            variant_to_pops = defaultdict(list)

            for pop in sorted(fw_pops):
                if pop not in pop_to_sample:
                    continue

                q_sample = pop_to_sample[pop]
                q_aln = seqs[q_sample]
                variants = get_variants(ref_aln, q_aln)

                summary_rows.append({
                    "category": category,
                    "config_category": config_category,
                    "gene": gene,
                    "complex": complex_,
                    "region": region,
                    "marine_ref": marine,
                    "pop": pop,
                    "sample": q_sample,
                    "n_variants": len(variants),
                    "variants": ",".join(variants),
                    "alignment_file": aln_path,
                    "seq_type": seq_type,
                    "genetic_code": code_name
                })

                for v in variants:
                    variant_to_pops[v].append(pop)

            for variant, pops in variant_to_pops.items():
                parsed = parse_variant(variant)
                if parsed is None:
                    continue

                ref_aa, pos, alt_aa = parsed

                aln_col = None
                for i, p in pos_map.items():
                    if p == pos:
                        aln_col = i
                        break
                if aln_col is None:
                    continue

                site_aas = [
                    s[aln_col] for s in seqs.values()
                    if aln_col < len(s) and s[aln_col] not in BAD_AA
                ]
                if not site_aas:
                    continue

                counts = Counter(site_aas)
                major_aa, major_count = counts.most_common(1)[0]
                total = len(site_aas)

                major_freq = major_count / total
                ref_freq = counts.get(ref_aa, 0) / total
                alt_freq = counts.get(alt_aa, 0) / total
                entropy = shannon_entropy(site_aas)
                blosum = BLOSUM.get((ref_aa, alt_aa), None)

                unique_pops = sorted(set(pops))
                n_pops = len(unique_pops)

                all_rows.append({
                    "category": category,
                    "config_category": config_category,
                    "gene": gene,
                    "complex": complex_,
                    "region": region,
                    "marine_ref": marine,
                    "variant": variant,
                    "ref_aa": ref_aa,
                    "position": pos,
                    "alt_aa": alt_aa,
                    "n_pops": n_pops,
                    "populations": ",".join(unique_pops),
                    "alignment_column_1based": aln_col + 1,
                    "n_sequences_scored": total,
                    "major_aa": major_aa,
                    "major_freq": round(major_freq, 3),
                    "ref_aa_freq": round(ref_freq, 3),
                    "alt_aa_freq": round(alt_freq, 3),
                    "entropy": round(entropy, 3),
                    "BLOSUM62": blosum,
                    "BLOSUM_class": (
                        "radical" if blosum is not None and blosum <= -2 else
                        "nonconservative" if blosum is not None and blosum <= 0 else
                        "conservative" if blosum is not None else
                        "NA"
                    ),
                    "impact_class": impact_class(major_freq, blosum, n_pops),
                    "alignment_file": aln_path,
                    "seq_type": seq_type,
                    "genetic_code": code_name
                })

df = pd.DataFrame(all_rows)
summary = pd.DataFrame(summary_rows)

out1 = os.path.join(OUTDIR, "AA_substitutions_all_categories_BLOSUM_conservation.tsv")
out2 = os.path.join(OUTDIR, "AA_substitution_per_population_summary.tsv")
out3 = os.path.join(OUTDIR, "category_AA_summary.tsv")

if not df.empty:
    df = df.sort_values(
        ["category", "gene", "region", "n_pops", "major_freq", "BLOSUM62"],
        ascending=[True, True, True, False, False, True]
    )

df.to_csv(out1, sep="\t", index=False)
summary.to_csv(out2, sep="\t", index=False)

print("\nDONE")
print("Main output:", out1)
print("Per-pop summary:", out2)

if not df.empty:
    cat_summary = df.groupby("category").agg(
        n_unique_substitutions=("variant", "count"),
        n_genes=("gene", "nunique"),
        recurrent_substitutions=("n_pops", lambda x: (x >= 2).sum()),
        radical_substitutions=("BLOSUM62", lambda x: (x <= -2).sum()),
        nonconservative_substitutions=("BLOSUM62", lambda x: (x <= 0).sum())
    ).reset_index()

    cat_summary.to_csv(out3, sep="\t", index=False)

    print("\nCategory summary:")
    print(cat_summary.to_string(index=False))
    print("\nSaved category summary:", out3)
else:
    print("No substitutions found.")
PY





#sift for287
/work/cyu/run_sift4g_287_all_categories_teleostDB.py
#!/usr/bin/env python3
import os
import re
import subprocess
import pandas as pd
from pathlib import Path

# ============================================================
# Paths
# ============================================================
BASE = Path("/work/cyu/AA_substitution_all_categories")

INFILE = BASE / "AA_substitutions_all_categories_BLOSUM_conservation.tsv"

OUTBASE = BASE / "sift4g_287_all_categories_teleostUniRef90"
OUTBASE.mkdir(parents=True, exist_ok=True)

SIFT4G = "/home/cyu/.conda/envs/provean_env/bin/sift4g"

# 这里换成你的 teleost SIFT database
# 如果你的 teleost db 是 fasta，就放 fasta 路径
# 如果是已经 build 好的 SIFT4G db，也放那个路径
DB = "/work/cyu/provean_db/teleost_uniref90.fa"

BAD_AA = {"-", "X", "N", "*", "?"}

STD_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
MT_CODE = STD_CODE.copy()
MT_CODE.update({'TGA':'W', 'ATA':'M', 'AGA':'*', 'AGG':'*'})

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
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs

def write_fasta(path, name, seq):
    with open(path, "w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")

def is_dna_alignment(seqs):
    chars = set("".join(seqs.values()).upper())
    return chars.issubset(set("ACGTN-?"))

def translate_aligned_codon(seq, genetic_code):
    code = MT_CODE if genetic_code == "vertebrate_mito" else STD_CODE
    clean = seq.upper().replace("?", "N")
    aa = []
    for i in range(0, len(clean) - 2, 3):
        codon = clean[i:i+3]
        if codon == "---":
            aa.append("-")
        elif "-" in codon:
            aa.append("-")
        elif "N" in codon:
            aa.append("X")
        else:
            aa.append(code.get(codon, "X"))
    return "".join(aa)

def ungap(seq):
    return seq.replace("-", "")

def infer_ref_sample(seqs, marine_ref):
    for name in seqs:
        parts = re.split(r"[_\.\-\|]+", name)
        if marine_ref in parts:
            return name
    for name in seqs:
        if marine_ref in name:
            return name
    return None

def parse_sift_prediction(pred_file, variant):
    if not pred_file.exists():
        return None, "NA"

    with open(pred_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split()

            if len(fields) >= 3 and fields[0] == variant:
                pred = fields[1].lower()
                try:
                    score = float(fields[2])
                except Exception:
                    score = None
                return score, pred

    return None, "NA"

# ============================================================
# Main
# ============================================================
df = pd.read_csv(INFILE, sep="\t")

required = ["category", "gene", "region", "marine_ref", "variant", "alignment_file", "genetic_code"]
missing = [x for x in required if x not in df.columns]
if missing:
    raise ValueError(f"Missing columns in input: {missing}")

# 每个 gene-region 跑一次，而不是每个 substitution 跑一次
groups = df.groupby(
    ["category", "gene", "region", "marine_ref", "alignment_file", "genetic_code"],
    dropna=False
)

rows = []

for (category, gene, region, marine_ref, alignment_file, genetic_code), sub in groups:
    print(f"\n[prepare] {category} {gene} {region}")

    if not os.path.exists(alignment_file):
        print(f"[skip] missing alignment file: {alignment_file}")
        for _, r in sub.iterrows():
            rr = r.to_dict()
            rr["SIFT_score"] = None
            rr["SIFT_prediction"] = "missing_alignment"
            rr["SIFT_prediction_file"] = ""
            rows.append(rr)
        continue

    seqs_raw = read_fasta(alignment_file)

    if is_dna_alignment(seqs_raw):
        seqs = {k: translate_aligned_codon(v, genetic_code) for k, v in seqs_raw.items()}
    else:
        seqs = seqs_raw

    ref_sample = infer_ref_sample(seqs, marine_ref)

    if ref_sample is None:
        print(f"[skip] missing marine ref {marine_ref}: {gene} {region}")
        for _, r in sub.iterrows():
            rr = r.to_dict()
            rr["SIFT_score"] = None
            rr["SIFT_prediction"] = "missing_marine_ref"
            rr["SIFT_prediction_file"] = ""
            rows.append(rr)
        continue

    ref_protein = ungap(seqs[ref_sample])

    # SIFT query protein 不要有 stop
    ref_protein = ref_protein.replace("*", "")

    if len(ref_protein) < 20:
        print(f"[skip] too short ref protein: {gene} {region}")
        for _, r in sub.iterrows():
            rr = r.to_dict()
            rr["SIFT_score"] = None
            rr["SIFT_prediction"] = "too_short_ref"
            rr["SIFT_prediction_file"] = ""
            rows.append(rr)
        continue

    variants = sorted(set(sub["variant"].dropna().astype(str)))

    work = OUTBASE / f"{category}__{gene}__{region}"
    subst_dir = work / "subst"
    outdir = work / "out"

    subst_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    header = f"{category}_{gene}_{region}_{marine_ref}"
    ref_fa = work / f"{header}.faa"
    subst_file = subst_dir / f"{header}.subst"

    write_fasta(ref_fa, header, ref_protein)

    with open(subst_file, "w") as out:
        for v in variants:
            out.write(v + "\n")

    cmd = [
        SIFT4G,
        "-q", str(ref_fa),
        "-d", DB,
        "--subst", str(subst_dir),
        "--out", str(outdir),
        "-t", "4"
    ]

    print("[run]", " ".join(cmd))

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print(f"[warn] SIFT failed: {category} {gene} {region}")
        for _, r in sub.iterrows():
            rr = r.to_dict()
            rr["SIFT_score"] = None
            rr["SIFT_prediction"] = "FAILED"
            rr["SIFT_prediction_file"] = ""
            rows.append(rr)
        continue

    pred_file = outdir / f"{header}.SIFTprediction"

    for _, r in sub.iterrows():
        variant = str(r["variant"])
        score, pred = parse_sift_prediction(pred_file, variant)

        rr = r.to_dict()
        rr["SIFT_score"] = score
        rr["SIFT_prediction"] = pred
        rr["SIFT_prediction_file"] = str(pred_file)
        rows.append(rr)

out = pd.DataFrame(rows)

outfile = BASE / "AA_substitutions_287_all_categories_BLOSUM_conservation_SIFT4G_teleostUniRef90.tsv"
out.to_csv(outfile, sep="\t", index=False)

print("\nDONE")
print("Saved:", outfile)

summary = (
    out.groupby(["category", "SIFT_prediction"])
    .size()
    .reset_index(name="n")
    .sort_values(["category", "SIFT_prediction"])
)

summary_file = BASE / "SIFT4G_287_category_summary_teleostUniRef90.tsv"
summary.to_csv(summary_file, sep="\t", index=False)

print("Summary:", summary_file)
print(summary.to_string(index=False))


