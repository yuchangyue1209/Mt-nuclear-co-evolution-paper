#!/usr/bin/env python3

import subprocess
import pandas as pd
from pathlib import Path

BASE = Path("/work/cyu/provean_inputs_OXPHOS72")
DB = "/work/cyu/provean_db/teleost_uniref90.fa"
SIFT4G = "/home/cyu/.conda/envs/provean_env/bin/sift4g"

INFILE = BASE / "unique_AA_substitutions_by_population.tsv"
OUTBASE = BASE / "sift4g_all_OXPHOS72_teleostUniRef90"
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
outfile = BASE / "SIFT4G_candidate_substitution_results_teleostUniRef90.tsv"
out.to_csv(outfile, sep="\t", index=False)

print("\nDONE")
print("Saved:", outfile)
print(out.to_string(index=False))
