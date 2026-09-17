import pandas as pd
import re

inp = "input/candidate_substitutions.tsv"
out = "results/candidate_priority_table.tsv"

df = pd.read_csv(inp, sep="\t")

def recurrence_score(x):
    x = str(x)
    nums = [int(n) for n in re.findall(r'\d+', x)]
    total = sum(nums) if nums else 1
    if total >= 8:
        return 3
    elif total >= 2:
        return 2
    else:
        return 1

def score_row(row):
    score = 0

    # recurrence
    score += recurrence_score(row["n_pops"])

    # conservation
    if float(row["Conserv"]) >= 0.95:
        score += 3
    elif float(row["Conserv"]) >= 0.90:
        score += 2
    else:
        score += 1

    # BLOSUM
    b = int(row["BLOSUM62"])
    if b <= -2:
        score += 3
    elif b <= 0:
        score += 2
    else:
        score += 1

    # SIFT
    if str(row["SIFT4G"]).lower() == "deleterious":
        score += 3
    else:
        score += 1

    # special pattern
    fc = str(row["Functional_class"]).lower()
    if "parallel" in fc or "recurrent" in fc or "reciprocal" in fc:
        score += 2

    return score

df["priority_score"] = df.apply(score_row, axis=1)

def tier(s):
    if s >= 13:
        return "Tier1"
    elif s >= 10:
        return "Tier2"
    else:
        return "Tier3"

df["priority_tier"] = df["priority_score"].apply(tier)

df = df.sort_values(["priority_tier", "priority_score"], ascending=[True, False])
df.to_csv(out, sep="\t", index=False)

print(df[["Gene","Variant","priority_score","priority_tier"]])
print(f"\n[write] {out}")
