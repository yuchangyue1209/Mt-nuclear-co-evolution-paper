#colabfold protein prediction
mkdir -p /work/cyu/OXPHOS_structural_validation/{input,protein_fasta,colabfold,pymol,foldx,results}
cd /work/cyu/OXPHOS_structural_validation
nano input/candidate_substitutions.tsv
Gene	Variant	Region_pattern	n_pops	Conserv	BLOSUM62	SIFT4G	Functional_class
ndufs7	F19S	AK+BC_recurrent	AK5+BC9	0.963	-2	deleterious	parallel_deleterious
ndufa7	P42A	AK_BC_parallel	2	0.926	-1	deleterious	parallel_constrained
ndufb10	P61S	AK_BC_parallel	2	0.923	-1	tolerated	parallel_tolerated
ndufs6	A22T	AK_BC_parallel	2	0.926	0	tolerated	parallel_tolerated
atp5pd	P80L	BC_specific	1	0.963	-3	deleterious	constrained_radical
ndufs1	P544L	BC_specific	1	0.963	-3	deleterious	constrained_radical
ndufs7	C33G	BC_recurrent	9	0.92	-3	tolerated	recurrent_radical
atp5f1b	G19R	BC_specific	1	0.963	-2	deleterious	constrained_radical
cox6b1	F84S	AK_specific	1	0.963	-2	deleterious	constrained_radical
ndufs7	G8R	BC_recurrent	2	0.926	-2	deleterious	constrained_radical
cyc1	T44M	BC_specific	1	0.963	-1	deleterious	constrained_radical
ndufa9a	T266M	AK_specific	1	0.963	-1	deleterious	constrained_radical
uqcrc2b	P45S	AK_specific	1	0.963	-1	deleterious	constrained_radical
ndufa6	T48A	AK_specific	1	0.963	0	deleterious	constrained_radical
ndufs2	T192A	AK_specific	1	0.963	0	deleterious	constrained_radical
uqcrc2b	A16V	AK_specific	1	0.963	0	deleterious	constrained_radical
uqcrq	G44A	BC_recurrent	11	0.963	0	tolerated	highly_recurrent
ndufa4	E3G	AK_recurrent	9	0.926	-2	tolerated	recurrent_radical
atp5f1b	A38S/S38A	AK_vs_BC_reciprocal	8	0.682	1	deleterious	reciprocal_parallelism
ndufs8	Y37H	AK+BC	2	0.917	2	deleterious	parallel_deleterious



cat > make_priority_table.py <<'PY'
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
PY

python3 make_priority_table.py


#aa seq  v 1.6.1 
mkdir -p /work/cyu/OXPHOS_structural_validation/colabfold

cd /work/cyu/provean_inputs_OXPHOS72

for gene in ndufs7 atp5f1b ndufa7 ndufa4 uqcrq ndufs8 atp5pd ndufs1 cox6b1; do
    ls ${gene}_AK/*_ref_*.faa ${gene}_BC/*_ref_*.faa 2>/dev/null
done > /work/cyu/OXPHOS_structural_validation/colabfold/top_candidate_ref_fasta_paths.txt

cat $(cat /work/cyu/OXPHOS_structural_validation/colabfold/top_candidate_ref_fasta_paths.txt) \
> /work/cyu/OXPHOS_structural_validation/colabfold/top_candidate_ref_proteins_for_colabfold.faa


cd /work/cyu/OXPHOS_structural_validation/colabfold
mkdir -p output

colabfold_batch \
  --num-recycle 3 \
  --model-type alphafold2_ptm \
  top_candidate_ref_proteins_for_colabfold.faa \
  output