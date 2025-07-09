# dn/ds calculation
#Merge the files into one fasta
cat /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/*.fasta > /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/all_mt_overlap.fasta

13. Realignment
mafft --auto all_mt_overlap.fasta > aligned_mt_overlap.fasta
muscle -in all_mt_overlap.fasta -out mc_aligned_mt_overlap.fasta

14. IQtree 
iqtree -s aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO
iqtree -s mc_aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO


/home/cyu/snpEff/data/Gasterosteus_aculeatus_MT/genes.gff


#去掉参考
from pathlib import Path
ref = "MH205729.1"
for phy in Path(".").glob("*_aligned_codon.seq.phy"):
    out = phy.with_suffix(".noRef.phy")
    with phy.open() as f:
        lines = f.readlines()
    # head 行
    nseq, seqlen = lines[0].split()
    nseq = str(int(nseq) - 1)
    new = [f"{nseq} {seqlen}\n"]
    # 过滤参考
    new += [ln for ln in lines[1:] if not ln.startswith(ref)]
    out.write_text("".join(new))
    print(f"✅ {out.name}")


    from pathlib import Path
tmpl = """seqfile = {phy}
treefile = freshwater_binary.tree
outfile  = {out}

noisy = 3      verbose = 1     runmode = 0
seqtype = 1    CodonFreq = 2
model  = 0     NSsites = {nss}
icode  = 2
fix_kappa = 0  kappa = 2
fix_omega = 0  omega = 1
cleandata = 1
"""

for phy in Path(".").glob("*_aligned_codon.seq.noRef.phy"):   # ← 改这里
    gene = phy.stem.split("_")[0]          # ND4, COX1 ...
    for tag,nss in [("M1",1),("M2",2)]:    # M1a vs M2a
        ctl = f"{gene}_{tag}.ctl"
        Path(ctl).write_text(
            tmpl.format(phy=phy.name,
                        out=f"{gene}_{tag}.out",
                        nss=nss)
        )
        print(f"✔ 写 {ctl}")




/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus

from Bio import SeqIO
from Bio.Seq import Seq

# 所有13个基因的参考坐标（0-based for Python）
genes = [
    ("ND1", 2849, 3824, "+"),
    ("ND2", 4034, 5081, "+"),
    ("COX1", 5463, 7014, "+"),
    ("COX2", 7171, 7862, "+"),
    ("ATP8", 7937, 8105, "+"),
    ("ATP6", 8095, 8779, "+"),
    ("COX3", 8778, 9563, "+"),
    ("ND3", 9633, 9982, "+"),
    ("ND4L", 10051, 10348, "+"),
    ("ND4", 10341, 11722, "+"),
    ("ND5", 11938, 13777, "+"),
    ("ND6", 13773, 14295, "-"),
    ("CYTB", 14370, 15511, "+"),
]

input_fasta = "aligned_all_pops_plus_ref.fasta"
records = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))

# 遍历每个基因
for gene, start, end, strand in genes:
    dna_out = open(f"{gene}_aligned_clean.fasta", "w")
    protein_out = open(f"{gene}_aligned_protein.fasta", "w")
    print(f"\nProcessing {gene}...")

    for record in records.values():
        name = record.id
        subseq = record.seq[start:end]

        # 去除 gaps
        clean_seq = str(subseq).replace("-", "")

        # 如果负链则取反向互补
        if strand == "-":
            clean_seq = str(Seq(clean_seq).reverse_complement())

        # 截断为3的倍数以便正确翻译
        if len(clean_seq) % 3 != 0:
            clean_seq = clean_seq[:len(clean_seq) - (len(clean_seq) % 3)]

        # 翻译为氨基酸（线粒体密码表）
        protein_seq = Seq(clean_seq).translate(table=2)

        # 写入文件
        dna_out.write(f">{name}\n{clean_seq}\n")
        protein_out.write(f">{name}\n{protein_seq}\n")

        print(f"{name} | DNA: {len(clean_seq)} bp | AA: {len(protein_seq)} aa")

    dna_out.close()
    protein_out.close()
    print(f"{gene} done: output to {gene}_aligned_clean.fasta and {gene}_aligned_protein.fasta")




步骤 2：对 13 个基因进行蛋白质比对（MUSCLE）

#!/bin/bash

genes=(ND1 ND2 COX1 COX2 ATP8 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 CYTB)

for gene in "${genes[@]}"; do
    echo "Running MUSCLE on $gene..."
    muscle -align /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/${gene}_aligned_protein.fasta -output /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/protein/${gene}_prot_muscle.fasta
done

cd /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus


#!/bin/bash
PROTDIR=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/protein

for f in "${PROTDIR}"/*_prot_muscle.fasta; do
  echo "Processing $f …"
  awk '
    /^>MH205729\.1$/ { skip=1; next }   # 遇到参考头则开始跳过
    /^>/ && skip==1 { skip=0 }           # 下一个序列头出现时取消跳过（不跳过这行）
    skip==0                            # 仅在 skip==0 时才打印
  ' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done

#!/bin/bash

NUCDIR=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus

for gene in ND1 ND2 COX1 COX2 ATP8 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 CYTB; do
  f="${NUCDIR}/${gene}_aligned_clean.fasta"
  echo "Filtering reference out of $f …"
  awk '
    /^>MH205729\.1$/ { skip=1; next }  # 遇到参考头开始跳过
    /^>/ && skip==1 { skip=0 }          # 下一个头出现时恢复打印
    skip==0                             # 只有 skip=0 时才输出
  ' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done


#codon phy

NUCDIR=/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus
PROTDIR=${NUCDIR}/protein
CODONDIR=${NUCDIR}/codon

for gene in ND1 ND2 COX1 COX2 ATP8 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 CYTB; do
  echo "Back‐translating $gene…"
  pal2nal.pl \
    "${PROTDIR}/${gene}_prot_muscle.fasta" \
    "${NUCDIR}/${gene}_aligned_clean.fasta" \
    -output paml \
    -codontable 2 \
    > "${CODONDIR}/${gene}_codon.phy"
done


#!/usr/bin/env python3
# concat_prot.py
# Concatenate protein alignments into a PHYLIP file for RAxML

from Bio import AlignIO
import glob
import os

# 1) Set your protein alignment directory
PROTDIR = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/protein"

# Change working directory
os.chdir(PROTDIR)

# 2) Find all *_prot_muscle.fasta files
files = sorted(glob.glob(os.path.join(PROTDIR, "*_prot_muscle.fasta")))
if not files:
    raise FileNotFoundError(f"No files matching *_prot_muscle.fasta in {PROTDIR}")

# 3) Read IDs from the first alignment
aln0 = AlignIO.read(files[0], "fasta")
ids = [rec.id for rec in aln0]

# 4) Initialize dict to hold concatenated sequences
d = {seq_id: '' for seq_id in ids}

# 5) Loop through each alignment and append
for fpath in files:
    print(f"Processing {os.path.basename(fpath)}...")
    aln = AlignIO.read(fpath, "fasta")
    seqdict = {rec.id: str(rec.seq) for rec in aln}
    missing = set(ids) - set(seqdict.keys())
    if missing:
        raise ValueError(f"Missing IDs {missing} in {fpath}")
    for seq_id in ids:
        d[seq_id] += seqdict[seq_id]

# 6) Write PHYLIP format
out_file = os.path.join(PROTDIR, "concatenated_prot.phy")
aln_len = len(next(iter(d.values())))
with open(out_file, "w") as fh:
    fh.write(f"{len(ids)} {aln_len}\n")
    for seq_id in ids:
        name = seq_id[:10].ljust(10)
        fh.write(name + d[seq_id] + "\n")

print(f"Written concatenated alignment to {out_file}")





#构建topology


cd /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/protein

raxmlHPC \
  -s concatenated_prot.phy \
  -n mito_concat \
  -m PROTGAMMAWAG \
  -f a \
  -x 12345 \
  -p 12345 \
  -# 100


sed -e 's/\(SAY\):/\1#1:/g' \
    -e 's/\(RS\):/\1#1:/g' \
    mito27_labeled.tre > mito_labeled.tre


cd /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/codon

cat > codeml.ctl << 'EOF'
seqfile = placeholder.phy    * alignment file, will be replaced by gene-specific .phy
treefile = mito27_labeled.tre * your Newick tree with #1 (foreground) and #0 (background) labels
outfile = placeholder.out     * output file, will be replaced by gene-specific .out

noisy = 3
verbose = 1

runmode = 0     * 0: user tree (with branch labels)

seqtype = 1     * 1 = codons
CodonFreq = 2   * F3x4

clock = 0       * no molecular clock
aaDist = 0

model = 2       * 2 = branch model (foreground/background)
NSsites = 0     * no site variation

icode = 2       * mitochondrial code
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.1     * initial omega

cleandata = 1   * remove sites with gaps
EOF


sed -i 's|^treefile = .*|treefile = mito_labeled.tre|' codeml.ctl

cd /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/protein



#!/usr/bin/env python3
# fix_codon_phylip_custom.py
# 将 *_codon.phy（interleaved）转换为严格 sequential PHYLIP (*_codon_fixed.phy)

import re
import glob
import os

def main():
    # 1) 设定目录
    CODONDIR = os.getcwd()  # 确保你已经 cd 到 codon/ 目录
    # 2) 遍历所有 .phy
    for fn in glob.glob("*_codon.phy"):
        if fn.endswith("_fixed.phy"):
            continue
        print("Reformatting", fn)
        # 读文件
        with open(fn) as f:
            header = f.readline().strip()
            seqs = {}
            order = []
            current = None
            for line in f:
                s = line.strip()
                if not s:
                    continue
                # 如果整行都是大写字母或数字，且长度<=10，则当 name
                if re.fullmatch(r"[A-Z0-9]+", s) and len(s) <= 10:
                    current = s
                    if current not in seqs:
                        seqs[current] = ""
                        order.append(current)
                else:
                    # 剩下的都是序列片段
                    if current:
                        seqs[current] += s
        # 写出新的 sequential PHYLIP
        out_fn = fn.replace(".phy", "_fixed.phy")
        with open(out_fn, "w") as out:
            out.write(f"{len(order)} {len(next(iter(seqs.values())))}\n")
            for name in order:
                out.write(f"{name.ljust(10)} {seqs[name]}\n")
        print("  -> wrote", out_fn)

if __name__ == "__main__":
    main()


model 2 example

cp codeml.ctl COX1_branch.ctl

sed -i "s|^seqfile = .*|seqfile = COX1_codon_fixed.phy|"  COX1_branch.ctl
sed -i "s|^treefile = .*|treefile = mito27_labeled.tre|"  COX1_branch.ctl
sed -i "s|^outfile = .*|outfile = COX1_branch.out|"      COX1_branch.ctl
sed -i "s|^model = .*|model = 2|"                       COX1_branch.ctl

codeml COX1_branch.ctl


cp codeml.ctl COX1_m0.ctl
sed -i "s|^seqfile = .*|seqfile = COX1_codon_fixed.phy|" COX1_m0.ctl
sed -i "s|^treefile = .*|treefile = mito_labeled.tre|" COX1_m0.ctl
sed -i "s|^outfile = .*|outfile = ATP6_m0.out|" ATP8_m0.ctl
sed -i "s|^model = .*|model = 0|" ATP8_m0.ctl
codeml ATP8_m0.ctl

#complex
# concat_codon.py
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os, glob

codon_dir = "/work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/codon"
os.chdir(codon_dir)

complexes = {
    "ComplexI":  ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"],
    "ComplexIII": ["CYTB"],
    "ComplexIV": ["COX1", "COX2", "COX3"],
    "ComplexV":  ["ATP6", "ATP8"]
}

for cname, genes in complexes.items():
    concat = {}
    ids = None
    for gene in genes:
        aln = AlignIO.read(f"{gene}_codon_fixed.phy", "phylip")
        if ids is None:
            ids = [rec.id for rec in aln]
            concat = {rec.id: "" for rec in aln}
        for rec in aln:
            concat[rec.id] += str(rec.seq)
    # format PHYLIP sequential
    out_file = f"{cname}_codon_concat.phy"
    with open(out_file, "w") as f:
        f.write(f"{len(ids)} {len(next(iter(concat.values())))}\n")
        for i in ids:
            f.write(f"{i.ljust(10)} {concat[i]}\n")
    print(f"✅ Written {out_file}")


cp codeml.ctl ComplexI_branch.ctl
sed -i "s|^seqfile = .*|seqfile = ComplexI_codon_concat.phy|" ComplexI_branch.ctl
sed -i "s|^outfile = .*|outfile = ComplexI_branch.out|" ComplexI_branch.ctl
codeml ComplexI_branch.ctl


cp codeml.ctl ComplexI_m0.ctl
sed -i "s|^seqfile = .*|seqfile = ComplexI_codon_concat.phy|" ComplexI_m0.ctl
sed -i "s|^outfile = .*|outfile = ComplexI_m0.out|" ComplexI_m0.ctl
sed -i "s|^model = .*|model = 0|" ComplexI_m0.ctl
codeml ComplexI_m0.ctl

