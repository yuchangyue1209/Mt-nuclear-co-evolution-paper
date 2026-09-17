#!/usr/bin/env python3
"""Prepare corrected mtPCG codon alignments and run codeml M0 on the 27-tip tree."""
import re, shutil, subprocess
from pathlib import Path

GENES = ["ATP6","ATP8","COX1","COX2","COX3","CYTB","ND1","ND2","ND3","ND4","ND4L","ND5","ND6"]
ALIGN_DIR = Path("/mnt/spareHD_2/mt_gene_tree/mt_consensus_aln")
TREE = Path("/work/cyu/Kuster2026_genomewide_reanalysis/genomewide_erc/00_inputs/erc_master_chr21_27_unrooted_paml.tree")
OUT = Path("/mnt/spareHD_2/genomewide_codeml_kuster/10_erc_mt13_codeml/results")
MT_STOPS = {"TAA","TAG","AGA","AGG"}

def fasta(path):
    out=[]; name=None; seq=[]
    for raw in path.open():
        line=raw.strip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: out.append((name,"".join(seq).upper()))
            name=line[1:].split()[0]; seq=[]
        else: seq.append(line)
    if name is not None: out.append((name,"".join(seq).upper()))
    return out

def tips(path):
    return set(re.findall(r"(?<=[(,])([A-Za-z0-9_.-]+)(?=[:),;])",path.read_text()))

def clean(gene, records):
    lengths={len(s) for _,s in records}
    if len(lengths)!=1: raise ValueError(f"{gene}: unequal lengths {lengths}")
    remainder=next(iter(lengths))%3
    if remainder:
        print(f"[trim] {gene}: removing {remainder} incomplete terminal nucleotide(s)")
        records=[(n,s[:-remainder]) for n,s in records]
    cleaned=[]
    for name,seq in records:
        # Historical sync counts A:T:C:G were mislabeled A:C:G:T.
        seq=seq.replace("U","T").translate(str.maketrans("ACGT","ATCG"))
        if gene=="ND6": seq=seq.translate(str.maketrans("ACGT","TGCA"))[::-1]
        if seq[-3:] in MT_STOPS: seq=seq[:-3]+"---"
        bad=[i//3+1 for i in range(0,len(seq)-3,3) if seq[i:i+3] in MT_STOPS]
        if bad: raise ValueError(f"{gene}/{name}: internal stops {bad[:10]}")
        cleaned.append((name,seq))
    return cleaned

def write_ctl(path):
    path.write_text("""seqfile = alignment.phy\ntreefile = tree.nwk\noutfile = mlc\nnoisy = 0\nverbose = 0\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\naaDist = 0\nmodel = 0\nNSsites = 0\nicode = 1\nfix_kappa = 0\nkappa = 2\nfix_omega = 0\nomega = 0.2\nfix_alpha = 1\nalpha = 0\nMalpha = 0\nncatG = 8\ngetSE = 0\nRateAncestor = 0\nSmall_Diff = 5e-7\ncleandata = 1\nmethod = 0\n""")

def main():
    if not shutil.which("codeml"): raise SystemExit("codeml unavailable")
    expected=tips(TREE)
    if len(expected)!=27: raise SystemExit(f"Tree has {len(expected)} tips")
    OUT.mkdir(parents=True,exist_ok=True)
    for gene in GENES:
        rec=fasta(ALIGN_DIR/f"{gene}.mt.aln.fasta")
        if len(rec)!=27 or {n for n,_ in rec}!=expected: raise ValueError(f"{gene}: 27-tip mismatch")
        rec=clean(gene,rec); d=OUT/gene; d.mkdir(parents=True,exist_ok=True)
        with (d/"alignment.phy").open("w") as h:
            h.write(f"27 {len(rec[0][1])}\n"); [h.write(f"{n}  {s}\n") for n,s in rec]
        shutil.copy2(TREE,d/"tree.nwk"); write_ctl(d/"codeml.ctl")
        print(f"[run] {gene}: samples=27, nt={len(rec[0][1])}",flush=True)
        with (d/"codeml.stdout.log").open("w") as log:
            p=subprocess.run(["codeml","codeml.ctl"],cwd=d,stdout=log,stderr=subprocess.STDOUT)
        mlc=d/"mlc"
        if p.returncode or not mlc.is_file() or "lnL(ntime" not in mlc.read_text(errors="replace"):
            raise RuntimeError(f"{gene}: codeml failed; inspect {d/'codeml.stdout.log'}")
        print(f"[OK] {gene}")
    print(f"[DONE] All 13 mtPCGs completed: {OUT}")
if __name__=="__main__": main()
