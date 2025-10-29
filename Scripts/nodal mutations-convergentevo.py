cd /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13

# 1) 收集每基因 codon 对齐
ls */*.codon.fas | sort > codon_align.list

# 2) 用小脚本拼接并生成“按基因×密码子位点”的分区文件
cat > concat_codon_make_partitions.py <<'PY'
import os
from collections import OrderedDict

def read_fasta(p):
    d=OrderedDict(); h=None; seq=[]
    for line in open(p):
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if h is not None: d[h]="".join(seq)
            h=line[1:].split()[0]; seq=[]
        else:
            seq.append(line)
    if h is not None: d[h]="".join(seq)
    return d

flist=[x.strip() for x in open("codon_align.list") if x.strip()]
genes=sorted({os.path.basename(os.path.dirname(f)) for f in flist})
gene_seqs={}; gene_len={}
samples=set()
for f in flist:
    g=os.path.basename(os.path.dirname(f))
    fa=read_fasta(f)
    L=set(len(s) for s in fa.values()); assert len(L)==1, f"length mismatch in {f}"
    L=L.pop()
    gene_len[g]=L
    gene_seqs[g]=fa
    samples|=set(fa.keys())
samples=sorted(samples, key=str.upper)

# 拼接
with open("mt_concat.codon.fas","w") as out:
    for s in samples:
        parts=[]
        for g in genes:
            seq=gene_seqs[g].get(s,"-"*gene_len[g])
            parts.append(seq)
        out.write(f">{s}\n{''.join(parts)}\n")

# 分区：DNA + 按 gene×(pos1,pos2,pos3)
start=1
with open("mt_concat.codon.partitions.txt","w") as p:
    for g in genes:
        L=gene_len[g]
        end=start+L-1
        # 连续区间的1/2/3位点（IQ-TREE/RAxML 风格）
        p.write(f"DNA, {g}_pos1 = {start}-{end}\\3\n")
        p.write(f"DNA, {g}_pos2 = {start+1}-{end}\\3\n")
        p.write(f"DNA, {g}_pos3 = {start+2}-{end}\\3\n")
        start=end+1
PY

python3 concat_codon_make_partitions.py

# 3) 构树（分区+自动并模）
iqtree2 \
  -s mt_concat.codon.fas \
  -p mt_concat.codon.partitions.txt \
  -m MFP+MERGE \
  -B 1000 -alrt 1000 \
  -nt AUTO \
  -pre mt_concat_nt

# 作为固定拓扑
cp mt_concat_nt.treefile mt_concat.tree


mkdir -p asr_mt
for F in */*.pep.aln.faa; do
  G=$(basename "$(dirname "$F")")
  iqtree2 -s "$F" -st AA -m MFP -asr -te mt_concat.tree -nt AUTO -pre "asr_mt/${G}"
done

#lineage_map.tsv
pop	lineage_id	region	status
SAY	Marine	BC	marine
RS	Marine	AK	marine

### Recent colonized
PACH	RC_Bamfield	BC	recent
FRED	RC_Bamfield	BC	recent
SC	RC_Alaska	AK	recent
CH	RC_Alaska	AK	recent

### Alaska long-established freshwater
FG	AK_MatSu	AK	freshwater
LG	AK_MatSu	AK	freshwater
SR	AK_MatSu	AK	freshwater
WB	AK_MatSu	AK	freshwater
LB	AK_MatSu	AK	freshwater

SL	AK_Kenai	AK	freshwater
TL	AK_Kenai	AK	freshwater
WT	AK_Kenai	AK	freshwater
WK	AK_Kenai	AK	freshwater

### British Columbia long-established freshwater
GOS	BC_Campbell	BC	freshwater
ROB	BC_Campbell	BC	freshwater
BOOT	BC_Campbell	BC	freshwater
ECHO	BC_Campbell	BC	freshwater
LAW	BC_Campbell	BC	freshwater

SWA	BC_NorthIsl	BC	freshwater
JOE	BC_NorthIsl	BC	freshwater
BEA	BC_NorthIsl	BC	freshwater
THE	BC_NorthIsl	BC	freshwater
MUC	BC_NorthIsl	BC	freshwater
PYE	BC_NorthIsl	BC	freshwater
AMO	BC_NorthIsl	BC	freshwater




#!/usr/bin/env python3
"""
rnm_directional_simple.py
----------------------------------------
Detects directionally-defined RNMs using RS+SAY as marine ancestors.
Counts parallel derived substitutions in freshwater and recent populations.

Usage:
  python3 rnm_directional_simple.py --aa_dir . --out_prefix rnm_dirA
"""

import os, re, argparse
from collections import Counter, defaultdict

# ===== helper =====
def read_fasta(fp):
    seqs = {}
    name = None
    with open(fp) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(seqs[name])
                name = re.sub(r'^\d+_', '', line[1:].split()[0])  # 去掉前缀数字
                seqs[name] = []
            else:
                seqs[name].append(line)
        if name is not None:
            seqs[name] = "".join(seqs[name])
    return seqs

# ===== main =====
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aa_dir", required=True, help="directory with per-gene fasta alignments (*.faa or *.fas)")
    ap.add_argument("--out_prefix", default="rnm_dirA")
    args = ap.parse_args()

    marine = ["RS", "SAY"]
    recent = ["PACH", "FRED", "SC", "CH"]
    # freshwater = everything else (inferred later)

    directional_hits = []

    for gdir in sorted(os.listdir(args.aa_dir)):
        gpath = os.path.join(args.aa_dir, gdir)
        if not os.path.isdir(gpath): continue
        fasta_files = [f for f in os.listdir(gpath) if f.endswith((".faa",".fas",".fasta"))]
        if not fasta_files: continue
        fpath = os.path.join(gpath, fasta_files[0])
        seqs = read_fasta(fpath)

        all_pops = list(seqs.keys())
        fw_pops = [p for p in all_pops if p not in marine + recent]
        L = len(next(iter(seqs.values())))

        # define ancestral base (RS+SAY consensus)
        if not all(m in seqs for m in marine):
            print(f"[WARN] Missing RS/SAY for {gdir}, skip.")
            continue
        anc_rs, anc_say = seqs["RS"], seqs["SAY"]
        ancestral = [anc_rs[i] if anc_rs[i] == anc_say[i] else "N" for i in range(L)]

        for pos in range(L):
            anc = ancestral[pos]
            if anc in "-N?": continue

            # look at derived in freshwater + recent
            derived_bases = {}
            for pop, seq in seqs.items():
                if pop in marine: continue
                if seq[pos] not in "-N?" and seq[pos] != anc:
                    derived_bases.setdefault(seq[pos], []).append(pop)

            for base, pops in derived_bases.items():
                if len(pops) >= 2:  # 平行突变
                    # classify habitat
                    habitat = "Recent" if any(p in recent for p in pops) else "Freshwater"
                    directional_hits.append((gdir, pos+1, anc, base, habitat, len(pops), ",".join(sorted(pops))))

        print(f"[{gdir}] directional RNM = {sum(1 for _ in directional_hits if _[0]==gdir)}")

    out_tsv = f"{args.out_prefix}.directional.tsv"
    with open(out_tsv, "w") as o:
        o.write("gene\tposition\tancestral\tderived\thabitat\tcount\tpops\n")
        for g,pos,anc,der,hab,n,pops in directional_hits:
            o.write(f"{g}\t{pos}\t{anc}\t{der}\t{hab}\t{n}\t{pops}\n")

    print(f"\n[OK] Wrote {out_tsv}")

if __name__ == "__main__":
    main()



python3 rnm_directional_simple.py --aa_dir . --out_prefix rnm_dirA




#haplotype newwork 根据mj 距离来看
 /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13/
python rnm_directional.py \
  --fasta /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13/mt_concat.codon.clean.fas \
  --lineage_map /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13/lineage_map.tsv \
  --sample_map  /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13/mt_concat.sample_map.tsv \
  --layout      /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13/mt_concat.layout.tsv \
  --out_prefix  /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13/rnm_dirA.directional


nano /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus/_mt_gene_align_13/

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Directional RNM/RCNM calls from SAY (marine ancestor) toward Freshwater/Recent.

Inputs
------
--fasta      : concatenated codon alignment FASTA (ACGT/U; gaps/others allowed, will become 'N')
--lineage_map: TSV with at least columns: pop, lineage_id, region, status
--sample_map : TSV with columns: sample, pop   (can be header-only; script will fallback to names)
--layout     : TSV with columns: gene,start,end,frame (1-based, closed; frame in {1,2,3})

Outputs
-------
<out_prefix>_nt.tsv : nucleotide-level events aggregated by gene:pos:ref>alt:habitat
<out_prefix>_aa.tsv : amino-acid-level events aggregated by gene:AApos:AAref>AAalt:habitat

Notes
-----
- Root is forced to the haplotype containing SAY (if present). If SAY not present, uses max-degree node.
- A directed edge parent->child is kept only if parent habitat == 'marine' and child in {'freshwater','recent'}.
- Haplotypes are defined as identical sequences in the concatenated FASTA (after cleaning to ACGT/N).
- Translation uses vertebrate mitochondrial code (NCBI table 2).
"""

import argparse, sys, re
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
import networkx as nx

ACGT = set("ACGT")
MT_TABLE = CodonTable.unambiguous_dna_by_id[2]  # vertebrate mitochondrial code

# --------------------------- Utilities ---------------------------

def clean_seq(s: str) -> str:
    s = s.upper().replace("U", "T")
    return "".join(ch if ch in ACGT else "N" for ch in s)

def seq_to_np(s: str) -> np.ndarray:
    return np.frombuffer(s.encode("ascii"), dtype=np.uint8)

def hamming_ignore_N(a: np.ndarray, b: np.ndarray) -> int:
    ok = (a != ord('N')) & (b != ord('N'))
    return int(np.sum(ok & (a != b)))

def diff_positions(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ok = (a != ord('N')) & (b != ord('N'))
    return np.where(ok & (a != b))[0]

def translate_codon(codon: str) -> str:
    codon = codon.replace("U","T").upper()
    if any(c not in "ACGT" for c in codon) or len(codon)!=3:
        return "?"
    if codon in MT_TABLE.stop_codons:
        return "*"
    return MT_TABLE.forward_table.get(codon, "?")

def parse_member_to_pop(header: str, num_to_pop: dict) -> str|None:
    """
    Try to parse a population name from sequence header.
    - If a leading number is found and present in num_to_pop (non-NA), return mapped pop.
    - Else try to extract a token like SAY/RS/FG/... directly.
    """
    h = header.strip()
    m = re.search(r'(^|\D)(\d{1,3})(?:[_\.\-]?)([A-Za-z][A-Za-z0-9]+)?', h)
    if m:
        num = m.group(2)
        if num in num_to_pop and pd.notna(num_to_pop[num]) and num_to_pop[num] != "":
            return num_to_pop[num].upper()
    m2 = re.search(r'([A-Za-z]{2,}[A-Za-z0-9]*)', h)
    return m2.group(1).upper() if m2 else None

# --------------------------- Loading helpers ---------------------------

def load_layout(layout_tsv: str) -> pd.DataFrame:
    df = pd.read_csv(layout_tsv, sep="\t", dtype={"gene":str,"start":int,"end":int,"frame":int})
    df["gene"] = df["gene"].str.upper()
    return df

def load_lineage_map(lineage_map_tsv: str) -> pd.DataFrame:
    lm = pd.read_csv(lineage_map_tsv, sep="\t", dtype=str).fillna("")
    lm["pop"] = lm["pop"].str.upper()
    lm["status"] = lm["status"].str.lower()
    return lm

def load_sample_map(sample_map_tsv: str) -> dict:
    """
    Returns dict num_to_pop; supports header-only (empty) files.
    """
    sm = pd.read_csv(sample_map_tsv, sep="\t", dtype=str)
    if sm.shape[0] == 0:
        return {}
    sm = sm.fillna("")
    assert set(["sample","pop"]).issubset(sm.columns), "sample_map must have columns: sample, pop"
    return {str(r.sample): (r.pop.upper() if isinstance(r.pop,str) else "") for r in sm.itertuples()}

# --------------------------- Haplotype construction ---------------------------

def build_haps(fasta_path: str):
    seqs = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        s = clean_seq(str(rec.seq))
        seqs.append((rec.id, s))
    if not seqs:
        raise RuntimeError(f"No sequences read from {fasta_path}")
    seq_to_members = defaultdict(list)
    for sid, s in seqs:
        seq_to_members[s].append(sid)
    hap_id = []
    hap_seq = []
    hap_members = {}
    for i, (s, members) in enumerate(seq_to_members.items(), start=1):
        name = f"H{i}"
        hap_id.append(name)
        hap_seq.append(s)
        hap_members[name] = members
    return hap_id, hap_seq, hap_members

def hap_graph_mst(hap_id, hap_seq):
    """
    Build MST on haplotypes using Hamming distance ignoring N.
    If multiple components, MST each and compose.
    """
    arr = [seq_to_np(s) for s in hap_seq]
    G = nx.Graph()
    for i, h in enumerate(hap_id):
        G.add_node(h, idx=i)
    for i in range(len(hap_id)):
        for j in range(i+1, len(hap_id)):
            d = hamming_ignore_N(arr[i], arr[j])
            if d == 0:
                continue
            G.add_edge(hap_id[i], hap_id[j], weight=d)
    if G.number_of_edges() == 0:
        return G, arr
    if not nx.is_connected(G):
        H = nx.Graph()
        for comp in nx.connected_components(G):
            sub = G.subgraph(comp).copy()
            if sub.number_of_edges() > 0:
                mst = nx.minimum_spanning_tree(sub, weight="weight")
                H = nx.compose(H, mst)
            else:
                H = nx.compose(H, sub)
        G = H
    else:
        G = nx.minimum_spanning_tree(G, weight="weight")
    return G, arr

# --------------------------- Mapping / aggregation ---------------------------

def pop_status(pop: str, lineage_map_df: pd.DataFrame) -> str|None:
    row = lineage_map_df.loc[lineage_map_df["pop"] == pop]
    if len(row) == 0:
        return None
    return row["status"].iloc[0]

def map_concat_to_gene(layout_df: pd.DataFrame, pos_concat: int):
    """
    Given 1-based concat position, return (gene, pos_in_gene, codon_pos(1..3), aa_pos, gene_start, frame)
    """
    w = (layout_df["start"] <= pos_concat) & (pos_concat <= layout_df["end"])
    if not w.any():
        return None
    r = layout_df[w].iloc[0]
    gene = r["gene"]
    pos_in_gene = int(pos_concat - r["start"] + 1)
    frame = int(r["frame"])
    codon_pos = ((pos_in_gene - frame) % 3) + 1
    aa_pos = (pos_in_gene - frame) // 3 + 1
    return gene, pos_in_gene, codon_pos, aa_pos, int(r["start"]), frame

def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Input columns required: gene, position, ancestral, derived, habitat, pops (comma-joined)
    Output: gene, position, ancestral, derived, habitat, count, pops
    """
    need = {"gene","position","ancestral","derived","habitat","pops"}
    if df is None or len(df)==0 or not need.issubset(df.columns):
        return pd.DataFrame(columns=["gene","position","ancestral","derived","habitat","count","pops"])

    keys = ["gene","position","ancestral","derived","habitat"]

    rows = []
    for key_vals, sub in df.groupby(keys, dropna=False):
        pop_set = set()
        for x in sub["pops"].dropna():
            if isinstance(x, str):
                for p in x.split(","):
                    p = p.strip()
                    if p:
                        pop_set.add(p)
        lst = sorted(pop_set)
        row = {k:v for k, v in zip(keys, key_vals)}
        row["count"] = len(lst)
        row["pops"]  = ",".join(lst)
        rows.append(row)

    if not rows:
        return pd.DataFrame(columns=["gene","position","ancestral","derived","habitat","count","pops"])

    out = pd.DataFrame(rows)
    out = out.sort_values(["gene","habitat","count","position"],
                          ascending=[True, True, False, True],
                          kind="mergesort")
    return out[["gene","position","ancestral","derived","habitat","count","pops"]]

# --------------------------- Rooting helpers ---------------------------

def choose_root_SAY(hap_members, num_to_pop) -> str|None:
    """
    Choose SAY-containing hap as root if present; else None.
    """
    for h, members in hap_members.items():
        pops = {parse_member_to_pop(m, num_to_pop) for m in members}
        if "SAY" in pops:
            return h
    return None

# --------------------------- Main ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Directional RNMs from SAY to FW/Recent")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--lineage_map", required=True)
    ap.add_argument("--sample_map", required=True)
    ap.add_argument("--layout", required=True)
    ap.add_argument("--out_prefix", default="rnm_dirA.directional")
    args = ap.parse_args()

    layout = load_layout(args.layout)
    lineage = load_lineage_map(args.lineage_map)
    num_to_pop = load_sample_map(args.sample_map)

    hap_id, hap_seq, hap_members = build_haps(args.fasta)
    G, hap_arr = hap_graph_mst(hap_id, hap_seq)

    if G.number_of_nodes() == 0:
        print("No haplotypes found.", file=sys.stderr)
        pd.DataFrame(columns=["gene","position","ancestral","derived","habitat","count","pops"]).to_csv(
            f"{args.out_prefix}_nt.tsv", sep="\t", index=False)
        pd.DataFrame(columns=["gene","position","ancestral","derived","habitat","count","pops"]).to_csv(
            f"{args.out_prefix}_aa.tsv", sep="\t", index=False)
        print(f"Done:\n  {args.out_prefix}_nt.tsv\n  {args.out_prefix}_aa.tsv")
        return

    root = choose_root_SAY(hap_members, num_to_pop)
    if root is None:
        if G.number_of_nodes() == 1:
            root = list(G.nodes())[0]
        else:
            root = max(G.degree, key=lambda x: x[1])[0]
        print("⚠️ SAY hap not found; using highest-degree hap as root:", root, file=sys.stderr)
    else:
        print("Root hap (SAY):", root, file=sys.stderr)

    parent = {root: None}
    order = [root]
    for u, v in nx.bfs_edges(G, root):
        parent[v] = u
        order.append(v)

    status_order = {"marine":0, "recent":1, "freshwater":2}
    hap2pops = {}
    hap2hab  = {}
    for h, members in hap_members.items():
        pops = []
        for m in members:
            p = parse_member_to_pop(m, num_to_pop)
            if p:
                pops.append(p)
        pops = sorted(set(pops))
        hap2pops[h] = pops
        st = None
        for p in pops:
            s = pop_status(p, lineage)
            if s is None:
                continue
            if st is None or status_order.get(s, 99) < status_order.get(st, 99):
                st = s
        hap2hab[h] = st

    records_nt = []
    records_aa = []
    id_index = {h:i for i,h in enumerate(hap_id)}

    for h in order:
        p = parent.get(h, None)
        if p is None:
            continue
        if not (hap2hab.get(p) == "marine" and hap2hab.get(h) in ("freshwater","recent")):
            continue

        a = hap_arr[id_index[p]]
        b = hap_arr[id_index[h]]
        diffs = diff_positions(a, b)  # 0-based concat
        for pos0 in diffs:
            pos1 = int(pos0 + 1)  # 1-based
            gi = map_concat_to_gene(layout, pos1)
            if gi is None:
                continue
            gene, pos_in_gene, codon_pos, aa_pos, gene_start, frame = gi

            anc_nt = chr(int(a[pos0])); der_nt = chr(int(b[pos0]))
            records_nt.append(dict(
                gene=gene,
                position=pos_in_gene,
                ancestral=anc_nt,
                derived=der_nt,
                habitat=("Freshwater" if hap2hab[h]=="freshwater" else "Recent"),
                pops=",".join(hap2pops[h])
            ))

            codon_start = gene_start + (aa_pos-1)*3 + (frame-1)  # concat 1-based
            codon_idx = [codon_start-1, codon_start, codon_start+1]  # 0-based
            try:
                fb = "".join(chr(int(a[i])) for i in codon_idx)
                tb = "".join(chr(int(b[i])) for i in codon_idx)
            except IndexError:
                continue
            if all(c in "ACGT" for c in fb) and all(c in "ACGT" for c in tb):
                AA_anc = translate_codon(fb)
                AA_der = translate_codon(tb)
                records_aa.append(dict(
                    gene=gene,
                    position=aa_pos,
                    ancestral=AA_anc,
                    derived=AA_der,
                    habitat=("Freshwater" if hap2hab[h]=="freshwater" else "Recent"),
                    pops=",".join(hap2pops[h])
                ))

    # 可选调试输出：
    # print("NT records:", len(records_nt), file=sys.stderr)
    # print("AA records:", len(records_aa), file=sys.stderr)

    out_nt = aggregate(pd.DataFrame.from_records(records_nt))
    out_aa = aggregate(pd.DataFrame.from_records(records_aa))

    out_nt.to_csv(f"{args.out_prefix}_nt.tsv", sep="\t", index=False)
    out_aa.to_csv(f"{args.out_prefix}_aa.tsv", sep="\t", index=False)
    print(f"Done:\n  {args.out_prefix}_nt.tsv\n  {args.out_prefix}_aa.tsv")

if __name__ == "__main__":
    main()





#37gene mt tree cinsensus
cd /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus

# 1) 去掉 D-loop（保留 1..15652）
for f in *_consensus.fasta; do
  base="${f%.fasta}"
  awk -v start=1 -v end=15652 '
    BEGIN{seq=""; header=""}
    /^>/ {header=$0; next}
    {seq = seq $0}
    END{
      print header
      s = substr(seq, start, end-start+1)
      for (i=1; i<=length(s); i+=60) print substr(s, i, 60)
    }' "$f" > "${base}.noDloop.fasta"
done

# 2) 检查每个样本是否都变成 15652 bp
for f in *_consensus.noDloop.fasta; do
  len=$(awk '/^>/ {next} {L+=length($0)} END {print L}' "$f")
  printf "%-35s %6d bp\n" "$f" "$len"
done

# 3) 合并为一个多序列（未对齐）FASTA
cat *_consensus.noDloop.fasta > all_mt_noDloop.fasta
# 4) MAFFT 自动对齐
mafft --auto all_mt_noDloop.fasta > aligned_mt_noDloop.fasta

# 5) IQ-TREE（SAY 外群；你给的是 iqtree 一代，这里用它）
iqtree -s aligned_mt_noDloop.fasta -m GTR+G -bb 1000 -alrt 1000 -o SAY -nt AUTO

# （如果你有另外一个“mc_aligned_mt_overlap.fasta”也要跑一次）
# iqtree -s mc_aligned_mt_overlap.fasta -m GTR+G -bb 1000 -alrt 1000 -o SAY -nt AUTO



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
call_nodal_rnm_scan_tree.py
Scan a rooted tree for Nodal & RNM directly from tip alignment (no .state needed).

Usage:
  python call_nodal_rnm_scan_tree.py \
    --aln   aligned_mt_noDloop.fasta \
    --tree  aligned_mt_noDloop.fasta.treefile \
    --outgroup SAY \
    --tau 0.80 \
    --mmin 3 \
    --group-by site+der \
    --independence subset \
    --prefix mt_noDloop_scan
"""

import argparse
from collections import defaultdict
from typing import List, Dict, Set, Tuple
import pandas as pd
from Bio import Phylo, SeqIO

VALID = set("ACGT")

def read_alignment(path:str)->Dict[str,str]:
    seqs = {}
    L = None
    for rec in SeqIO.parse(path, "fasta"):
        s = str(rec.seq).upper()
        if L is None: L = len(s)
        elif len(s)!=L: raise ValueError(f"Alignment not equal length at {rec.id}")
        seqs[rec.id]=s
    if not seqs: raise ValueError("No sequences in alignment")
    return seqs

def root_tree(tree_path:str, outgroup:str):
    T = Phylo.read(tree_path, "newick")
    tips = {t.name for t in T.get_terminals()}
    if outgroup not in tips:
        raise ValueError(f"Outgroup {outgroup} not in tree tips")
    T.root_with_outgroup(outgroup)
    return T

def descendants(clade)->List[str]:
    return [t.name for t in clade.get_terminals()]

def bases_at(seqs:Dict[str,str], tips:List[str], site:int)->List[str]:
    idx = site-1
    out=[]
    for t in tips:
        s = seqs.get(t, "")
        if idx<0 or idx>=len(s): out.append("N")
        else: out.append(s[idx])
    return out

def prop_of(bases:List[str], base:str)->Tuple[float,int]:
    ok=[x for x in bases if x in VALID]
    n=len(ok)
    if n==0: return (float("nan"),0)
    return (sum(1 for x in ok if x==base)/n, n)

def independent(sets:List[Set[str]], mode:str="subset")->bool:
    if len(sets)<2: return False
    for i in range(len(sets)-1):
        for j in range(i+1,len(sets)):
            a,b=sets[i],sets[j]
            if mode=="subset":
                if a.issubset(b) or b.issubset(a): return False
            else: # no-overlap
                if len(a & b)>0: return False
    return True

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--aln", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--tau", type=float, default=0.80)
    ap.add_argument("--mmin", type=int, default=3)
    ap.add_argument("--group-by", choices=["site","site+der"], default="site+der")
    ap.add_argument("--independence", choices=["subset","no-overlap"], default="subset")
    ap.add_argument("--prefix", default="nodal_scan")
    args=ap.parse_args()

    seqs = read_alignment(args.aln)
    L = len(next(iter(seqs.values())))
    T = root_tree(args.tree, args.outgroup)

    # Precompute parent -> descendants
    parent_desc = {}
    for cl in T.find_clades(order="level"):
        parent_desc[id(cl)] = set(descendants(cl))

    # Collect all internal clades (exclude tips and root’s parent None)
    internal = [cl for cl in T.find_clades(order="preorder") if not cl.is_terminal()]

    records=[]
    # Scan each internal clade as a candidate Nodal node
    for cl in internal:
        child_tips = list(parent_desc[id(cl)])
        # parent clade = actual parent in the rooted tree
        parent = None
        for p in T.find_clades(order="preorder"):
            if cl in getattr(p,"clades",[]):
                parent = p; break
        parent_tips = list(parent_desc[id(parent)]) if parent is not None else list(parent_desc[id(cl)])

        for site in range(1, L+1):
            child_bases = bases_at(seqs, child_tips, site)
            parent_bases = bases_at(seqs, parent_tips, site)
            # try each base as candidate derived state
            for b in "ACGT":
                prop_c, n_c = prop_of(child_bases, b)
                if not (n_c>=args.mmin and prop_c>=args.tau): 
                    continue
                prop_p, n_p = prop_of(parent_bases, b)
                # require parent clade not already high for b (prevents trivial propagation)
                if n_p>0 and prop_p>=args.tau:
                    continue
                records.append({
                    "node_id": id(cl),
                    "site": site,
                    "der": b,
                    "prop_child": prop_c,
                    "n_child": n_c,
                    "prop_parent": prop_p,
                    "n_parent": n_p,
                    "descendants": ",".join(child_tips)
                })

    nodal_df = pd.DataFrame.from_records(records)
    nodal_df.to_csv(f"{args.prefix}.nodal_calls.tsv", sep="\t", index=False)

    # RNM grouping
    if nodal_df.empty:
        pd.DataFrame(columns=["site","der_key","n_events","events_descendants"]).to_csv(
            f"{args.prefix}.rnm_calls.tsv", sep="\t", index=False)
        print("Done. Nodal events: 0 ; RNM groups: 0")
        return

    if args.group_by=="site+der":
        nodal_df["grp"] = list(zip(nodal_df["site"].astype(int), nodal_df["der"].astype(str)))
    else:
        nodal_df["grp"] = nodal_df["site"].astype(int)

    rnm_rows=[]
    for grp, g in nodal_df.groupby("grp"):
        sets=[set(x.split(",")) for x in g["descendants"].tolist()]
        if len(sets)>=2 and independent(sets, args.independence):
            rnm_rows.append({
                "site": int(g.iloc[0]["site"]),
                "der_key": None if args.group_by=="site" else str(g.iloc[0]["der"]),
                "n_events": len(sets),
                "events_descendants": "|".join(g["descendants"].tolist())
            })
    rnm_df = pd.DataFrame.from_records(rnm_rows)
    if rnm_df.empty:
        pd.DataFrame(columns=["site","der_key","n_events","events_descendants"]).to_csv(
            f"{args.prefix}.rnm_calls.tsv", sep="\t", index=False)
    else:
        rnm_df.to_csv(f"{args.prefix}.rnm_calls.tsv", sep="\t", index=False)

    print(f"Done. Nodal events: {len(nodal_df)} ; RNM groups: {0 if rnm_df is None else len(rnm_rows)}")

if __name__=="__main__":
    main()
cd /work/cyu/poolseq/PPalign_output/overlap.vcf/consensus
python call_nodal_rnm_scan_tree.py \
  --aln    aligned_mt_noDloop.fasta \
  --tree   aligned_mt_noDloop.fasta.treefile \
  --outgroup SAY \
  --tau 0.9 \
  --mmin 2 \
  --group-by site+der \
  --independence subset \
  --prefix mt_noDloop



python annotate_mt_sites.py \
  mt_noDloop.nodal_calls.tsv  mt_partitions.tsv  mt_noDloop.nodal_calls.ann.tsv
python annotate_mt_sites.py \
  mt_noDloop.rnm_calls.tsv    mt_partitions.tsv  mt_noDloop.rnm_calls.ann.tsv


python label_syn_gapaware.py \
  aligned_mt_noDloop.fasta  mt_partitions.tsv \
  mt_noDloop.rnm_calls.ann.tsv  mt_noDloop.rnm_calls.syn_gap.tsv  SAY

