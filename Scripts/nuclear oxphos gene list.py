nuclear oxphos gene list
# 必备输入
GFF=/work/cyu/stickleback_v5.gff3.gz
FA=/work/cyu/stickleback_nuclear_only.fa
GENES=/work/cyu/oxphos_from_ref_no_biomart/human_oxphos_symbols.txt

# 映射/配置
MAP=/work/cyu/oxphos_from_ref_no_biomart/06_igv/ortholog_map.subunits.try1.tsv   # 你逐步补充的映射表（含别名/转录本）
OUTDIR=/work/cyu/oxphos_from_ref_no_biomart/06_igv

# 脚本（建议放 repo 的 scripts/ 下）
PY_BUILD=/work/cyu/oxphos_from_ref_no_biomart/build_oxphos_bed.py
PY_ANNOT=/work/cyu/oxphos_from_ref_no_biomart/annotate_span_bed.py
PY_DIFF=/work/cyu/oxphos_from_ref_no_biomart/diff_oxphos_sets_mapped.py

# 索引
[ -f "${FA}.fai" ] || samtools faidx "$FA"


python "$PY_BUILD" \
  --gff "$GFF" \
  --genes "$GENES" \
  --out-prefix "$OUTDIR/nuOXPHOS_v5.v3" \
  --feature cds \
  --paralog_policy all \
  --faidx "${FA}.fai" \
  --ortholog-map "$MAP"

bedtools sort -faidx "${FA}.fai" \
  -i "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.bed" \
  > "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.sorted.bed"

awk 'BEGIN{FS=OFS="\t"}{g=$4; if(!(k[g]++)){chr[g]=$1;s[g]=$2;e[g]=$3;st[g]=$6}else{if($2<s[g])s[g]=$2;if($3>e[g])e[g]=$3}} END{for(g in chr) print chr[g],s[g],e[g],g,".",st[g]}' \
  "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.sorted.bed" \
| bedtools sort -faidx "${FA}.fai" \
> "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.SPAN.bed"

python "$PY_ANNOT" \
  --span  "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.SPAN.bed" \
  --parts "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.sorted.bed" \
  --out   "$OUTDIR/nuOXPHOS_v5.v3.annot.bed"
#qc
python "$PY_DIFF" \
  --span-bed "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.SPAN.bed" \
  --ortholog-map "$MAP" \
  --out-prefix "$OUTDIR/nuOXPHOS_vs_human.v3map.nuclear_only" \
  --ignore-mt --keep-ab

# 汇总 & 缺失清单（subunits 目前 67/89；缺 22）
cat "$OUTDIR/nuOXPHOS_vs_human.v3map.nuclear_only.summary.tsv"
cat "$OUTDIR/nuOXPHOS_vs_human.v3map.nuclear_only.OXPHOS_subunits.missing.txt"


# 产一个 exon 版
python "$PY_BUILD" --gff "$GFF" --genes "$GENES" \
  --out-prefix "$OUTDIR/nuOXPHOS_v5.v3_exon" \
  --feature exon --paralog_policy all --faidx "${FA}.fai" --ortholog-map "$MAP"

# 找 zero-CDS 基因并用 exon 替换它们的坐标
awk -F'\t' 'NR>1 && $3==0 && $4!~/NOT_FOUND|SKIPPED/' \
  "$OUTDIR/nuOXPHOS_v5.v3.mapping_report.tsv" | cut -f2 | sort -u > "$OUTDIR/zero_CDS.genes.txt"

awk 'BEGIN{FS=OFS="\t"} NR==FNR{u[toupper($0)]=1;next} !(toupper($4) in u)' \
  "$OUTDIR/zero_CDS.genes.txt" "$OUTDIR/nuOXPHOS_v5.v3.cds_by_gene.bed" > "$OUTDIR/_no_zeroCDS.bed"

awk 'BEGIN{FS=OFS="\t"} NR==FNR{u[toupper($0)]=1;next} (toupper($4) in u)' \
  "$OUTDIR/zero_CDS.genes.txt" "$OUTDIR/nuOXPHOS_v5.v3_exon.cds_by_gene.bed" > "$OUTDIR/_only_zeroCDS_from_exon.bed"

cat "$OUTDIR/_no_zeroCDS.bed" "$OUTDIR/_only_zeroCDS_from_exon.bed" > "$OUTDIR/nuOXPHOS_v5.v3.mix_cds_exon.by_gene.bed"

# SPAN & 标注
awk 'BEGIN{FS=OFS="\t"}{g=$4; if(!(k[g]++)){chr[g]=$1;s[g]=$2;e[g]=$3;st[g]=$6}else{if($2<s[g])s[g]=$2;if($3>e[g])e[g]=$3}} END{for(g in chr) print chr[g],s[g],e[g],g,".",st[g]}' \
  "$OUTDIR/nuOXPHOS_v5.v3.mix_cds_exon.by_gene.bed" \
| bedtools sort -faidx "${FA}.fai" > "$OUTDIR/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.bed"

python "$PY_ANNOT" \
  --span  "$OUTDIR/nuOXPHOS_v5.v3.mix_cds_exon.SPAN.bed" \
  --parts "$OUTDIR/nuOXPHOS_v5.v3.mix_cds_exon.by_gene.bed" \
  --out   "$OUTDIR/nuOXPHOS_v5.v3.mix.annot.bed"



#build_oxphos_bed.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build OXPHOS BED from a human gene list + stickleback GFF/GTF.

Features
- Read .gff3/.gtf (plain or .gz)
- Liberal attribute parsing (gene_name/gene/Name/gene_id/ID/transcript_id/Parent)
- Optional MT-* mapping to fish names (COX1, ND1, ...)
- Paralog policy: all | prefer_a | prefer_b (teleost a/b copies)
- Optional ortholog map (human -> stickleback gene OR transcript name)
- Optional --faidx to restrict outputs to contigs present in a FASTA index (e.g., nuclear-only)

Outputs (with --out-prefix PREFIX):
  - PREFIX.cds_parts.bed          # raw parts (CDS/exon/gene)
  - PREFIX.cds_by_gene.bed        # merged by *gene* (multi-block merged)
  - PREFIX.cds_by_tx.bed          # merged by *transcript*
  - PREFIX.mapping_report.tsv     # per query symbol: chosen names & segment counts
  - PREFIX.missing_genes.txt      # NOT_FOUND / SKIPPED_MT queries
"""

import argparse, csv, gzip, re
from collections import defaultdict

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gff", required=True, help="Stickeback GFF3/GTF (.gz ok)")
    p.add_argument("--genes", required=True,
                   help="Human gene symbols file (one per line) or a comma-separated string")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--feature", choices=["cds", "exon", "gene"], default="cds",
                   help="Build from which feature type (default: cds)")
    p.add_argument("--include-mt", action="store_true",
                   help="Include mtDNA genes via MT-* → fish mapping (COX1/ND1/ATP6...)")
    p.add_argument("--mt-chr-name", default=None,
                   help="Mitochondrial chr name in reference (e.g., chrM/MT); optional")
    p.add_argument("--paralog_policy", choices=["all", "prefer_a", "prefer_b"], default="all",
                   help="Handle a/b duplicates (default: all)")
    p.add_argument("--ortholog-map", default=None,
                   help="TSV with columns: human_symbol<TAB>stickleback_gene_or_transcript (header ok)")
    p.add_argument("--faidx", default=None,
                   help=".fai file to restrict output to present contigs (e.g., nuclear-only)")
    return p.parse_args()

def read_gene_list(path_or_csv):
    try:
        with open(path_or_csv) as f:
            txt = f.read()
    except FileNotFoundError:
        txt = path_or_csv
    raw = re.split(r"[,\n\r\t ]+", txt.strip())
    return [x for x in (s.strip() for s in raw) if x]

def parse_attr(attr_text):
    d = {}
    for chunk in re.split(r";\s*", attr_text.strip().strip(";")):
        if not chunk: continue
        if "=" in chunk:
            k, v = chunk.split("=", 1)
        elif " " in chunk:
            k, v = chunk.split(" ", 1)
        else:
            continue
        d[k.strip()] = v.strip().strip('"')
    return d

def read_fai_contigs(fai_path):
    if not fai_path: return None
    contigs = set()
    with open(fai_path) as f:
        for line in f:
            if not line.strip(): continue
            contigs.add(line.split("\t", 1)[0])
    return contigs

def normalize_symbol(s): return s.lower()
def is_a_copy(name): return bool(re.search(r"[aA]$", name))
def is_b_copy(name): return bool(re.search(r"[bB]$", name))

MT_MAP = {
    "MT-CO1":"COX1","MT-CO2":"COX2","MT-CO3":"COX3",
    "MT-CYB":"CYTB","MT-ATP6":"ATP6","MT-ATP8":"ATP8",
    "MT-ND1":"ND1","MT-ND2":"ND2","MT-ND3":"ND3",
    "MT-ND4":"ND4","MT-ND4L":"ND4L","MT-ND5":"ND5","MT-ND6":"ND6"
}

def merge_intervals(blocks):
    if not blocks: return []
    blocks = sorted(blocks, key=lambda x: (x[0], x[1], x[2]))
    out = []
    cur_chr, cur_s, cur_e, cur_st = blocks[0]
    for c, s, e, st in blocks[1:]:
        if c == cur_chr and s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            out.append((cur_chr, cur_s, cur_e, cur_st))
            cur_chr, cur_s, cur_e, cur_st = c, s, e, st
    out.append((cur_chr, cur_s, cur_e, cur_st))
    return out

def load_ortholog_map(tsv_path):
    m = defaultdict(set)
    if not tsv_path: return m
    with open(tsv_path, newline="") as f:
        for i, row in enumerate(f):
            row = row.rstrip("\n")
            if not row: continue
            parts = row.split("\t")
            if len(parts) < 2: continue
            h, s = parts[0].strip(), parts[1].strip()
            if i == 0 and ("human" in h.lower() or "symbol" in h.lower()):
                continue
            if h and s:
                m[normalize_symbol(h)].add(s)
    return m

def main():
    args = parse_args()
    want_feat = {"cds":"CDS","exon":"exon","gene":"gene"}[args.feature]
    contig_filter = read_fai_contigs(args.faidx)

    query_symbols = read_gene_list(args.genes)
    q_mt = [q for q in query_symbols if q.upper().startswith("MT-")]
    q_nu = [q for q in query_symbols if not q.upper().startswith("MT-")]

    allowed = set(normalize_symbol(x) for x in q_nu)
    if args.include_mt:
        for mt_h in q_mt:
            fish = MT_MAP.get(mt_h.upper())
            if fish:
                allowed.add(normalize_symbol(fish))

    # Pass 1: collect gene/transcript relations
    gene_id_to_name = {}
    tx_id_to_gene_id = {}
    with open_text(args.gff) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9: continue
            chrom, src, ftype, start, end, score, strand, phase, attr = cols
            if ftype not in ("gene","mRNA","transcript"): continue
            ad = parse_attr(attr)
            if ftype == "gene":
                gid = ad.get("ID") or ad.get("gene_id")
                gname = ad.get("gene_name") or ad.get("Name") or ad.get("gene") or gid
                if gid: gene_id_to_name[gid] = gname
            else:
                tid = ad.get("ID") or ad.get("transcript_id")
                parent = ad.get("Parent")
                if tid and parent:
                    parent = parent.split(",")[0]
                    tx_id_to_gene_id[tid] = parent

    # Pass 2: collect parts
    gene_parts = defaultdict(list)
    tx_parts   = defaultdict(list)
    name_index = defaultdict(set)
    tx_name_index = defaultdict(set)

    with open_text(args.gff) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9: continue
            chrom, src, ftype, start, end, score, strand, phase, attr = cols
            if ftype != want_feat: continue
            ad = parse_attr(attr)
            s0 = int(start) - 1
            e1 = int(end)

            gene_name = None
            if "gene_name" in ad: gene_name = ad["gene_name"]
            elif "gene" in ad:   gene_name = ad["gene"]
            elif "Name" in ad:   gene_name = ad["Name"]

            gid = ad.get("gene_id")
            tid = ad.get("transcript_id")
            parent = ad.get("Parent")

            if not gene_name and gid:
                gene_name = gene_id_to_name.get(gid)
            if not gene_name and tid:
                gid2 = tx_id_to_gene_id.get(tid)
                if gid2: gene_name = gene_id_to_name.get(gid2)
            if not gene_name and parent:
                par0 = parent.split(",")[0]
                gene_name = gene_id_to_name.get(par0) or gene_id_to_name.get(tx_id_to_gene_id.get(par0, ""))
            if not gene_name:
                gene_name = gid or ad.get("ID") or parent or "UNKNOWN"

            if contig_filter is not None and chrom not in contig_filter:
                continue

            gene_parts[gene_name].append((chrom, s0, e1, strand))
            name_index[normalize_symbol(gene_name)].add(gene_name)

            tx_id = tid or ad.get("ID") or parent
            if tx_id:
                tx_parts[tx_id].append((chrom, s0, e1, strand))
                tx_name_index[normalize_symbol(tx_id)].add(tx_id)

    # txid -> gene_name
    txid_to_gene_name = {}
    for txid in tx_parts.keys():
        gid = tx_id_to_gene_id.get(txid)
        if gid:
            gname = gene_id_to_name.get(gid)
            if gname:
                txid_to_gene_name[txid] = gname

    # Direct matches
    ci_map = {}
    for sym in set(q_nu):
        ci_map[sym] = sorted(name_index.get(normalize_symbol(sym), set()))

    # MT mapping
    if args.include_mt:
        for mt_h in q_mt:
            fish = MT_MAP.get(mt_h.upper())
            if not fish: continue
            ci_map.setdefault(fish, [])
            ci_map[fish] = sorted(set(ci_map[fish]) |
                                  set(name_index.get(normalize_symbol(fish), set())))

    # Apply ortholog map
    ortho = load_ortholog_map(args.ortholog_map)
    if ortho:
        for human_sym in list(set(q_nu + ([MT_MAP.get(h.upper(), h) for h in q_mt] if args.include_mt else []))):
            low_h = normalize_symbol(human_sym)
            aliases = ortho.get(low_h, set())
            for alias in aliases:
                low_s = normalize_symbol(alias)
                real_gene_names = name_index.get(low_s)
                if real_gene_names:
                    ci_map[human_sym] = sorted(set(ci_map.get(human_sym, [])) | set(real_gene_names))
                    continue
                possible_tx = tx_name_index.get(low_s)
                added_any = False
                if possible_tx:
                    for txid in possible_tx:
                        gname = txid_to_gene_name.get(txid)
                        if gname:
                            ci_map[human_sym] = sorted(set(ci_map.get(human_sym, [])) | {gname})
                            added_any = True
                if added_any:
                    continue
                for txid in tx_parts.keys():
                    if normalize_symbol(txid) == low_s:
                        gname = txid_to_gene_name.get(txid)
                        if gname:
                            ci_map[human_sym] = sorted(set(ci_map.get(human_sym, [])) | {gname})
                        break

    def choose_candidates(cands):
        if args.paralog_policy == "all":
            return cands
        aN = [n for n in cands if is_a_copy(n)]
        bN = [n for n in cands if is_b_copy(n)]
        if args.paralog_policy == "prefer_a":
            return aN if aN else cands
        else:
            return bN if bN else cands

    chosen = {}
    for q in query_symbols:
        if q.upper().startswith("MT-") and not args.include_mt:
            chosen[q] = []
            continue
        base = MT_MAP.get(q.upper(), q) if args.include_mt and q.upper().startswith("MT-") else q
        cands = ci_map.get(base, [])
        chosen[q] = choose_candidates(cands) if cands else []

    cds_parts_bed = args.out_prefix + ".cds_parts.bed"
    gene_bed      = args.out_prefix + ".cds_by_gene.bed"
    tx_bed        = args.out_prefix + ".cds_by_tx.bed"
    report_tsv    = args.out_prefix + ".mapping_report.tsv"
    missing_txt   = args.out_prefix + ".missing_genes.txt"

    with open(cds_parts_bed, "w") as fo:
        for q, names in chosen.items():
            if q.upper().startswith("MT-") and not args.include_mt:
                continue
            for gname in names:
                for chrom, s, e, st in gene_parts.get(gname, []):
                    fo.write(f"{chrom}\t{s}\t{e}\t{gname}\t.\t{st}\n")

    with open(gene_bed, "w") as fo:
        for q, names in chosen.items():
            if q.upper().startswith("MT-") and not args.include_mt:
                continue
            for gname in names:
                merged = merge_intervals(gene_parts.get(gname, []))
                for chrom, s, e, st in merged:
                    fo.write(f"{chrom}\t{s}\t{e}\t{gname}\t.\t{st}\n")

    with open(tx_bed, "w") as fo:
        for txid, parts in tx_parts.items():
            merged = merge_intervals(parts)
            for chrom, s, e, st in merged:
                fo.write(f"{chrom}\t{s}\t{e}\t{txid}\t.\t{st}\n")

    missing = []
    with open(report_tsv, "w", newline="") as fo:
        wr = csv.writer(fo, delimiter="\t")
        wr.writerow(["query_symbol","chosen_gene_name","n_segments_raw","note"])
        for q in query_symbols:
            names = chosen.get(q, [])
            if q.upper().startswith("MT-") and not args.include_mt:
                wr.writerow([q, "", 0, "SKIPPED_MT"])
                continue
            if not names:
                missing.append(q)
                wr.writerow([q, "", 0, "NOT_FOUND"])
            else:
                note = []
                if q.upper().startswith("MT-"): note.append("MT_map")
                if args.paralog_policy != "all": note.append(args.paralog_policy)
                for gname in names:
                    nseg = len(gene_parts.get(gname, []))
                    wr.writerow([q, gname, nseg, ";".join(note) if note else ""])

    with open(missing_txt, "w") as fo:
        for x in missing:
            fo.write(x + "\n")

    print("✅ Done.")
    print(" -", cds_parts_bed)
    print(" -", gene_bed)
    print(" -", tx_bed)
    print(" -", report_tsv)
    print(" -", missing_txt)

if __name__ == "__main__":
    main()




#diff_oxphos_sets_mapped.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, re
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description="Compare human OXPHOS sets vs stickleback genes (supports ortholog map)")
    p.add_argument("--span-bed", required=True, help="SPAN bed (one gene per row)")
    p.add_argument("--out-prefix", required=True, help="Output prefix")
    p.add_argument("--ignore-mt", action="store_true", help="Ignore mitochondrial genes (MT-*) in human sets")
    p.add_argument("--keep-ab", action="store_true", help="Do NOT collapse teleost a/b copies")
    p.add_argument("--ortholog-map", help="TSV: human_symbol<TAB>stickleback_gene_or_transcript (header ok)")
    return p.parse_args()

def S(s):  # CSV -> UPPER set
    return set(x.strip().upper() for x in s.split(",") if x.strip())

# ---- Human sets ----
HUMAN = {
"OXPHOS (all)": S("""ACAD9, AIFM1, ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E, ATP5IF1, ATP5MC1, ATP5MC2, ATP5MC3, ATP5MD, ATP5ME, ATP5MF, ATP5MG, ATP5MPL, ATP5PB, ATP5PD, ATP5PF, ATP5PO, ATPAF1, ATPAF2, ATPSCKMT, BCS1L, CEP89, CMC1, CMC2, COA1, COA3, COA4, COA5, COA6, COA7, COA8, COX10, COX11, COX14, COX15, COX16, COX17, COX18, COX19, COX20, COX4I1, COX4I2, COX5A, COX5B, COX6A1, COX6A2, COX6B1, COX6B2, COX6C, COX7A1, COX7A2, COX7A2L, COX7B, COX7B2, COX7C, COX8A, COX8C, CYC1, CYCS, DMAC1, DMAC2, DMAC2L, ECSIT, FMC1, FOXRED1, HCCS, HIGD1A, HIGD2A, LYRM2, LYRM7, MT-ATP6, MT-ATP8, MT-CO1, MT-CO2, MT-CO3, MT-CYB, MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6, NDUFA1, NDUFA10, NDUFA11, NDUFA12, NDUFA13, NDUFA2, NDUFA3, NDUFA4, NDUFA5, NDUFA6, NDUFA7, NDUFA8, NDUFA9, NDUFAB1, NDUFAF1, NDUFAF2, NDUFAF3, NDUFAF4, NDUFAF5, NDUFAF6, NDUFAF7, NDUFAF8, NDUFB1, NDUFB10, NDUFB11, NDUFB2, NDUFB3, NDUFB4, NDUFB5, NDUFB6, NDUFB7, NDUFB8, NDUFB9, NDUFC1, NDUFC2, NDUFS1, NDUFS2, NDUFS3, NDUFS4, NDUFS5, NDUFS6, NDUFS7, NDUFS8, NDUFV1, NDUFV2, NDUFV3, NUBPL, PET100, PET117, PNKD, RAB5IF, SCO1, SCO2, SDHA, SDHAF1, SDHAF2, SDHAF3, SDHAF4, SDHB, SDHC, SDHD, SMIM20, SURF1, TACO1, TIMM21, TIMMDC1, TMEM126A, TMEM126B, TMEM177, TMEM186, TMEM70, TTC19, UQCC1, UQCC2, UQCC3, UQCR10, UQCR11, UQCRB, UQCRC1, UQCRC2, UQCRFS1, UQCRH, UQCRQ"""),
"OXPHOS subunits": S("""ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E, ATP5IF1, ATP5MC1, ATP5MC2, ATP5MC3, ATP5MD, ATP5ME, ATP5MF, ATP5MG, ATP5MPL, ATP5PB, ATP5PD, ATP5PF, ATP5PO, COX4I1, COX4I2, COX5A, COX5B, COX6A1, COX6A2, COX6B1, COX6B2, COX6C, COX7A1, COX7A2, COX7A2L, COX7B, COX7B2, COX7C, COX8A, COX8C, CYC1, CYCS, DMAC2L, HCCS, MT-ATP6, MT-ATP8, MT-CO1, MT-CO2, MT-CO3, MT-CYB, MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6, NDUFA1, NDUFA10, NDUFA11, NDUFA12, NDUFA13, NDUFA2, NDUFA3, NDUFA4, NDUFA5, NDUFA6, NDUFA7, NDUFA8, NDUFA9, NDUFAB1, NDUFB1, NDUFB10, NDUFB11, NDUFB2, NDUFB3, NDUFB4, NDUFB5, NDUFB6, NDUFB7, NDUFB8, NDUFB9, NDUFC1, NDUFC2, NDUFS1, NDUFS2, NDUFS3, NDUFS4, NDUFS5, NDUFS6, NDUFS7, NDUFS8, NDUFV1, NDUFV2, NDUFV3, SDHA, SDHB, SDHC, SDHD, UQCR10, UQCR11, UQCRB, UQCRC1, UQCRC2, UQCRFS1, UQCRH, UQCRQ"""),
"OXPHOS assembly factors": S("""ACAD9, AIFM1, ATPAF1, ATPAF2, ATPSCKMT, BCS1L, CEP89, CMC1, CMC2, COA1, COA3, COA4, COA5, COA6, COA7, COA8, COX10, COX11, COX14, COX15, COX16, COX17, COX18, COX19, COX20, COX7A2L, DMAC1, DMAC2, ECSIT, FMC1, FOXRED1, HIGD1A, HIGD2A, LYRM2, LYRM7, NDUFAF1, NDUFAF2, NDUFAF3, NDUFAF4, NDUFAF5, NDUFAF6, NDUFAF7, NDUFAF8, NUBPL, PET100, PET117, PNKD, RAB5IF, SCO1, SCO2, SDHAF1, SDHAF2, SDHAF3, SDHAF4, SMIM20, SURF1, TACO1, TIMM21, TIMMDC1, TMEM126A, TMEM126B, TMEM177, TMEM186, TMEM70, TTC19, UQCC1, UQCC2, UQCC3"""),
"Complex I": S("""ACAD9, AIFM1, COA1, DMAC1, DMAC2, ECSIT, FOXRED1, LYRM2, MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6, NDUFA1, NDUFA10, NDUFA11, NDUFA12, NDUFA13, NDUFA2, NDUFA3, NDUFA5, NDUFA6, NDUFA7, NDUFA8, NDUFA9, NDUFAB1, NDUFAF1, NDUFAF2, NDUFAF3, NDUFAF4, NDUFAF5, NDUFAF6, NDUFAF7, NDUFAF8, NDUFB1, NDUFB10, NDUFB11, NDUFB2, NDUFB3, NDUFB4, NDUFB5, NDUFB6, NDUFB7, NDUFB8, NDUFB9, NDUFC1, NDUFC2, NDUFS1, NDUFS2, NDUFS3, NDUFS4, NDUFS5, NDUFS6, NDUFS7, NDUFS8, NDUFV1, NDUFV2, NDUFV3, NUBPL, TIMMDC1, TMEM126A, TMEM126B, TMEM186, TMEM70"""),
"CI subunits": S("""MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6, NDUFA1, NDUFA10, NDUFA11, NDUFA12, NDUFA13, NDUFA2, NDUFA3, NDUFA5, NDUFA6, NDUFA7, NDUFA8, NDUFA9, NDUFAB1, NDUFB1, NDUFB10, NDUFB11, NDUFB2, NDUFB3, NDUFB4, NDUFB5, NDUFB6, NDUFB7, NDUFB8, NDUFB9, NDUFC1, NDUFC2, NDUFS1, NDUFS2, NDUFS3, NDUFS4, NDUFS5, NDUFS6, NDUFS7, NDUFS8, NDUFV1, NDUFV2, NDUFV3"""),
"CI assembly factors": S("""ACAD9, AIFM1, COA1, DMAC1, DMAC2, ECSIT, FOXRED1, LYRM2, NDUFAF1, NDUFAF2, NDUFAF3, NDUFAF4, NDUFAF5, NDUFAF6, NDUFAF7, NDUFAF8, NUBPL, TIMMDC1, TMEM126A, TMEM126B, TMEM186, TMEM70"""),
"Complex II": S("""SDHA, SDHAF1, SDHAF2, SDHAF3, SDHAF4, SDHB, SDHC, SDHD"""),
"CII subunits": S("""SDHA, SDHB, SDHC, SDHD"""),
"CII assembly factors": S("""SDHAF1, SDHAF2, SDHAF3, SDHAF4"""),
"Complex III": S("""BCS1L, CYC1, LYRM7, MT-CYB, TTC19, UQCC1, UQCC2, UQCC3, UQCR10, UQCR11, UQCRB, UQCRC1, UQCRC2, UQCRFS1, UQCRH, UQCRQ"""),
"CIII subunits": S("""CYC1, MT-CYB, UQCR10, UQCR11, UQCRB, UQCRC1, UQCRC2, UQCRFS1, UQCRH, UQCRQ"""),
"CIII assembly factors": S("""BCS1L, LYRM7, TTC19, UQCC1, UQCC2, UQCC3"""),
"Complex IV": S("""CEP89, CMC1, CMC2, COA1, COA3, COA4, COA5, COA6, COA7, COA8, COX10, COX11, COX14, COX15, COX16, COX17, COX18, COX19, COX20, COX4I1, COX4I2, COX5A, COX5B, COX6A1, COX6A2, COX6B1, COX6B2, COX6C, COX7A1, COX7A2, COX7A2L, COX7B, COX7B2, COX7C, COX8A, COX8C, HIGD1A, MT-CO1, MT-CO2, MT-CO3, NDUFA4, PET100, PET117, PNKD, SCO1, SCO2, SMIM20, SURF1, TACO1, TIMM21, TMEM177"""),
"CIV subunits": S("""COX4I1, COX4I2, COX5A, COX5B, COX6A1, COX6A2, COX6B1, COX6B2, COX6C, COX7A1, COX7A2, COX7A2L, COX7B, COX7B2, COX7C, COX8A, COX8C, MT-CO1, MT-CO2, MT-CO3, NDUFA4"""),
"CIV assembly factors": S("""CEP89, CMC1, CMC2, COA1, COA3, COA4, COA5, COA6, COA7, COX10, COX11, COX14, COX15, COX16, COX17, COX18, COX19, COX20, HIGD1A, PET100, PET117, PNKD, SCO1, SCO2, SMIM20, SURF1, TACO1, TIMM21, TMEM177"""),
"Complex V": S("""ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E, ATP5IF1, ATP5MC1, ATP5MC2, ATP5MC3, ATP5MD, ATP5ME, ATP5MF, ATP5MG, ATP5MPL, ATP5PB, ATP5PD, ATP5PF, ATP5PO, ATPAF1, ATPAF2, ATPSCKMT, DMAC2L, FMC1, MT-ATP6, MT-ATP8, TMEM70"""),
"CV subunits": S("""ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E, ATP5IF1, ATP5MC1, ATP5MC2, ATP5MC3, ATP5MD, ATP5ME, ATP5MF, ATP5MG, ATP5MPL, ATP5PB, ATP5PD, ATP5PF, ATP5PO, DMAC2L, MT-ATP6, MT-ATP8"""),
"CV assembly factors": S("""ATPAF1, ATPAF2, ATPSCKMT, FMC1, TMEM70"""),
}
ALL_HUMAN = set().union(*HUMAN.values())

def is_mt(g): return g.startswith("MT-")

def load_ortholog_map(tsv):
    if not tsv: return {}
    fish2human = defaultdict(set)
    with open(tsv) as f:
        first = f.readline()
        if first and ("\t" in first):
            cols = first.rstrip("\n").split("\t")
            if not (("human" in cols[0].lower()) or ("symbol" in cols[0].lower())):
                h, s = cols[0].strip().upper(), cols[1].strip().upper() if len(cols) > 1 else ""
                if h and s: fish2human[s].add(h)
        else:
            return {}
        for line in f:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2: continue
            h, s = parts[0].strip().upper(), parts[1].strip().upper()
            if h and s: fish2human[s].add(h)
    return fish2human

def canon_name_for_human_match(x, keep_ab=False):
    out = set()
    X = x.upper()
    if X in ALL_HUMAN: out.add(X)
    if (not keep_ab) and re.search(r"[AB]$", X) and len(X) > 1:
        base = X[:-1]
        if base in ALL_HUMAN: out.add(base)
    return out

def read_fish_genes_as_humans(span_bed, fish2human_map, keep_ab=False):
    present = set()
    with open(span_bed) as f:
        for line in f:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4: continue
            fish = parts[3].strip().upper()
            if fish in fish2human_map:
                present |= fish2human_map[fish]
            if (not keep_ab) and re.search(r"[AB]$", fish) and len(fish) > 1:
                base = fish[:-1]
                if base in fish2human_map:
                    present |= fish2human_map[base]
            present |= canon_name_for_human_match(fish, keep_ab=keep_ab)
    return present

def main():
    args = parse_args()
    human_sets = {k:set(v) for k,v in HUMAN.items()}
    if args.ignore_mt:
        for k in human_sets:
            human_sets[k] = {g for g in human_sets[k] if not is_mt(g)}
    fish2human = load_ortholog_map(args.ortholog_map) if args.ortholog_map else {}
    present_humans = read_fish_genes_as_humans(args.span_bed, fish2human, keep_ab=args.keep_ab)

    summary = args.out_prefix + ".summary.tsv"
    with open(summary, "w") as fo:
        fo.write("set\thuman_n\tfish_present_n\tmissing_n\n")
        for set_name, H in human_sets.items():
            present = sorted(H & present_humans)
            missing = sorted(H - present_humans)
            fo.write(f"{set_name}\t{len(H)}\t{len(present)}\t{len(missing)}\n")
            tag = re.sub(r"[^A-Za-z0-9]+","_", set_name)
            with open(args.out_prefix + f".{tag}.missing.txt", "w") as fm:
                fm.write("\n".join(missing) + ("\n" if missing else ""))

    print("✅ Wrote:", summary)
    print("   + per-set missing lists alongside the summary.")

if __name__ == "__main__":
    main()



#scan_alias_hits.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, gzip, re, csv

def open_text(p): return gzip.open(p,"rt") if p.endswith(".gz") else open(p)

def parse_attr(s):
    d={}
    for x in re.split(r";\s*", s.strip().strip(";")):
        if not x: continue
        if "=" in x: k,v=x.split("=",1)
        elif " " in x: k,v=x.split(" ",1)
        else: continue
        d[k.strip()]=v.strip().strip('"')
    return d

def main():
    ap=argparse.ArgumentParser(description="Scan GFF for alias hits of missing genes")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--aliases", required=True, help="TSV: human<TAB>alias1,alias2")
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    A={}
    with open(args.aliases) as f:
        for ln in f:
            if not ln.strip() or "\t" not in ln: continue
            h, al = ln.rstrip("\n").split("\t",1)
            A[h.strip().upper()] = [x.strip().upper() for x in al.split(",") if x.strip()]

    hits=[]
    with open_text(args.gff) as fh:
        for ln in fh:
            if not ln or ln.startswith("#"): continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9: continue
            chrom,_,ftype,start,end,_,strand,_,attr = c
            if ftype not in ("gene","mRNA","transcript","CDS","exon"): continue
            ad = parse_attr(attr)
            vals = []
            for k in ("gene_name","gene","Name","ID","transcript_id","Parent"):
                v = ad.get(k)
                if v: vals.append(v.upper())
            text = " ".join(vals)

            for h, aliases in A.items():
                for al in aliases:
                    if re.search(rf"\b{re.escape(al)}(\b|-201\b)", text, flags=re.I):
                        hits.append([h, al, ftype, chrom, start, end, strand,
                                     ad.get("ID",""), ad.get("Parent",""),
                                     ad.get("gene_name","") or ad.get("Name","") or ad.get("gene","")])

    with open(args.out,"w",newline="") as fo:
        wr=csv.writer(fo, delimiter="\t")
        wr.writerow(["human","alias_matched","feature","chrom","start","end","strand","ID","Parent","display_name"])
        wr.writerows(hits)

if __name__ == "__main__":
    main()

