#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build OXPHOS BED from a human gene list + stickleback GFF/GTF.

Outputs (with --out-prefix PREFIX):
  - PREFIX.cds_parts.bed          # 原始 CDS 片段
  - PREFIX.cds_by_gene.bed        # gene 级（CDS 并集合并后的多段）
  - PREFIX.cds_by_tx.bed          # transcript 级（CDS 并集合并后的多段）
  - PREFIX.mapping_report.tsv     # 每个查询基因的命中详情
  - PREFIX.missing_genes.txt      # 未命中的（按策略）基因

Notes
- 兼容 .gff3/.gtf，支持 .gz 压缩
- 宽松匹配：大小写不敏感；可选择 a/b 拷贝策略(all|prefer_a|prefer_b)
- 可选 --include-mt 把 MT-*（人类命名）映射为鱼类常用命名（COX1、ND1 等）
- 可选 --ortholog-map human\tstickleback 两列表 TSV 扩充同源映射
- 可选 --faidx 仅保留在该 .fai 中存在的染色体（方便 nuclear-only 参考）
"""

import argparse, csv, gzip, re, sys
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
                   help="TSV with columns: human_symbol<TAB>stickleback_symbol (header ok)")
    p.add_argument("--faidx", default=None,
                   help=".fai file to restrict output to present contigs (e.g., nuclear-only)")
    return p.parse_args()

def read_gene_list(path_or_csv):
    # allow file path or CSV string
    try:
        with open(path_or_csv) as f:
            txt = f.read()
    except FileNotFoundError:
        txt = path_or_csv
    raw = re.split(r"[,\n\r\t ]+", txt.strip())
    return [x for x in (s.strip() for s in raw) if x]

def normalize_symbol(s): return s.lower()
def is_a_copy(name): return bool(re.search(r"[aA]$", name))
def is_b_copy(name): return bool(re.search(r"[bB]$", name))

MT_MAP = {
    "MT-CO1":"COX1","MT-CO2":"COX2","MT-CO3":"COX3",
    "MT-CYB":"CYTB","MT-ATP6":"ATP6","MT-ATP8":"ATP8",
    "MT-ND1":"ND1","MT-ND2":"ND2","MT-ND3":"ND3",
    "MT-ND4":"ND4","MT-ND4L":"ND4L","MT-ND5":"ND5","MT-ND6":"ND6"
}

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

def load_ortholog_map(tsv_path):
    m = defaultdict(set)
    if not tsv_path: return m
    with open(tsv_path, newline="") as f:
        rd = csv.reader(f, delimiter="\t")
        # try to read header; if looks wrong, rewind
        pos = f.tell()
        header = next(rd, None)
        if not header or len(header) < 2 or ("human" not in header[0].lower() and "stick" not in header[1].lower()):
            f.seek(0); rd = csv.reader(f, delimiter="\t")
        else:
            f.seek(pos)
            rd = csv.reader(f, delimiter="\t")
            next(rd, None)
        for row in rd:
            if len(row) < 2: continue
            h, s = row[0].strip(), row[1].strip()
            if h and s:
                m[normalize_symbol(h)].add(s)
    return m

def read_fai_contigs(fai_path):
    if not fai_path: return None
    contigs = set()
    with open(fai_path) as f:
        for line in f:
            if not line.strip(): continue
            contigs.add(line.split("\t", 1)[0])
    return contigs

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

def main():
    args = parse_args()
    want_feat = {"cds":"CDS","exon":"exon","gene":"gene"}[args.feature]
    contig_filter = read_fai_contigs(args.faidx)

    # split query list into nuclear & MT
    query_symbols = read_gene_list(args.genes)
    query_set = [q for q in query_symbols if q]  # keep order
    q_mt = [q for q in query_set if q.upper().startswith("MT-")]
    q_nu = [q for q in query_set if not q.upper().startswith("MT-")]

    # allowed matching names (lowercased)
    allowed = set(normalize_symbol(x) for x in q_nu)
    if args.include_mt:
        for mt_h in q_mt:
            fish = MT_MAP.get(mt_h.upper())
            if fish: allowed.add(normalize_symbol(fish))

    # -------- First pass: build maps gene/mRNA relationships --------
    gene_id_to_name = {}          # gene ID -> gene symbol (prefer gene_name/Name)
    tx_id_to_gene_id = {}         # transcript ID -> parent gene ID
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
            else:  # mRNA/transcript
                tid = ad.get("ID") or ad.get("transcript_id")
                parent = ad.get("Parent")
                if tid and parent:
                    # Parent might be a list (comma-separated) in some GFFs
                    parent = parent.split(",")[0]
                    tx_id_to_gene_id[tid] = parent

    # -------- Second pass: collect parts for features --------
    gene_parts = defaultdict(list)   # gene_name -> [(chr,s,e,strand)]
    tx_parts   = defaultdict(list)   # transcript_id -> [(chr,s,e,strand)]
    name_index = defaultdict(set)    # lower(gene_name) -> set(real names seen)

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

            # Identify gene_name for this feature
            gene_name = None
            # Direct attributes (GTF often has these)
            if "gene_name" in ad: gene_name = ad["gene_name"]
            elif "gene" in ad:   gene_name = ad["gene"]
            elif "Name" in ad:   gene_name = ad["Name"]

            # If not, try mapping via IDs/Parent chain
            gid = ad.get("gene_id")
            tid = ad.get("transcript_id")
            parent = ad.get("Parent")

            if not gene_name and gid:
                gene_name = gene_id_to_name.get(gid)
            if not gene_name and tid:
                gid2 = tx_id_to_gene_id.get(tid)
                if gid2: gene_name = gene_id_to_name.get(gid2)
            if not gene_name and parent:
                # parent may be transcript or gene
                par0 = parent.split(",")[0]
                gene_name = gene_id_to_name.get(par0)
                if not gene_name:
                    gid2 = tx_id_to_gene_id.get(par0)
                    if gid2: gene_name = gene_id_to_name.get(gid2)

            # If still missing, last resort: use gene_id/ID itself
            if not gene_name:
                gene_name = gid or ad.get("ID") or parent or "UNKNOWN"

            # contig filter (e.g., nuclear-only)
            if contig_filter is not None and chrom not in contig_filter:
                continue

            gene_parts[gene_name].append((chrom, s0, e1, strand))
            name_index[normalize_symbol(gene_name)].add(gene_name)

            # transcript parts
            tx_id = tid or ad.get("ID") or parent
            if tx_id:
                tx_parts[tx_id].append((chrom, s0, e1, strand))

    # -------- Build candidate mapping for queries --------
    ortho = load_ortholog_map(args.ortholog_map)  # human(lower) -> set(stickle symbols)
    ci_map = {}  # query symbol (possibly fish name for MT) -> candidate observed names
    # nuclear queries direct
    for sym in set(q_nu):
        ci_map[sym] = sorted(name_index.get(normalize_symbol(sym), set()))
    # add ortholog expansions
    for sym in list(ci_map.keys()):
        low = normalize_symbol(sym)
        for stick in ortho.get(low, []):
            ci_map[sym] = sorted(set(ci_map.get(sym, [])) |
                                 set(name_index.get(normalize_symbol(stick), set())))
    # MT mapping if requested
    if args.include_mt:
        for mt_h in q_mt:
            fish = MT_MAP.get(mt_h.upper())
            if not fish: continue
            ci_map.setdefault(fish, [])
            ci_map[fish] = sorted(set(ci_map[fish]) |
                                  set(name_index.get(normalize_symbol(fish), set())))

    # -------- Choose per query by paralog policy --------
    chosen = {}  # query_symbol -> list of chosen gene_name(s)
    def choose_candidates(cands):
        if args.paralog_policy == "all":
            return cands
        aN = [n for n in cands if is_a_copy(n)]
        bN = [n for n in cands if is_b_copy(n)]
        if args.paralog_policy == "prefer_a":
            return aN if aN else cands
        else:
            return bN if bN else cands

    for q in query_set:
        if q.upper().startswith("MT-") and not args.include_mt:
            chosen[q] = []  # will be marked SKIPPED_MT
            continue
        base = q
        if args.include_mt and q.upper().startswith("MT-"):
            base = MT_MAP.get(q.upper(), q)
        cands = ci_map.get(base, [])
        chosen[q] = choose_candidates(cands) if cands else []

    # -------- Write outputs --------
    cds_parts_bed = args.out_prefix + ".cds_parts.bed"
    gene_bed      = args.out_prefix + ".cds_by_gene.bed"
    tx_bed        = args.out_prefix + ".cds_by_tx.bed"
    report_tsv    = args.out_prefix + ".mapping_report.tsv"
    missing_txt   = args.out_prefix + ".missing_genes.txt"

    # 1) raw CDS/exon/gene parts
    with open(cds_parts_bed, "w") as fo:
        for q, names in chosen.items():
            # skip MT queries if not included
            if q.upper().startswith("MT-") and not args.include_mt:
                continue
            for gname in names:
                for chrom, s, e, st in gene_parts.get(gname, []):
                    fo.write(f"{chrom}\t{s}\t{e}\t{gname}\t.\t{st}\n")

    # 2) merged by gene
    with open(gene_bed, "w") as fo:
        for q, names in chosen.items():
            if q.upper().startswith("MT-") and not args.include_mt:
                continue
            for gname in names:
                merged = merge_intervals(gene_parts.get(gname, []))
                for chrom, s, e, st in merged:
                    fo.write(f"{chrom}\t{s}\t{e}\t{gname}\t.\t{st}\n")

    # 3) merged by transcript
    with open(tx_bed, "w") as fo:
        for txid, parts in tx_parts.items():
            merged = merge_intervals(parts)
            for chrom, s, e, st in merged:
                fo.write(f"{chrom}\t{s}\t{e}\t{txid}\t.\t{st}\n")

    # 4) report + missing
    missing = []
    with open(report_tsv, "w", newline="") as fo:
        wr = csv.writer(fo, delimiter="\t")
        wr.writerow(["query_symbol","chosen_gene_name_or_id","n_segments_raw","note"])
        for q in query_set:
            names = chosen.get(q, [])
            if q.upper().startswith("MT-") and not args.include_mt:
                wr.writerow([q, "", 0, "SKIPPED_MT"])
                continue
            if not names:
                missing.append(q)
                wr.writerow([q, "", 0, "NOT_FOUND"])
            else:
                note = []
                if q.upper().startswith("MT-"):
                    note.append("MT_map")
                if args.paralog_policy != "all":
                    note.append(args.paralog_policy)
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



BASE="/work/cyu/oxphos_from_ref_no_biomart"
OUTDIR="$BASE/06_igv"
GENE_LIST="$BASE/human_oxphos_symbols.txt"      # 之前我给你的那份
GFF="/work/cyu/stickleback_v5.gff3.gz"
FA="/work/cyu/stickleback_nuclear_only.fa"
mkdir -p "$OUTDIR"
[ -f "${FA}.fai" ] || samtools faidx "$FA"

python "$BASE/build_oxphos_bed.py" \
  --gff "$GFF" \
  --genes "$GENE_LIST" \
  --out-prefix "$OUTDIR/nuOXPHOS_v5" \
  --feature cds \
  --paralog_policy all \
  --faidx "${FA}.fai"

# 排序（按 nuclear-only 的 .fai 顺序）
bedtools sort -faidx "${FA}.fai" \
  -i "$OUTDIR/nuOXPHOS_v5.cds_by_gene.bed" \
  > "$OUTDIR/nuOXPHOS_v5.cds_by_gene.sorted.bed"

# 小体检（MT-* 不会算缺失；会标 SKIPPED_MT）
echo "✔️ 覆盖到的核 OXPHOS 基因数："
cut -f4 "$OUTDIR/nuOXPHOS_v5.cds_by_gene.sorted.bed" | sort -u | wc -l

echo "✔️ 未命中的核基因："
awk 'NR>1 && $3==0 && $4=="NOT_FOUND"{print $1}' "$OUTDIR/nuOXPHOS_v5.mapping_report.tsv" \
  | grep -v '^MT-' | sort -u | tee "$OUTDIR/nuOXPHOS_v5.missing_nuclear.txt"

echo "✅ 输出主文件：$OUTDIR/nuOXPHOS_v5.cds_by_gene.sorted.bed"
