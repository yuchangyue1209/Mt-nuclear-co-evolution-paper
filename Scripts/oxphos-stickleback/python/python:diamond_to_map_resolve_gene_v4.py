#!/usr/bin/env python3
import argparse, gzip, re, sys
from collections import defaultdict

def opengz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_attr(s):
    d={}
    for x in re.split(r";\s*", s.strip().strip(";")):
        if not x: continue
        k,v = (x.split("=",1)+[""])[:2] if "=" in x else (x.split(" ",1)+[""])[:2]
        d[k.strip()] = v.strip().strip('"')
    return d

def build_maps(gff):
    geneid2name={}
    tx2gene={}
    prot2gene={}
    any2gene={}

    with opengz(gff) as fh:
        for ln in fh:
            if not ln or ln.startswith("#"): continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9: continue
            ftype = c[2]
            a = parse_attr(c[8])
            if ftype == "gene":
                gid = a.get("ID") or a.get("gene_id")
                gname = a.get("Name") or a.get("gene_name") or a.get("gene") or a.get("locus_tag") or gid
                if gid:
                    geneid2name[gid] = gname
                    any2gene[gid] = gname
                    if gname: any2gene[gname] = gname
            elif ftype in ("mRNA","transcript"):
                tid = a.get("ID")
                par = a.get("Parent")
                gname = None
                if par and par in geneid2name: gname = geneid2name[par]
                if tid and gname:
                    tx2gene[tid] = gname
                    any2gene[tid] = gname
            elif ftype == "CDS":
                pid = a.get("protein_id") or a.get("ID")
                par = a.get("Parent")
                gname = None
                if par in tx2gene: gname = tx2gene[par]
                elif par in geneid2name: gname = geneid2name[par]
                if pid and gname:
                    prot2gene[pid] = gname
                    any2gene[pid] = gname
                    # also strip version
                    pid2 = pid.split(".")[0]
                    any2gene[pid2] = gname
    return any2gene

def read_span_names(span_bed):
    names=set()
    with open(span_bed) as f:
        for ln in f:
            if not ln.strip(): continue
            c=ln.split("\t")
            if len(c)>=4: names.add(c[3])
    return names

def best_hits(diam_path):
    best={}
    with open(diam_path) as f:
        for ln in f:
            if not ln.strip(): continue
            q,s,*rest = ln.rstrip("\n").split("\t", 2)
            bits = float(ln.rstrip("\n").split("\t")[-1])
            if q not in best or bits > best[q][1]:
                best[q] = (s, bits)
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diamond", required=True, help="DIAMOND outfmt 6 (qseqid sseqid ... bitscore)")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--span", required=True)
    ap.add_argument("--out-map", required=True)
    args = ap.parse_args()

    any2gene = build_maps(args.gff)
    span_names = read_span_names(args.span)
    hits = best_hits(args.diamond)

    out = []
    for q,(s,_) in hits.items():
        s_norm = s.split()[0]
        s_norm2 = s_norm.split(".")[0]
        g = any2gene.get(s_norm) or any2gene.get(s_norm2)
        if not g and s_norm in span_names:
            g = s_norm
        if g:
            out.append((q, g))

    with open(args.out_map, "w") as w:
        w.write("human_symbol\tstickleback_name\n")
        for q,g in sorted(set(out)):
            w.write(f"{q}\t{g}\n")

if __name__ == "__main__":
    main()
