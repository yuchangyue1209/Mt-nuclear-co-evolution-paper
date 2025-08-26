#!/usr/bin/env python3
import argparse, sys
from collections import defaultdict

def parse_tblastn(tbl):
    # columns: qseqid sseqid pident length qlen sstart send evalue bitscore
    rows=[]
    with open(tbl) as f:
        for ln in f:
            if not ln.strip(): continue
            toks = ln.rstrip("\n").split("\t")
            q = toks[0]; s = toks[1]
            pident = float(toks[2]); length = int(toks[3]); qlen = int(toks[4])
            sstart = int(toks[5]); send = int(toks[6])
            evalue = float(toks[7]); bits = float(toks[8])
            rows.append((q,s,pident,length,qlen,sstart,send,evalue,bits))
    return rows

def cluster_intervals(ints, max_gap):
    if not ints: return []
    ints = sorted(ints)
    out=[list(ints[0])]
    for s,e in ints[1:]:
        if s - out[-1][1] <= max_gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s,e])
    return [tuple(x) for x in out]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tblastn", required=True)
    ap.add_argument("--span_out", required=True)
    ap.add_argument("--map_out", required=True)
    ap.add_argument("--gap", type=int, default=3000)
    ap.add_argument("--pad", type=int, default=200)
    ap.add_argument("--min_qcov", type=float, default=0.45)
    ap.add_argument("--min_pident", type=float, default=0.20)
    ap.add_argument("--mode", choices=["best","all"], default="best")
    args = ap.parse_args()

    rows = parse_tblastn(args.tblastn)

    # filter by coverage and identity
    keep=[r for r in rows if (r[3]/max(1,r[4])>=args.min_qcov and r[2]>=args.min_pident)]

    # group by (qseqid, sseqid)
    grp=defaultdict(list)
    strand_vote=defaultdict(int)
    score_sum=defaultdict(float)

    for q,s,pident,length,qlen,sstart,send,evalue,bits in keep:
        a=min(sstart,send); b=max(sstart,send)
        grp[(q,s)].append((a,b))
        score_sum[(q,s)] += bits
        strand_vote[(q,s)] += 1 if send>=sstart else -1

    # choose best group per query by summed bitscore
    best_by_q={}
    for (q,s),intervals in grp.items():
        if (q not in best_by_q) or (score_sum[(q,s)] > best_by_q[q][2]):
            best_by_q[q]=(s, intervals, score_sum[(q,s)])

    # write outputs
    with open(args.span_out,"w") as sb, open(args.map_out,"w") as mb:
        mb.write("human_symbol\tstickleback_name\n")
        out_count=0
        for q,(s,intervals,sumscore) in sorted(best_by_q.items()):
            merged = cluster_intervals(intervals, args.gap)
            # pick the longest merged block
            merged.sort(key=lambda x: x[1]-x[0], reverse=True)
            if not merged: continue
            a,b = merged[0]
            a = max(0, a-args.pad); b = b+args.pad
            strand = "+" if strand_vote[(q,s)]>=0 else "-"
            name = f"{q.lower()}_putative"
            sb.write(f"{s}\t{a}\t{b}\t{name}\t.\t{strand}\n")
            mb.write(f"{q}\t{name}\n")
            out_count+=1
        # print small summary to stderr
        print(f"[tblastn_to_spans] wrote {out_count} putative spans", file=sys.stderr)

if __name__ == "__main__":
    main()
