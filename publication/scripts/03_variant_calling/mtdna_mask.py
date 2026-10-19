#!/usr/bin/env python3
"""Build coordinate-preserving mt masks from full pileup and caller candidates."""
import argparse
import json
from collections import Counter


def records(path):
    with open(path) as handle:
        for line in handle:
            if line.strip():
                yield line.rstrip("\n").split("\t")


def build(args):
    lengths = {r[0]: int(r[1]) for r in records(args.fai)}
    n = lengths[args.contig]
    if not 0.5 < args.major_af <= 1 or not 0 < args.min_dp <= args.max_dp or args.indel_pad < 0:
        raise ValueError("Invalid depth, majority fraction, or indel padding")
    # Every reference position starts unknown, including sites absent from pileup.
    reasons = [set(["no_usable_coverage"]) for _ in range(n)]
    passed = set()
    for chrom, pos, ref, alt, gt in records(args.passed):
        if chrom != args.contig:
            continue
        if gt != "1" or len(ref) != 1 or alt not in "ACGT" or len(alt) != 1:
            raise ValueError("Expected biallelic haploid ALT SNPs in passed VCF")
        passed.add((int(pos), ref.upper(), alt.upper()))
    seen = set()
    for chrom, pos, ref, alt, dp, ad in records(args.pileup):
        if chrom != args.contig:
            continue
        pos = int(pos)
        if not 1 <= pos <= n or pos in seen:
            raise ValueError("Invalid or duplicate pileup position")
        seen.add(pos)
        why = reasons[pos - 1]
        why.clear()
        depth = 0 if dp == "." else int(dp)
        alleles = [ref.upper()] + alt.upper().split(",")
        counts = [0 if x == "." else int(x) for x in ad.split(",")]
        if len(counts) != len(alleles):
            raise ValueError("AD and allele dimensions differ")
        # Unknown and symbolic alleles cannot establish a confident base.
        base_counts = Counter()
        for allele, count in zip(alleles, counts):
            if len(allele) == 1 and allele in "ACGT":
                base_counts[allele] += count
        total = sum(counts)
        if depth < args.min_dp or total < args.min_dp:
            why.add("low_depth")
        if depth > args.max_dp:
            why.add("high_depth")
        if not base_counts or not total:
            why.add("no_usable_coverage")
            continue
        major, count = base_counts.most_common(1)[0]
        if count / total < args.major_af:
            why.add("mixed_or_uncertain")
        if major != ref.upper() and (pos, ref.upper(), major) not in passed:
            why.add("unsupported_major_alt")
        if any(p == pos and a != major for p, r, a in passed):
            why.add("pass_call_disagrees_with_pileup")
    # Conservative: mask every ALT candidate from either caller that did not
    # survive exact intersection/quality filtering. Mask indels and flanks.
    for chrom, pos, ref, alt, gt in records(args.candidates):
        if chrom != args.contig:
            continue
        pos = int(pos)
        indel = len(ref) != len(alt)
        if (pos, ref.upper(), alt.upper()) in passed and not indel:
            continue
        pad = args.indel_pad if indel else 0
        start = max(0, pos - 1 - pad)
        end = min(n, pos - 1 + len(ref) + pad)
        for i in range(start, end):
            reasons[i].add("indel_or_flank" if indel else "rejected_candidate")
    with open(args.bed, "w") as out:
        start = None
        for i in range(n + 1):
            masked = i < n and bool(reasons[i])
            if masked and start is None:
                start = i
            if not masked and start is not None:
                out.write(f"{args.contig}\t{start}\t{i}\n")
                start = None
    counts = Counter(reason for site in reasons for reason in site)
    with open(args.qa, "w") as out:
        json.dump(dict(reference_bp=n, masked_bp=sum(bool(x) for x in reasons),
                       callable_bp=sum(not x for x in reasons),
                       reason_counts=dict(counts), min_dp=args.min_dp,
                       max_dp=args.max_dp, major_af=args.major_af,
                       indel_pad=args.indel_pad), out, indent=2)
        out.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fai", "contig", "pileup", "passed", "candidates", "bed", "qa"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--min-dp", type=int, default=10)
    parser.add_argument("--max-dp", type=int, default=5000)
    parser.add_argument("--major-af", type=float, default=0.8)
    parser.add_argument("--indel-pad", type=int, default=5)
    build(parser.parse_args())
