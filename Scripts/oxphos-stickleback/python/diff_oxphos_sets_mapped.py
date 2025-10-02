#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, re, sys
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

# ---- Human sets (same as你之前那份) ----
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
"CIV assembly factors": S("""CEP89, CMC1, CMC2, COA1, COA3, COA4, COA5, COA6, COA7, COA8, COX10, COX11, COX14, COX15, COX16, COX17, COX18, COX19, COX20, HIGD1A, PET100, PET117, PNKD, SCO1, SCO2, SMIM20, SURF1, TACO1, TIMM21, TMEM177"""),
"Complex V": S("""ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E, ATP5IF1, ATP5MC1, ATP5MC2, ATP5MC3, ATP5MD, ATP5ME, ATP5MF, ATP5MG, ATP5MPL, ATP5PB, ATP5PD, ATP5PF, ATP5PO, ATPAF1, ATPAF2, ATPSCKMT, DMAC2L, FMC1, MT-ATP6, MT-ATP8, TMEM70"""),
"CV subunits": S("""ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E, ATP5IF1, ATP5MC1, ATP5MC2, ATP5MC3, ATP5MD, ATP5ME, ATP5MF, ATP5MG, ATP5MPL, ATP5PB, ATP5PD, ATP5PF, ATP5PO, DMAC2L, MT-ATP6, MT-ATP8"""),
"CV assembly factors": S("""ATPAF1, ATPAF2, ATPSCKMT, FMC1, TMEM70"""),
}
ALL_HUMAN = set().union(*HUMAN.values())

def is_mt(g): return g.startswith("MT-")

def load_ortholog_map(tsv):
    """human -> fish aliases (upper), then invert to fish->set(human)"""
    if not tsv: return {}
    fish2human = defaultdict(set)
    with open(tsv) as f:
        first = f.readline()
        if first and ("\t" in first):
            cols = first.rstrip("\n").split("\t")
            if not (("human" in cols[0].lower()) or ("symbol" in cols[0].lower())):
                # first line is data
                h, s = cols[0].strip().upper(), cols[1].strip().upper() if len(cols) > 1 else ""
                if h and s:
                    fish2human[s].add(h)
        else:
            # empty file or no tabs
            return {}
        for line in f:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2: continue
            h, s = parts[0].strip().upper(), parts[1].strip().upper()
            if h and s:
                fish2human[s].add(h)
    return fish2human

def canon_name_for_human_match(x, keep_ab=False):
    """Return a set of human symbols we can match by 'same name' logic."""
    out = set()
    X = x.upper()
    if X in ALL_HUMAN:
        out.add(X)
    if (not keep_ab) and re.search(r"[AB]$", X) and len(X) > 1:
        base = X[:-1]
        if base in ALL_HUMAN:
            out.add(base)
    return out

def read_fish_genes_as_humans(span_bed, fish2human_map, keep_ab=False):
    """Read fish gene names and convert to 'present human symbols' via mapping or same-name match."""
    present = set()
    with open(span_bed) as f:
        for line in f:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4: continue
            fish = parts[3].strip().upper()
            # 1) mapping fish->human
            if fish in fish2human_map:
                present |= fish2human_map[fish]
            # 2) also尝试去掉尾部 a/b 后再查映射（谨慎）
            if (not keep_ab) and re.search(r"[AB]$", fish) and len(fish) > 1:
                base = fish[:-1]
                if base in fish2human_map:
                    present |= fish2human_map[base]
            # 3) same-name 匹配（fish 名字本身就是人类符号，或安全合并 a/b）
            present |= canon_name_for_human_match(fish, keep_ab=keep_ab)
    return present

def main():
    args = parse_args()
    # 可选忽略 MT-：先复制一份集合
    human_sets = {k:set(v) for k,v in HUMAN.items()}
    if args.ignore_mt:
        for k in human_sets:
            human_sets[k] = {g for g in human_sets[k] if not is_mt(g)}
    fish2human = load_ortholog_map(args.ortholog_map) if args.ortholog_map else {}
    present_humans = read_fish_genes_as_humans(args.span_bed, fish2human, keep_ab=args.keep_ab)

    # 输出
    summary = args.out_prefix + ".summary.tsv"
    with open(summary, "w") as fo:
        fo.write("set\thuman_n\tfish_present_n\tmissing_n\n")
        for set_name, H in human_sets.items():
            present = sorted(H & present_humans)
            missing = sorted(H - present_humans)
            fo.write(f"{set_name}\t{len(H)}\t{len(present)}\t{len(missing)}\n")
            # per-set missing list
            tag = re.sub(r"[^A-Za-z0-9]+","_", set_name)
            with open(args.out_prefix + f".{tag}.missing.txt", "w") as fm:
                fm.write("\n".join(missing) + ("\n" if missing else ""))

    print("✅ Wrote:", summary)
    print("   + per-set missing lists alongside the summary.")

if __name__ == "__main__":
    main()
