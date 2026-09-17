#!/usr/bin/env python3
"""Extract gene-level M0 dN, dS, and omega from the corrected mt13 mlc files."""
import csv,re
from pathlib import Path
ROOT=Path("/mnt/spareHD_2/genomewide_codeml_kuster/10_erc_mt13_codeml/results")
OUT=ROOT.parent/"mt13_M0_updated.tsv"
GENES=["ATP6","ATP8","COX1","COX2","COX3","CYTB","ND1","ND2","ND3","ND4","ND4L","ND5","ND6"]
CX={**{g:"CI" for g in ["ND1","ND2","ND3","ND4","ND4L","ND5","ND6"]},"CYTB":"CIII","COX1":"CIV","COX2":"CIV","COX3":"CIV","ATP6":"CV","ATP8":"CV"}
def val(pattern,text,gene,label):
    m=re.search(pattern,text,re.M)
    if not m: raise ValueError(f"{gene}: missing {label}")
    return float(m.group(1))
rows=[]
for g in GENES:
    p=ROOT/g/"mlc"; t=p.read_text(errors="replace")
    rows.append(dict(gene=g,gene_id=g,symbol=g,role="mt",complex=CX[g],model="M0",
      dN=val(r"tree length for dN:\s*([-+0-9.eE]+)",t,g,"dN"),
      dS=val(r"tree length for dS:\s*([-+0-9.eE]+)",t,g,"dS"),
      omega=val(r"omega\s+\(dN/dS\)\s*=\s*([-+0-9.eE]+)",t,g,"omega"),mlc=str(p)))
with OUT.open("w",newline="") as h:
    w=csv.DictWriter(h,delimiter="\t",fieldnames=rows[0].keys(),lineterminator="\n");w.writeheader();w.writerows(rows)
print(f"[DONE] Wrote 13 genes: {OUT}")
