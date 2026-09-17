#!/usr/bin/env python3
"""Corrected pooled piN/piS for the 13 mitochondrial protein-coding genes."""
import argparse, glob, gzip, math, os, re, sys
from collections import defaultdict
import pysam

CODE = {
 "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
 "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"W","TGG":"W",
 "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
 "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
 "ATT":"I","ATC":"I","ATA":"M","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
 "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"*","AGG":"*",
 "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
 "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
BASES=("A","C","G","T"); COMP={"A":"T","T":"A","C":"G","G":"C","N":"N"}
COMPLEX={"ND1":"CI","ND2":"CI","ND3":"CI","ND4":"CI","ND4L":"CI","ND5":"CI","ND6":"CI",
         "CYTB":"CIII","COX1":"CIV","COX2":"CIV","COX3":"CIV","ATP6":"CV","ATP8":"CV"}

def attrs(s):
 d={}
 for item in s.rstrip(";").split(";"):
  if "=" in item:
   k,v=item.split("=",1); d[k]=v
 return d

def read_map(path):
 rows=[]
 with open(path) as h:
  hd=h.readline().rstrip().split("\t"); ix={v:i for i,v in enumerate(hd)}
  for line in h:
   f=line.rstrip().split("\t")
   if f[ix["include"]].lower()!="yes": continue
   rows.append((int(f[ix["sync_column"]]),f[ix["population"]],int(f[ix["n_individuals"]])))
 rows.sort()
 return rows

def models(gff, contig):
 out=defaultdict(lambda:{"strand":None,"blocks":[]})
 with open(gff) as h:
  for line in h:
   if line.startswith("#"): continue
   f=line.rstrip().split("\t")
   if len(f)!=9 or f[2]!="CDS": continue
   a=attrs(f[8]); gene=a.get("gene") or a.get("Name")
   if gene not in COMPLEX: continue
   out[gene]["strand"]=f[6]; out[gene]["blocks"].append((int(f[3]),int(f[4])))
 for gene in out: out[gene]["contig"]=contig
 return out

def positions(m):
 blocks=sorted(m["blocks"],reverse=m["strand"]=="-"); z=[]
 for s,e in blocks: z.extend(range(s,e+1) if m["strand"]=="+" else range(e,s-1,-1))
 return z

def aa(c): return CODE.get(c) if len(c)==3 and set(c)<=set(BASES) else None

def opp(c,p):
 ref=aa(c)
 if ref in (None,"*"): return None
 n=s=0.0
 for b in BASES:
  if b==c[p]: continue
  x=list(c); x[p]=b; alt=aa("".join(x))
  if alt in (None,"*"): continue
  if alt==ref: s+=1/3
  else: n+=1/3
 return n,s

def cnt(field,minus):
 try: x=[int(v) for v in field.split(":")[:4]]
 except: return None
 g={"A":x[0],"T":x[1],"C":x[2],"G":x[3]}
 return ({"A":g["T"],"C":g["G"],"G":g["C"],"T":g["A"]} if minus else {b:g[b] for b in BASES})

def div(counts,codon,p,corr):
 dp=sum(counts.values()); fr={b:counts[b]/dp for b in BASES}; pn=ps=0.0
 for i,b1 in enumerate(BASES):
  for b2 in BASES[i+1:]:
   if not fr[b1] or not fr[b2]: continue
   c1=list(codon); c2=list(codon); c1[p]=b1; c2[p]=b2
   a1,a2=aa("".join(c1)),aa("".join(c2))
   if a1 in (None,"*") or a2 in (None,"*"): continue
   v=2*fr[b1]*fr[b2]*corr
   if a1==a2: ps+=v
   else: pn+=v
 return pn,ps

def fmt(x): return "NA" if x is None or not math.isfinite(x) else f"{x:.10g}"

def main():
 p=argparse.ArgumentParser(); p.add_argument("--sync-dir",required=True); p.add_argument("--gff",required=True)
 p.add_argument("--reference",required=True); p.add_argument("--contig",required=True); p.add_argument("--column-map",required=True)
 p.add_argument("--output",required=True); p.add_argument("--min-coverage",type=int,default=10); a=p.parse_args()
 pops=read_map(a.column_map); mod=models(a.gff,a.contig); fa=pysam.FastaFile(a.reference); seq=fa.fetch(a.contig).upper()
 header=["gene_id","symbol","role","complex","population","piN","piS","piN_piS","piN_sum","piS_sum",
         "N_opportunities","S_opportunities","callable_sites","CDS_length","callable_fraction","ratio_status"]
 with open(a.output,"w") as out:
  out.write("\t".join(header)+"\n")
  for gene in sorted(mod):
   m=mod[gene]; pos=positions(m)
   if len(pos)%3: print(f"[warn non-triplet] {gene} {len(pos)}",file=sys.stderr)
   site={}
   for start in range(0,len(pos)-2,3):
    gp=pos[start:start+3]; bases=[seq[x-1] for x in gp]
    codon="".join(bases if m["strand"]=="+" else [COMP.get(x,"N") for x in bases])
    if aa(codon) in (None,"*"): continue
    for cp,x in enumerate(gp):
     o=opp(codon,cp)
     if o: site[x]=(codon,cp,o[0],o[1])
   candidates=glob.glob(os.path.join(a.sync_dir,gene+".sync"))+glob.glob(os.path.join(a.sync_dir,gene.lower()+".sync"))
   if not candidates: raise FileNotFoundError(f"No sync for {gene}")
   acc={name:[0.,0.,0.,0.,0] for _,name,_ in pops}
   with open(candidates[0]) as h:
    for line in h:
     f=line.rstrip().split("\t"); position=int(f[1]); rec=site.get(position)
     if not rec: continue
     codon,cp,no,so=rec
     for col,name,nind in pops:
      counts=cnt(f[2+col],m["strand"]=="-")
      if counts is None or sum(counts.values())<a.min_coverage: continue
      pn,ps=div(counts,codon,cp,(2*nind)/(2*nind-1)); z=acc[name]
      z[0]+=pn; z[1]+=ps; z[2]+=no; z[3]+=so; z[4]+=1
   for _,name,_ in pops:
    pn0,ps0,no,so,nsite=acc[name]; pn=pn0/no if no else None; ps=ps0/so if so else None
    if pn is not None and ps is not None and ps>0: ratio,status=pn/ps,"finite"
    elif pn==0 and ps==0: ratio,status=None,"piN0_piS0"
    elif pn is not None and ps==0 and pn>0: ratio,status=None,"piNpositive_piS0"
    else: ratio,status=None,"insufficient_data"
    frac=nsite/len(pos) if pos else 0
    out.write("\t".join(map(str,[gene,gene,"mtOXPHOS",COMPLEX[gene],name,fmt(pn),fmt(ps),fmt(ratio),
      fmt(pn0),fmt(ps0),fmt(no),fmt(so),nsite,len(pos),fmt(frac),status]))+"\n")
   print(f"[done] {gene} strand={m['strand']} CDS={len(pos)} sites={len(site)}",file=sys.stderr)
 print(f"[output] {a.output}",file=sys.stderr)

if __name__=="__main__": main()
